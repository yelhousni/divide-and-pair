package curve25519

import (
	"math/big"
	"sync"

	fp "github.com/yelhousni/divide-and-pair/curve25519/fp"
)

var (
	subgroupInitOnce sync.Once
	aPornin          fp.Element // A = 2(a+d)
	bPornin          fp.Element // B = (a-d)^2
	aPrimePorn       fp.Element // A' = -2*A
	bPrimePorn       fp.Element // B' = A^2 - 4*B
	aMinusD          fp.Element // a - d
	sqrtMinOne       fp.Element // i = sqrt(-1) mod q
	twoI             fp.Element // 2i
	sqrt2Bp          fp.Element // sqrt(2*B')
	subgroupOrder    big.Int
)

func initSubgroupConstants() {
	initOnce.Do(initCurveParams)

	aMinusD.Sub(&curveParams.A, &curveParams.D)

	var aPlusD fp.Element
	aPlusD.Add(&curveParams.A, &curveParams.D)
	aPornin.Double(&aPlusD)

	bPornin.Square(&aMinusD)

	aPrimePorn.Double(&aPornin)
	aPrimePorn.Neg(&aPrimePorn)

	var a2, b4 fp.Element
	a2.Square(&aPornin)
	b4.Double(&bPornin)
	b4.Double(&b4)
	bPrimePorn.Sub(&a2, &b4)

	sqrtMinOne.SetString("19681161376707505956807079304988542015446066515923890162744021073123829784752")
	twoI.Double(&sqrtMinOne)

	// sqrt(2*B') — 2 is NQR for this field, but 2*B' is QR
	var twoBp fp.Element
	twoBp.Double(&bPrimePorn)
	sqrt2Bp.Sqrt(&twoBp)

	subgroupOrder.Set(&curveParams.Order)
}

// isLowOrder checks if an affine point (X, Y) is a low-order point.
// For curve25519 (cofactor 8), low-order points in E[8] are exactly
// the points with X=0, Y=0, or X = ±i*Y.
func isLowOrder(X, Y *fp.Element) bool {
	if X.IsZero() || Y.IsZero() {
		return true
	}
	var iY fp.Element
	iY.Mul(&sqrtMinOne, Y)
	if X.Equal(&iY) {
		return true
	}
	iY.Neg(&iY)
	return X.Equal(&iY)
}

// edwardsToPorninMontgomery converts an affine Edwards point (x, y) to
// Pornin's Montgomery (u, w) coordinates, batching the two inversions
// into one using Montgomery's trick: 1/(a*b) then recover 1/a and 1/b.
func edwardsToPorninMontgomery(p *PointAffine) (u, w fp.Element) {
	// u = (a-d)(1+y)/(1-y), w = 2/x
	// Need 1/(1-y) and 1/x → batch: inv = 1/(x*(1-y))
	var oneMinusY, prod, inv fp.Element
	var one fp.Element
	one.SetOne()
	oneMinusY.Sub(&one, &p.Y)
	prod.Mul(&p.X, &oneMinusY) // x * (1-y)
	inv.Inverse(&prod)         // 1 / (x * (1-y))

	// 1/(1-y) = inv * x,  1/x = inv * (1-y)
	var invOneMinusY, invX fp.Element
	invOneMinusY.Mul(&inv, &p.X)
	invX.Mul(&inv, &oneMinusY)

	var onePlusY fp.Element
	onePlusY.Add(&one, &p.Y)
	u.Mul(&aMinusD, &onePlusY)
	u.Mul(&u, &invOneMinusY)

	w.Double(&invX) // 2/x

	return
}

// sqrtFn computes dst = sqrt(src), returning nil if src is not a QR.
type sqrtFn func(dst, src *fp.Element) *fp.Element

func sqrtNative(dst, src *fp.Element) *fp.Element  { return dst.Sqrt(src) }
func sqrtFilippo(dst, src *fp.Element) *fp.Element { return dst.SqrtFilippo(src) }

// halvePornin performs one division-free halving on Pornin's Montgomery
// coordinates (u, w) with scaling factor e. Returns false if not halvable.
//
// Following crrl, the halvings are performed division-free by tracking
// a scaling factor e such that we work on the isomorphic curve
// Curve(A*e², B*e⁴). When up is not a QR in inverse psi1, we switch
// to an isomorphic curve instead of computing B'/up.
func halvePornin(u, w, e *fp.Element, sqrt sqrtFn) bool {
	// Inverse iso: (u, w) -> (4u, 2w)
	var us, ws fp.Element
	us.Double(u)
	us.Double(&us)
	ws.Double(w)

	// Inverse psi2: wp = sqrt(us)
	var wp fp.Element
	if sqrt(&wp, &us) == nil {
		return false
	}

	// up = (us - A'*e² - wp*ws) / 2
	var up, tmp, e2 fp.Element
	e2.Square(e)
	tmp.Mul(&e2, &aPrimePorn)
	up.Sub(&us, &tmp)
	tmp.Mul(&wp, &ws)
	up.Sub(&up, &tmp)
	up.Halve()

	// Inverse psi1.
	var sqrtUp fp.Element
	if sqrt(&sqrtUp, &up) != nil {
		// QR case
		w.Set(&sqrtUp)
		tmp.Mul(&e2, &aPornin)
		u.Square(w)
		u.Sub(u, &tmp)
		tmp.Mul(w, &wp)
		u.Sub(u, &tmp)
		u.Halve()
	} else {
		// NQR case: compute tt = sqrt(2*up)
		var twoUp, tt fp.Element
		twoUp.Double(&up)
		sqrt(&tt, &twoUp)

		wp.Mul(&wp, &tt)
		w.Mul(&sqrt2Bp, &e2)
		e.Mul(e, &tt)

		e2.Square(e)
		tmp.Mul(&e2, &aPornin)
		u.Square(w)
		u.Sub(u, &tmp)
		tmp.Mul(w, &wp)
		u.Sub(u, &tmp)
		u.Halve()

		w.Neg(w)
	}
	return true
}

// quarticCriterion computes the quartic test: χ₄(uR * (wR + 2i)²) using the given symbol function.
func quarticCriterion(uR, wR *fp.Element, symbolFn func(*fp.Element) uint8) bool {
	var wShift, f fp.Element
	wShift.Add(wR, &twoI)
	f.Square(&wShift)
	f.Mul(&f, uR)
	return symbolFn(&f) == 0
}

// quarticCriterionScaled computes χ₄(u * (w + 2i·e)²) using the given symbol function.
func quarticCriterionScaled(u, w, e *fp.Element, symbolFn func(*fp.Element) uint8) bool {
	var twoIe, wShift, f fp.Element
	twoIe.Mul(&twoI, e)
	wShift.Add(w, &twoIe)
	f.Square(&wShift)
	f.Mul(&f, u)
	return symbolFn(&f) == 0
}

// isInSubGroupGeneric is the shared driver for all halving-based subgroup tests.
func isInSubGroupGeneric(p *PointAffine, halvings int, sqrt sqrtFn, finalCheck func(u, w, e *fp.Element) bool) bool {
	subgroupInitOnce.Do(initSubgroupConstants)

	if isLowOrder(&p.X, &p.Y) {
		return p.IsZero()
	}

	u, w := edwardsToPorninMontgomery(p)
	var e fp.Element
	e.SetOne()

	for range halvings {
		if !halvePornin(&u, &w, &e, sqrt) {
			return false
		}
	}

	return finalCheck(&u, &w, &e)
}

// isInSubGroupNaive tests subgroup membership by scalar multiplication by ℓ.
func (p *PointAffine) isInSubGroupNaive() bool {
	subgroupInitOnce.Do(initSubgroupConstants)
	var res PointAffine
	res.ScalarMultiplication(p, &subgroupOrder)
	return res.IsZero()
}

// isInSubGroupPornin tests subgroup membership using Pornin's method:
// 2 halvings + 1 Legendre symbol (3rd halving check).
func (p *PointAffine) isInSubGroupPornin() bool {
	return isInSubGroupGeneric(p, 2, sqrtNative, func(u, w, e *fp.Element) bool {
		return u.Legendre() == 1
	})
}

// isInSubGroupQuartic tests subgroup membership using our improved method:
// 1 halving + 1 quartic symbol (Weilert GCD).
func (p *PointAffine) isInSubGroupQuartic() bool {
	return isInSubGroupGeneric(p, 1, sqrtNative, func(u, w, e *fp.Element) bool {
		return quarticCriterionScaled(u, w, e, (*fp.Element).QuarticSymbol)
	})
}

// isInSubGroupQuarticExp tests subgroup membership using:
// 1 halving + 1 quartic symbol (addition chain exponentiation).
func (p *PointAffine) isInSubGroupQuarticExp() bool {
	return isInSubGroupGeneric(p, 1, sqrtNative, func(u, w, e *fp.Element) bool {
		return quarticCriterionScaled(u, w, e, (*fp.Element).QuarticSymbolExp)
	})
}

// isInSubGroupPorninFilippo is isInSubGroupPornin but uses filippo's
// 5×51-bit field for sqrt and Legendre (faster squaring).
func (p *PointAffine) isInSubGroupPorninFilippo() bool {
	return isInSubGroupGeneric(p, 2, sqrtFilippo, func(u, w, e *fp.Element) bool {
		return u.LegendreFilippo() == 1
	})
}

// isInSubGroupQuarticExpFilippo is isInSubGroupQuarticExp but uses filippo's
// field for sqrt and quartic symbol.
func (p *PointAffine) isInSubGroupQuarticExpFilippo() bool {
	return isInSubGroupGeneric(p, 1, sqrtFilippo, func(u, w, e *fp.Element) bool {
		return quarticCriterionScaled(u, w, e, (*fp.Element).QuarticSymbolExpFilippo)
	})
}

// isInSubGroupQuarticFilippo is isInSubGroupQuartic but uses filippo's
// field for sqrt and quartic symbol (Weilert GCD with filippo fallback).
func (p *PointAffine) isInSubGroupQuarticFilippo() bool {
	return isInSubGroupGeneric(p, 1, sqrtFilippo, func(u, w, e *fp.Element) bool {
		return quarticCriterionScaled(u, w, e, (*fp.Element).QuarticSymbolFilippo)
	})
}

// IsInSubGroup tests subgroup membership using the fastest available method.
func (p *PointAffine) IsInSubGroup() bool {
	return p.isInSubGroupQuartic()
}
