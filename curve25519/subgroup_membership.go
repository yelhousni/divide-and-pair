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

// quarticCriterion computes the quartic test: χ₄(uR * (wR + 2i)²) using the given symbol function.
func quarticCriterion(uR, wR *fp.Element, symbolFn func(*fp.Element) uint8) bool {
	var wShift, f fp.Element
	wShift.Add(wR, &twoI)
	f.Square(&wShift)
	f.Mul(&f, uR)
	return symbolFn(&f) == 0
}

// isInSubGroupNaive tests subgroup membership by scalar multiplication by ℓ.
func (p *PointAffine) isInSubGroupNaive() bool {
	subgroupInitOnce.Do(initSubgroupConstants)
	var res PointAffine
	res.ScalarMultiplication(p, &subgroupOrder)
	return res.IsZero()
}

// isInSubGroupPornin tests subgroup membership using Pornin's method
// (https://eprint.iacr.org/2022/1164):
// 2 halvings + 1 Legendre symbol (3rd halving check).
//
// Following crrl, the halvings are performed division-free by tracking
// a scaling factor e such that we work on the isomorphic curve
// Curve(A*e², B*e⁴). When up is not a QR in inverse psi1, we switch
// to an isomorphic curve instead of computing B'/up.
func (p *PointAffine) isInSubGroupPornin() bool {
	subgroupInitOnce.Do(initSubgroupConstants)

	if isLowOrder(&p.X, &p.Y) {
		return p.IsZero()
	}

	u, w := edwardsToPorninMontgomery(p)

	// Scaling factor e (starts at 1 since we converted to affine above).
	var e fp.Element
	e.SetOne()

	for range 2 {
		// Inverse iso: (u, w) -> (4u, 2w) on Curve(As*e², Bs*e⁴)
		var us, ws fp.Element
		us.Double(&u)
		us.Double(&us) // 4u
		ws.Double(&w)  // 2w

		// Inverse psi2: wp = sqrt(us)
		// If us is not a square, the point cannot be halved.
		var wp fp.Element
		if wp.Sqrt(&us) == nil {
			return false
		}
		// up = (us - A'*e² - wp*ws) / 2
		var up, tmp, e2 fp.Element
		e2.Square(&e)
		tmp.Mul(&e2, &aPrimePorn)
		up.Sub(&us, &tmp)
		tmp.Mul(&wp, &ws)
		up.Sub(&up, &tmp)
		up.Halve()

		// Inverse psi1.
		// If up is a QR:
		//   w = sqrt(up)
		//   u = (w² - A*e² - w*wp) / 2
		// If up is NQR (2 is NQR for this field, so 2*up is QR):
		//   tt = sqrt(2*up)
		//   Switch to isomorphic curve:
		//     up' = 2*up²,  wp' = wp*tt,  e' = e*tt
		//     w = -sqrt(2*B')*e²
		//   Then: u = (w² - A*e'² - w*wp') / 2; w = -w
		var sqrtUp fp.Element
		if sqrtUp.Sqrt(&up) != nil {
			// QR case
			w.Set(&sqrtUp)
			// u = (w² - A*e² - w*wp) / 2
			tmp.Mul(&e2, &aPornin)
			u.Square(&w)
			u.Sub(&u, &tmp)
			tmp.Mul(&w, &wp)
			u.Sub(&u, &tmp)
			u.Halve()
		} else {
			// NQR case: compute tt = sqrt(2*up)
			var twoUp fp.Element
			twoUp.Double(&up)
			var tt fp.Element
			tt.Sqrt(&twoUp)

			// Update isomorphism
			wp.Mul(&wp, &tt)
			w.Mul(&sqrt2Bp, &e2)
			e.Mul(&e, &tt)

			// u = (w² - A*e² - w*wp) / 2
			e2.Square(&e)
			tmp.Mul(&e2, &aPornin)
			u.Square(&w)
			u.Sub(&u, &tmp)
			tmp.Mul(&w, &wp)
			u.Sub(&u, &tmp)
			u.Halve()

			// Negate w for next iteration (sign from using Bp/up preimage)
			w.Neg(&w)
		}
	}

	// Third halving: only check whether u is a QR (Legendre symbol).
	return u.Legendre() == 1
}

// isInSubGroupQuartic tests subgroup membership using our improved method:
// 1 halving + 1 quartic symbol (Weilert GCD).
func (p *PointAffine) isInSubGroupQuartic() bool {
	subgroupInitOnce.Do(initSubgroupConstants)

	if isLowOrder(&p.X, &p.Y) {
		return p.IsZero()
	}

	u, w := edwardsToPorninMontgomery(p)

	// Single halving (division-free with scaling factor e).
	var e fp.Element
	e.SetOne()

	// Inverse iso
	var us, ws fp.Element
	us.Double(&u)
	us.Double(&us)
	ws.Double(&w)

	// Inverse psi2
	var wp fp.Element
	if wp.Sqrt(&us) == nil {
		return false
	}
	var up, tmp, e2 fp.Element
	e2.Square(&e)
	tmp.Mul(&e2, &aPrimePorn)
	up.Sub(&us, &tmp)
	tmp.Mul(&wp, &ws)
	up.Sub(&up, &tmp)
	up.Halve()

	// Inverse psi1
	var sqrtUp fp.Element
	if sqrtUp.Sqrt(&up) != nil {
		w.Set(&sqrtUp)
		tmp.Mul(&e2, &aPornin)
		u.Square(&w)
		u.Sub(&u, &tmp)
		tmp.Mul(&w, &wp)
		u.Sub(&u, &tmp)
		u.Halve()
	} else {
		var twoUp fp.Element
		twoUp.Double(&up)
		var tt fp.Element
		tt.Sqrt(&twoUp)

		wp.Mul(&wp, &tt)
		w.Mul(&sqrt2Bp, &e2)
		e.Mul(&e, &tt)

		e2.Square(&e)
		tmp.Mul(&e2, &aPornin)
		u.Square(&w)
		u.Sub(&u, &tmp)
		tmp.Mul(&w, &wp)
		u.Sub(&u, &tmp)
		u.Halve()

		w.Neg(&w)
	}

	// Quartic criterion on (u, w) — note: u and w are on an isomorphic
	// curve scaled by e. The quartic symbol is scale-invariant for the
	// criterion χ₄(u * (w + 2i*e)²) since the e factors cancel in the
	// exponentiation. But we need to adjust 2i by e.
	var wShift, f fp.Element
	var twoIe fp.Element
	twoIe.Mul(&twoI, &e)
	wShift.Add(&w, &twoIe)
	f.Square(&wShift)
	f.Mul(&f, &u)
	return f.QuarticSymbol() == 0
}

// isInSubGroupQuarticExp tests subgroup membership using:
// 1 halving + 1 quartic symbol (addition chain exponentiation).
func (p *PointAffine) isInSubGroupQuarticExp() bool {
	subgroupInitOnce.Do(initSubgroupConstants)

	if isLowOrder(&p.X, &p.Y) {
		return p.IsZero()
	}

	u, w := edwardsToPorninMontgomery(p)

	// Single halving (division-free with scaling factor e).
	var e fp.Element
	e.SetOne()

	// Inverse iso
	var us, ws fp.Element
	us.Double(&u)
	us.Double(&us)
	ws.Double(&w)

	// Inverse psi2
	var wp fp.Element
	if wp.Sqrt(&us) == nil {
		return false
	}
	var up, tmp, e2 fp.Element
	e2.Square(&e)
	tmp.Mul(&e2, &aPrimePorn)
	up.Sub(&us, &tmp)
	tmp.Mul(&wp, &ws)
	up.Sub(&up, &tmp)
	up.Halve()

	// Inverse psi1
	var sqrtUp fp.Element
	if sqrtUp.Sqrt(&up) != nil {
		w.Set(&sqrtUp)
		tmp.Mul(&e2, &aPornin)
		u.Square(&w)
		u.Sub(&u, &tmp)
		tmp.Mul(&w, &wp)
		u.Sub(&u, &tmp)
		u.Halve()
	} else {
		var twoUp fp.Element
		twoUp.Double(&up)
		var tt fp.Element
		tt.Sqrt(&twoUp)

		wp.Mul(&wp, &tt)
		w.Mul(&sqrt2Bp, &e2)
		e.Mul(&e, &tt)

		e2.Square(&e)
		tmp.Mul(&e2, &aPornin)
		u.Square(&w)
		u.Sub(&u, &tmp)
		tmp.Mul(&w, &wp)
		u.Sub(&u, &tmp)
		u.Halve()

		w.Neg(&w)
	}

	var wShift, f fp.Element
	var twoIe fp.Element
	twoIe.Mul(&twoI, &e)
	wShift.Add(&w, &twoIe)
	f.Square(&wShift)
	f.Mul(&f, &u)
	return f.QuarticSymbolExp() == 0
}

// isInSubGroupPorninFilippo is isInSubGroupPornin but uses filippo's
// 5×51-bit field for sqrt and Legendre (faster squaring).
func (p *PointAffine) isInSubGroupPorninFilippo() bool {
	subgroupInitOnce.Do(initSubgroupConstants)

	if isLowOrder(&p.X, &p.Y) {
		return p.IsZero()
	}

	u, w := edwardsToPorninMontgomery(p)

	var e fp.Element
	e.SetOne()

	for range 2 {
		var us, ws fp.Element
		us.Double(&u)
		us.Double(&us)
		ws.Double(&w)

		var wp fp.Element
		if wp.SqrtFilippo(&us) == nil {
			return false
		}
		var up, tmp, e2 fp.Element
		e2.Square(&e)
		tmp.Mul(&e2, &aPrimePorn)
		up.Sub(&us, &tmp)
		tmp.Mul(&wp, &ws)
		up.Sub(&up, &tmp)
		up.Halve()

		var sqrtUp fp.Element
		if sqrtUp.SqrtFilippo(&up) != nil {
			w.Set(&sqrtUp)
			tmp.Mul(&e2, &aPornin)
			u.Square(&w)
			u.Sub(&u, &tmp)
			tmp.Mul(&w, &wp)
			u.Sub(&u, &tmp)
			u.Halve()
		} else {
			var twoUp fp.Element
			twoUp.Double(&up)
			var tt fp.Element
			tt.SqrtFilippo(&twoUp)

			wp.Mul(&wp, &tt)
			w.Mul(&sqrt2Bp, &e2)
			e.Mul(&e, &tt)

			e2.Square(&e)
			tmp.Mul(&e2, &aPornin)
			u.Square(&w)
			u.Sub(&u, &tmp)
			tmp.Mul(&w, &wp)
			u.Sub(&u, &tmp)
			u.Halve()

			w.Neg(&w)
		}
	}

	return u.LegendreFilippo() == 1
}

// isInSubGroupQuarticExpFilippo is isInSubGroupQuarticExp but uses filippo's
// field for sqrt and quartic symbol.
func (p *PointAffine) isInSubGroupQuarticExpFilippo() bool {
	subgroupInitOnce.Do(initSubgroupConstants)

	if isLowOrder(&p.X, &p.Y) {
		return p.IsZero()
	}

	u, w := edwardsToPorninMontgomery(p)

	var e fp.Element
	e.SetOne()

	var us, ws fp.Element
	us.Double(&u)
	us.Double(&us)
	ws.Double(&w)

	var wp fp.Element
	if wp.SqrtFilippo(&us) == nil {
		return false
	}
	var up, tmp, e2 fp.Element
	e2.Square(&e)
	tmp.Mul(&e2, &aPrimePorn)
	up.Sub(&us, &tmp)
	tmp.Mul(&wp, &ws)
	up.Sub(&up, &tmp)
	up.Halve()

	var sqrtUp fp.Element
	if sqrtUp.SqrtFilippo(&up) != nil {
		w.Set(&sqrtUp)
		tmp.Mul(&e2, &aPornin)
		u.Square(&w)
		u.Sub(&u, &tmp)
		tmp.Mul(&w, &wp)
		u.Sub(&u, &tmp)
		u.Halve()
	} else {
		var twoUp fp.Element
		twoUp.Double(&up)
		var tt fp.Element
		tt.SqrtFilippo(&twoUp)

		wp.Mul(&wp, &tt)
		w.Mul(&sqrt2Bp, &e2)
		e.Mul(&e, &tt)

		e2.Square(&e)
		tmp.Mul(&e2, &aPornin)
		u.Square(&w)
		u.Sub(&u, &tmp)
		tmp.Mul(&w, &wp)
		u.Sub(&u, &tmp)
		u.Halve()

		w.Neg(&w)
	}

	var wShift, f fp.Element
	var twoIe fp.Element
	twoIe.Mul(&twoI, &e)
	wShift.Add(&w, &twoIe)
	f.Square(&wShift)
	f.Mul(&f, &u)
	return f.QuarticSymbolExpFilippo() == 0
}

// isInSubGroupQuarticFilippo is isInSubGroupQuartic but uses filippo's
// field for sqrt and quartic symbol (Weilert GCD with filippo fallback).
func (p *PointAffine) isInSubGroupQuarticFilippo() bool {
	subgroupInitOnce.Do(initSubgroupConstants)

	if isLowOrder(&p.X, &p.Y) {
		return p.IsZero()
	}

	u, w := edwardsToPorninMontgomery(p)

	var e fp.Element
	e.SetOne()

	var us, ws fp.Element
	us.Double(&u)
	us.Double(&us)
	ws.Double(&w)

	var wp fp.Element
	if wp.SqrtFilippo(&us) == nil {
		return false
	}
	var up, tmp, e2 fp.Element
	e2.Square(&e)
	tmp.Mul(&e2, &aPrimePorn)
	up.Sub(&us, &tmp)
	tmp.Mul(&wp, &ws)
	up.Sub(&up, &tmp)
	up.Halve()

	var sqrtUp fp.Element
	if sqrtUp.SqrtFilippo(&up) != nil {
		w.Set(&sqrtUp)
		tmp.Mul(&e2, &aPornin)
		u.Square(&w)
		u.Sub(&u, &tmp)
		tmp.Mul(&w, &wp)
		u.Sub(&u, &tmp)
		u.Halve()
	} else {
		var twoUp fp.Element
		twoUp.Double(&up)
		var tt fp.Element
		tt.SqrtFilippo(&twoUp)

		wp.Mul(&wp, &tt)
		w.Mul(&sqrt2Bp, &e2)
		e.Mul(&e, &tt)

		e2.Square(&e)
		tmp.Mul(&e2, &aPornin)
		u.Square(&w)
		u.Sub(&u, &tmp)
		tmp.Mul(&w, &wp)
		u.Sub(&u, &tmp)
		u.Halve()

		w.Neg(&w)
	}

	var wShift, f fp.Element
	var twoIe fp.Element
	twoIe.Mul(&twoI, &e)
	wShift.Add(&w, &twoIe)
	f.Square(&wShift)
	f.Mul(&f, &u)
	return f.QuarticSymbolFilippo() == 0
}

// IsInSubGroup tests subgroup membership using the fastest available method.
func (p *PointAffine) IsInSubGroup() bool {
	return p.isInSubGroupQuartic()
}
