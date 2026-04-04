package curve25519

import (
	"math/big"
	"sync"

	fp "github.com/yelhousni/divide-and-pair/curve25519/fp"
)

var (
	subgroupInitOnce sync.Once
	aPornin          fp.Element // A_p = 2(a+d)
	bPornin          fp.Element // B_p = (a-d)^2
	aPrimePorn       fp.Element // A' = -2*A_p
	bPrimePorn       fp.Element // B' = A_p^2 - 4*B_p
	aMinusD          fp.Element // a - d
	sqrtMinOne       fp.Element // i = sqrt(-1) mod q
	twoI             fp.Element // 2i
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

	subgroupOrder.Set(&curveParams.Order)
}

// isLowOrder checks if an affine point (X, Y) is a low-order point.
// Low-order iff X=0, Y=0, or X = ±i*Y.
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
	prod.Mul(&p.X, &oneMinusY)   // x * (1-y)
	inv.Inverse(&prod)            // 1 / (x * (1-y))

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

// halvePornin performs one halving on Pornin's Montgomery curve C(A_p, B_p).
func halvePornin(u, w *fp.Element) (uOut, wOut fp.Element, ok bool) {
	var upp, wpp fp.Element
	upp.Double(u)
	upp.Double(&upp) // 4u
	wpp.Double(w)    // 2w

	// w' = sqrt(u'')
	var wp fp.Element
	if wp.Sqrt(&upp) == nil {
		return fp.Element{}, fp.Element{}, false
	}

	// u' = (w'^2 - A' - w'*w'') / 2
	var up, tmp fp.Element
	up.Square(&wp)
	up.Sub(&up, &aPrimePorn)
	tmp.Mul(&wp, &wpp)
	up.Sub(&up, &tmp)
	up.Halve()

	// Check if u' is QR; if not, switch to other preimage
	var sqrtUp fp.Element
	if sqrtUp.Sqrt(&up) == nil {
		var bDivU fp.Element
		bDivU.Inverse(&up)
		bDivU.Mul(&bDivU, &bPrimePorn)
		up.Set(&bDivU)
		wp.Neg(&wp)
		if sqrtUp.Sqrt(&up) == nil {
			return fp.Element{}, fp.Element{}, false
		}
	}

	wOut.Set(&sqrtUp)
	uOut.Square(&wOut)
	uOut.Sub(&uOut, &aPornin)
	tmp.Mul(&wOut, &wp)
	uOut.Sub(&uOut, &tmp)
	uOut.Halve()

	return uOut, wOut, true
}

// quarticCriterion computes the quartic test: χ₄(uR * (wR + 2i)²) using the given symbol function.
func quarticCriterion(uR, wR *fp.Element, symbolFn func(*fp.Element) uint8) bool {
	var wShift, f fp.Element
	wShift.Add(wR, &twoI)
	f.Square(&wShift)
	f.Mul(&f, uR)
	return symbolFn(&f) == 0
}

// IsInSubGroupNaive tests subgroup membership by scalar multiplication by ℓ.
func (p *PointAffine) IsInSubGroupNaive() bool {
	subgroupInitOnce.Do(initSubgroupConstants)
	var res PointAffine
	res.ScalarMultiplication(p, &subgroupOrder)
	return res.IsZero()
}

// IsInSubGroupPornin tests subgroup membership using Pornin's method:
// 2 halvings + 1 Legendre symbol.
func (p *PointAffine) IsInSubGroupPornin() bool {
	subgroupInitOnce.Do(initSubgroupConstants)

	if isLowOrder(&p.X, &p.Y) {
		return p.IsZero()
	}

	u, w := edwardsToPorninMontgomery(p)

	u1, w1, ok := halvePornin(&u, &w)
	if !ok {
		return false
	}
	u2, _, ok := halvePornin(&u1, &w1)
	if !ok {
		return false
	}

	return u2.Legendre() == 1
}

// IsInSubGroupQuartic tests subgroup membership using our improved method:
// 1 halving + 1 quartic symbol (Weilert GCD).
func (p *PointAffine) IsInSubGroupQuartic() bool {
	subgroupInitOnce.Do(initSubgroupConstants)

	if isLowOrder(&p.X, &p.Y) {
		return p.IsZero()
	}

	u, w := edwardsToPorninMontgomery(p)

	uR, wR, ok := halvePornin(&u, &w)
	if !ok {
		return false
	}

	return quarticCriterion(&uR, &wR, (*fp.Element).QuarticSymbol)
}

// IsInSubGroupQuarticExp tests subgroup membership using:
// 1 halving + 1 quartic symbol (addition chain exponentiation).
func (p *PointAffine) IsInSubGroupQuarticExp() bool {
	subgroupInitOnce.Do(initSubgroupConstants)

	if isLowOrder(&p.X, &p.Y) {
		return p.IsZero()
	}

	u, w := edwardsToPorninMontgomery(p)

	uR, wR, ok := halvePornin(&u, &w)
	if !ok {
		return false
	}

	return quarticCriterion(&uR, &wR, (*fp.Element).QuarticSymbolExp)
}
