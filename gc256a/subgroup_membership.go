package gc256a

import (
	"math/big"
	"sync"

	fp "github.com/yelhousni/divide-and-pair/gc256a/fp"
)

var (
	subgroupInitOnce sync.Once
	aMinusD          fp.Element // a - d
	aPornin          fp.Element // A = 2(a+d)
	bPornin          fp.Element // B = (a-d)^2
	aPrimePorn       fp.Element // A' = -2*A
	bPrimePorn       fp.Element // B' = A^2 - 4*B
	minusAPornin     fp.Element // -A
	sqrtMinBp        fp.Element // sqrt(-B')
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

	minusAPornin.Neg(&aPornin)

	var negBp fp.Element
	negBp.Neg(&bPrimePorn)
	sqrtMinBp.Sqrt(&negBp)

	subgroupOrder.Set(&curveParams.Order)
}

// isLowOrder checks if an affine point (X, Y) is a low-order point.
// For GC256A (a=1, cofactor 4), the low-order points are exactly
// the points with X=0 or Y=0.
func isLowOrder(X, Y *fp.Element) bool {
	return X.IsZero() || Y.IsZero()
}

// edwardsToPorninMontgomery converts an affine Edwards point (x, y) to
// Pornin's Montgomery (u, w) coordinates.
func edwardsToPorninMontgomery(p *PointAffine) (u, w fp.Element) {
	var oneMinusY, prod, inv fp.Element
	var one fp.Element
	one.SetOne()
	oneMinusY.Sub(&one, &p.Y)
	prod.Mul(&p.X, &oneMinusY)
	inv.Inverse(&prod)

	var invOneMinusY, invX fp.Element
	invOneMinusY.Mul(&inv, &p.X)
	invX.Mul(&inv, &oneMinusY)

	var onePlusY fp.Element
	onePlusY.Add(&one, &p.Y)
	u.Mul(&aMinusD, &onePlusY)
	u.Mul(&u, &invOneMinusY)

	w.Double(&invX)

	return
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
// 1 halving + 1 Legendre symbol, division-free with scaling factor e.
func (p *PointAffine) isInSubGroupPornin() bool {
	subgroupInitOnce.Do(initSubgroupConstants)

	if isLowOrder(&p.X, &p.Y) {
		return p.IsZero()
	}

	u, w := edwardsToPorninMontgomery(p)

	var e fp.Element
	e.SetOne()

	// === First halving ===

	var us, ws fp.Element
	us.Double(&u)
	us.Double(&us)
	ws.Double(&w)

	var wp fp.Element
	if wp.Sqrt(&us) == nil {
		return false
	}
	var up, tmp fp.Element
	tmp.Square(&e)
	tmp.Mul(&tmp, &aPrimePorn)
	up.Sub(&us, &tmp)
	tmp.Mul(&wp, &ws)
	up.Sub(&up, &tmp)
	up.Halve()

	var sqrtResult fp.Element
	if sqrtResult.Sqrt(&up) != nil {
		w.Set(&sqrtResult)
	} else {
		var negUp fp.Element
		negUp.Neg(&up)
		sqrtResult.Sqrt(&negUp)
		wp.Mul(&wp, &up)
		wp.Neg(&wp)
		var e2 fp.Element
		e2.Square(&e)
		w.Mul(&sqrtResult, &sqrtMinBp)
		w.Mul(&w, &e2)
		e.Mul(&e, &up)
	}
	var e2 fp.Element
	e2.Square(&e)
	tmp.Mul(&e2, &minusAPornin)
	u.Square(&w)
	u.Add(&u, &tmp)
	tmp.Mul(&w, &wp)
	u.Sub(&u, &tmp)
	u.Halve()

	// === Second halving ===
	return u.Legendre() == 1
}

// IsInSubGroup tests subgroup membership using the fastest available method.
func (p *PointAffine) IsInSubGroup() bool {
	return p.isInSubGroupPornin()
}
