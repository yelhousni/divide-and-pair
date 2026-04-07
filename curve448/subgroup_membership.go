package curve448

import (
	"math/big"
	"sync"

	fp "github.com/yelhousni/divide-and-pair/curve448/fp"
)

var (
	subgroupInitOnce sync.Once
	aMinusD          fp.Element // a - d = 39082
	aPornin          fp.Element // A = 2(a+d)
	bPornin          fp.Element // B = (a-d)^2
	aPrimePorn       fp.Element // A' = -2*A
	bPrimePorn       fp.Element // B' = A^2 - 4*B
	minusAPornin     fp.Element // -A = -2(a+d)
	sqrtMinBp        fp.Element // sqrt(-B')
	subgroupOrder    big.Int
)

func initSubgroupConstants() {
	initOnce.Do(initCurveParams)

	aMinusD.Sub(&curveParams.A, &curveParams.D)

	var aPlusD fp.Element
	aPlusD.Add(&curveParams.A, &curveParams.D)
	aPornin.Double(&aPlusD) // A = 2(a+d)

	bPornin.Square(&aMinusD) // B = (a-d)^2

	aPrimePorn.Double(&aPornin)
	aPrimePorn.Neg(&aPrimePorn) // A' = -2A

	var a2, b4 fp.Element
	a2.Square(&aPornin)
	b4.Double(&bPornin)
	b4.Double(&b4)
	bPrimePorn.Sub(&a2, &b4) // B' = A^2 - 4B

	minusAPornin.Neg(&aPornin) // -A

	// sqrt(-B') — exists because p ≡ 3 mod 4 and -B' is a QR
	var negBp fp.Element
	negBp.Neg(&bPrimePorn)
	sqrtMinBp.Sqrt(&negBp)

	subgroupOrder.Set(&curveParams.Order)
}

// isLowOrder checks if an affine point (X, Y) is a low-order point.
// For curve448 (a=1), the low-order points are exactly:
//
//	(0, 1) order 1, (0, -1) order 2, (1, 0) and (-1, 0) order 4.
//
// These are exactly the points with X=0 or Y=0.
func isLowOrder(X, Y *fp.Element) bool {
	return X.IsZero() || Y.IsZero()
}

// edwardsToPorninMontgomery converts an affine Edwards point (x, y) to
// Pornin's Montgomery (u, w) coordinates, batching the two inversions
// into one using Montgomery's trick.
func edwardsToPorninMontgomery(p *PointAffine) (u, w fp.Element) {
	// u = (a-d)(1+y)/(1-y), w = 2/x
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
// 2 halvings on the isogenous Montgomery curve + 1 Legendre symbol.
//
// The algorithm uses isogenies psi1, psi2 and an isomorphism iso such that
// 2*P = iso(psi2(psi1(P))). To halve, we invert iso, psi2, psi1.
// A point is in the prime-order subgroup iff it can be halved twice
// (checked via QR tests) and the final u-coordinate is a QR.
//
// Following crrl (Pornin's implementation), we track a scaling factor e
// to avoid divisions when up is not a QR.
func (p *PointAffine) isInSubGroupPornin() bool {
	subgroupInitOnce.Do(initSubgroupConstants)

	if isLowOrder(&p.X, &p.Y) {
		return p.IsZero()
	}

	u, w := edwardsToPorninMontgomery(p)

	// We are on Curve(A, B) with e = 1 (affine).
	// Actual curve is Curve(A*e^2, B*e^4).
	var e fp.Element
	e.SetOne()

	// === First halving ===

	// Inverse iso: (u, w) -> (4u, 2w) on Curve(As*e^2, Bs*e^4)
	var us, ws fp.Element
	us.Double(&u)
	us.Double(&us) // 4u
	ws.Double(&w)  // 2w

	// Inverse psi2: need sqrt(us). If us is not a QR, point cannot be halved.
	var wp fp.Element
	if wp.Sqrt(&us) == nil {
		return false
	}
	// up = (us - A'*e^2 - wp*ws) / 2
	var up, tmp fp.Element
	tmp.Square(&e)
	tmp.Mul(&tmp, &aPrimePorn)
	up.Sub(&us, &tmp)
	tmp.Mul(&wp, &ws)
	up.Sub(&up, &tmp)
	up.Halve()

	// Inverse psi1.
	// If up is a QR: w = sqrt(up), proceed normally.
	// If up is NQR (p ≡ 3 mod 4 → -up is QR):
	//   Switch to isomorphic curve:
	//     wp' = -wp*up, e' = e*up
	//     w = sqrt(-up) * sqrt(-Bp) * e^2
	var sqrtResult fp.Element
	if sqrtResult.Sqrt(&up) != nil {
		// up is QR
		w.Set(&sqrtResult)
	} else {
		// up is NQR: sqrt(-up) exists
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
	// u = (w^2 + (-A)*e^2 - w*wp) / 2
	var e2 fp.Element
	e2.Square(&e)
	tmp.Mul(&e2, &minusAPornin)
	u.Square(&w)
	u.Add(&u, &tmp)
	tmp.Mul(&w, &wp)
	u.Sub(&u, &tmp)
	u.Halve()

	// === Second halving ===
	// Only need Legendre symbol of u to check if it can be halved again.
	return u.Legendre() == 1
}

// IsInSubGroup tests subgroup membership using the fastest available method.
func (p *PointAffine) IsInSubGroup() bool {
	return p.isInSubGroupPornin()
}
