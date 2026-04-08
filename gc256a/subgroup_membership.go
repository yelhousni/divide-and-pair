package gc256a

import (
	"math/big"
	"sync"

	fp "github.com/yelhousni/divide-and-pair/gc256a/fp"
	fp2 "github.com/yelhousni/divide-and-pair/gc256a/fp2"
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

	twoFp fp.Element

	quarticLam    fp2.E2
	quarticC      fp2.E2
	quarticConjRe fp.Element
	quarticConjIm fp.Element
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

	twoFp.SetUint64(2)

	quarticLam.A0.SetString("98978386755764401450203974662446294613320410218429152052753462557826881889291")
	quarticLam.A1.SetString("21937508357931127973100225998923042910701061142013484653662311276571639819534")
	quarticC.A0.SetString("103296244185588961271874791050375166492124013234565164938412440660178304685659")
	quarticC.A1.SetString("43637342742974478728501216132433420977611717423370731437025443950094608957872")
	quarticConjRe.SetString("113067675126841589491736716507523150207271257787167487230024979784498777963931")
	quarticConjIm.SetString("3525713647912769663176164827195583175754606956748302963316409412254261816257")
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

// halvePornin performs one division-free halving on Pornin's Montgomery
// coordinates (u, w) with scaling factor e. Returns false if not halvable.
// For GC256A (p ≡ 3 mod 4), when up is NQR, -up is QR.
func halvePornin(u, w, e *fp.Element) bool {
	var us, ws fp.Element
	us.Double(u)
	us.Double(&us)
	ws.Double(w)

	var wp fp.Element
	if wp.Sqrt(&us) == nil {
		return false
	}

	var up, tmp fp.Element
	tmp.Square(e)
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
		e2.Square(e)
		w.Mul(&sqrtResult, &sqrtMinBp)
		w.Mul(w, &e2)
		e.Mul(e, &up)
	}

	var e2 fp.Element
	e2.Square(e)
	tmp.Mul(&e2, &minusAPornin)
	u.Square(w)
	u.Add(u, &tmp)
	tmp.Mul(w, &wp)
	u.Sub(u, &tmp)
	u.Halve()

	return true
}

// isInSubGroupNaive tests subgroup membership by scalar multiplication by ℓ.
func (p *PointAffine) isInSubGroupNaive() bool {
	subgroupInitOnce.Do(initSubgroupConstants)
	var res PointAffine
	res.ScalarMultiplication(p, &subgroupOrder)
	return res.IsZero()
}

// isInSubGroupPornin tests subgroup membership using Pornin's method:
// 1 halving + 1 Legendre symbol, division-free with scaling factor e.
func (p *PointAffine) isInSubGroupPornin() bool {
	subgroupInitOnce.Do(initSubgroupConstants)

	if isLowOrder(&p.X, &p.Y) {
		return p.IsZero()
	}

	u, w := edwardsToPorninMontgomery(p)

	var e fp.Element
	e.SetOne()

	if !halvePornin(&u, &w, &e) {
		return false
	}

	return u.Legendre() == 1
}

// pracOpsQuartic encodes the PRAC differential addition chain for (p+1)/4.
// Cost: 470 field ops (vs 506 for binary ladder, 7% saving).
var pracOpsQuartic = [...]byte{
	10, 3, 10, 3, 0, 3, 10, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, // 0-19
	0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, // 20-39
	0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, // 40-59
	0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, // 60-79
	1, 3, 0, 7, 3, 0, 3, 2, 0, 9, 4, 4, 3, 3, 3, 0, 3, 0, 1, 3, // 80-99
	0, 3, 0, 4, 5, 4, 5, 4, 3, 3, 0, 3, 0, 5, 5, 4, 4, 5, 4, 4, // 100-119
	3, 3, 0, 3, 3, 0, 3, 0, 3, 3, 3, 0, 5, 3, 3, 3, 0, 3, 0, 3, // 120-139
	2, 0, 4, 4, 3, 0, 3, 0, 3, 3, 0, 5, 3, 3, 0, 3, 3, 0, 3, 0, // 140-159
	3, 3, 0, 3, 0, 3, 3, 0, 3, 3, 3, 0, 3, 0, 3, 2, 0, 4, 4, 3, // 160-179
	0, 3, 0, 3, 3, 3, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 4, 4, 5, // 180-199
	5, 5, 4, 4, 5, 5, 4, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 5, // 200-219
	4, 5, 4, 5, 4, 5, 5, 5, 4, 4, 4, 4, 5, 4, 5, 4, 5, 5, 4, 4, // 220-239
	5, 4, 5, 4, 5, 4, 5, 4, 5, 4, 4, 4, 4, 5, 5, 4, 5, 5, 4, 4, // 240-259
	5, 5, 4, 5, 4, 4, 4, 4, 5, 5, 5, 5, 4, 5, 4, 5, 4, 4, 4, 5, // 260-279
	5, 4, 5, 4, 4, 5, 5, 4, 5, 5, 5, 5, 5, 3, 3, 3, 0, 3, 0, 3, // 280-299
	0, 3, 3, 3, 0, 4, 3, 0, 3, 0, 5, 4, 3, 3, 0, 4, 3, 3, 0, 1, // 300-319
	1, 3, 0, 3, 0, 2, 0, 9, 4, 4, 3, 3, 3, 0, 3, 3, 0, 3, 0, 3, // 320-339
	0, 3, 0, 1, 1, 3, 0, 4, 3, 3, 3, 0, 3, 0, 5, 5, 3, 3, 0, 3, // 340-359
	1, 3, 0, 5, 4, 5, 5, 5, 4, 5, 5, 5, 3, 3, 0, 3, 0, 6, 3, 3, // 360-379
	0, 3, 3, 0, 3, 0, 3, 3, 0, 3, 3, 0, 3, 0, 3, 3, 10, // 380-396
}

// pracLucasV evaluates the Lucas V-sequence using a PRAC differential
// addition chain. Returns V_n where n is encoded in pracOps.
func pracLucasV(t *fp.Element, pracOps []byte) fp.Element {
	var A, B, C fp.Element
	var T1, T2, T3 fp.Element
	A.Set(t)
	B.Set(t)
	C.Set(&twoFp)

	for _, op := range pracOps {
		switch op {
		case 0:
			A, B = B, A
		case 1:
			T1.Mul(&A, &B)
			T1.Sub(&T1, &C)
			T2.Mul(&T1, &A)
			T2.Sub(&T2, &B)
			B.Mul(&T1, &B)
			B.Sub(&B, &A)
			A.Set(&T2)
		case 2, 4:
			T1.Mul(&A, &B)
			T1.Sub(&T1, &C)
			A.Square(&A)
			A.Sub(&A, &twoFp)
			B.Set(&T1)
		case 3:
			T1.Mul(&A, &B)
			T1.Sub(&T1, &C)
			C.Set(&B)
			B.Set(&T1)
		case 5:
			T1.Mul(&A, &C)
			T1.Sub(&T1, &B)
			A.Square(&A)
			A.Sub(&A, &twoFp)
			C.Set(&T1)
		case 6:
			T1.Mul(&A, &B)
			T1.Sub(&T1, &C)
			T2.Square(&A)
			T2.Sub(&T2, &twoFp)
			T3.Mul(&T2, &A)
			T3.Sub(&T3, &A)
			A.Set(&T3)
			T3.Mul(&T2, &T1)
			T3.Sub(&T3, &C)
			C.Set(&B)
			B.Set(&T3)
		case 7:
			T1.Mul(&A, &B)
			T1.Sub(&T1, &C)
			T2.Square(&A)
			T2.Sub(&T2, &twoFp)
			T3.Mul(&T2, &A)
			T3.Sub(&T3, &A)
			T2.Mul(&T1, &A)
			T2.Sub(&T2, &B)
			A.Set(&T3)
			B.Set(&T2)
		case 8:
			T1.Mul(&A, &B)
			T1.Sub(&T1, &C)
			T2.Square(&A)
			T2.Sub(&T2, &twoFp)
			T3.Mul(&T2, &A)
			T3.Sub(&T3, &A)
			C.Mul(&A, &C)
			C.Sub(&C, &B)
			A.Set(&T3)
			B.Set(&T1)
		case 9:
			T1.Mul(&C, &B)
			T1.Sub(&T1, &A)
			B.Square(&B)
			B.Sub(&B, &twoFp)
			C.Set(&T1)
		case 10:
			A.Mul(&A, &B)
			A.Sub(&A, &C)
			B.Set(&A)
			C.Set(&twoFp)
		}
	}

	return A
}

// isInSubGroupQuartic tests subgroup membership using the torus approach
// with a PRAC differential addition chain (Montgomery 1992) for the Lucas
// V-sequence evaluation.
func (p *PointAffine) isInSubGroupQuartic() bool {
	subgroupInitOnce.Do(initSubgroupConstants)

	if isLowOrder(&p.X, &p.Y) {
		return p.IsZero()
	}

	u, w := edwardsToPorninMontgomery(p)

	var l1 fp2.E2
	l1.A0.Mul(&quarticLam.A0, &u)
	l1.A1.Mul(&quarticLam.A1, &u)
	l1.A0.Sub(&w, &l1.A0)
	l1.A0.Sub(&l1.A0, &quarticC.A0)
	l1.A1.Neg(&l1.A1)
	l1.A1.Sub(&l1.A1, &quarticC.A1)

	var conjV1 fp2.E2
	conjV1.A0.Sub(&u, &quarticConjRe)
	conjV1.A1.Set(&quarticConjIm)

	var alpha fp2.E2
	alpha.Square(&l1)
	alpha.Mul(&alpha, &conjV1)

	var a2, b2, T, N fp.Element
	a2.Square(&alpha.A0)
	b2.Square(&alpha.A1)
	T.Sub(&a2, &b2)
	T.Double(&T)
	N.Add(&a2, &b2)

	var nInv, t fp.Element
	nInv.Inverse(&N)
	t.Mul(&T, &nInv)

	result := pracLucasV(&t, pracOpsQuartic[:])
	return result.Equal(&twoFp)
}

// IsInSubGroup tests subgroup membership using the fastest available method.
func (p *PointAffine) IsInSubGroup() bool {
	return p.isInSubGroupQuartic()
}
