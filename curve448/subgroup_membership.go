package curve448

import (
	"math/big"
	"sync"

	fp "github.com/yelhousni/divide-and-pair/curve448/fp"
	"github.com/yelhousni/divide-and-pair/curve448/fp2"
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

	twoFp fp.Element // precomputed constant 2 ∈ Fp

	// Quartic test constants (degree-4 Tate pairing over Fp2).
	// T4 = (uT4, wT4) is a 4-torsion point on the Montgomery curve over Fp2.
	// lam is the tangent slope at T4.
	// C = wT4 + lam*uT4 (precomputed to simplify line evaluation).
	// conjUT2Im = 2*sqrt(39081) (imaginary part of conj(uT2), the conjugate 2-torsion u-coordinate).
	quarticLam    fp2.E2
	quarticC      fp2.E2
	quarticConjRe fp.Element // 39080
	quarticConjIm fp.Element // 2*sqrt(39081) = imaginary part of conj(u_{T2})
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

	twoFp.SetUint64(2)

	// Quartic test constants.
	// These come from a 4-torsion point T4 ∈ E(Fp2) \ E(Fp) on the Montgomery curve,
	// obtained by halving the non-rational 2-torsion T2 = (39080 + 2√39081·i, 0).
	quarticLam.A0.SetString("409466785180518956961611265803290831348604928765886125521636546093507230070268408209394576752377465036003676898038424039336503348623066")
	quarticLam.A1.SetString("607666065788055156040088647550522819510877211007230551732676526742583204812484502734497541805514130443469320956493054479303320776889788")

	// C = wT4 - lam*uT4 (precomputed constant in the tangent line evaluation).
	quarticC.A0.SetString("407637470893696240970274054411012998449963411119155090599355826322254189874134214156025763704645782904595902611200522707074365055469475")
	quarticC.A1.SetString("506581704541428460739451059287942017578380987414684991673455715920297891990135472012827895268542298278523614843293368590301287473065690")

	quarticConjRe.SetUint64(39080)
	// Note: conj(uT2) = (39080, -2s). We compute uP - conj(uT2) = (uP-39080, +2s).
	// So quarticConjIm stores +2*sqrt(39081), the imaginary part of (uP - conj(uT2)).
	var sqrt39081 fp.Element
	sqrt39081.SetString("627894490647874670780146803011075515225223784391788159207390309582568626050729514829594252134780030556161172229750791477826575600769225")
	quarticConjIm.Double(&sqrt39081) // 2*sqrt(39081)
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

// edwardsToPorninMontgomeryScaled maps an affine Edwards point to scaled
// Montgomery coordinates without the initial inversion.
func edwardsToPorninMontgomeryScaled(p *PointAffine) (u, w, e fp.Element) {
	var oneMinusY, onePlusY fp.Element
	var one fp.Element
	one.SetOne()
	oneMinusY.Sub(&one, &p.Y)
	onePlusY.Add(&one, &p.Y)
	e.Mul(&p.X, &oneMinusY)
	u.Mul(&aMinusD, &onePlusY)
	u.Mul(&u, &p.X)
	u.Mul(&u, &e)
	w.Double(&oneMinusY)
	return
}

// edwardsToPorninMontgomery keeps the affine Montgomery coordinates used by
// the quartic torus path.
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

func sqrtExtCurve448(z, x *fp.Element) bool {
	if z.Sqrt(x) != nil {
		return true
	}
	var negX fp.Element
	negX.Neg(x)
	z.Sqrt(&negX)
	return false
}

// halvePornin performs one division-free halving on Pornin's Montgomery
// coordinates (u, w) with scaling factor e. Returns false if not halvable.
// For curve448 (p ≡ 3 mod 4), when up is NQR, -up is QR.
func halvePornin(u, w, e *fp.Element) bool {
	// Inverse iso: (u, w) -> (4u, 2w)
	var us, ws fp.Element
	us.Double(u)
	us.Double(&us)
	ws.Double(w)

	// Inverse psi2: need sqrt(us). If us is not a QR, point cannot be halved.
	var wp fp.Element
	if wp.Sqrt(&us) == nil {
		return false
	}

	// up = (us - A'*e^2 - wp*ws) / 2
	var up, tmp fp.Element
	tmp.Square(e)
	tmp.Mul(&tmp, &aPrimePorn)
	up.Sub(&us, &tmp)
	tmp.Mul(&wp, &ws)
	up.Sub(&up, &tmp)
	up.Halve()

	// Inverse psi1.
	var sqrtResult fp.Element
	if sqrtExtCurve448(&sqrtResult, &up) {
		w.Set(&sqrtResult)
	} else {
		wp.Mul(&wp, &up)
		wp.Neg(&wp)
		var e2 fp.Element
		e2.Square(e)
		w.Mul(&sqrtResult, &sqrtMinBp)
		w.Mul(w, &e2)
		e.Mul(e, &up)
	}

	// u = (w^2 + (-A)*e^2 - w*wp) / 2
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
	var proj PointProj
	proj.FromAffine(p)
	var res PointProj
	res.ScalarMultiplication(&proj, &subgroupOrder)
	return res.IsZero()
}

// isInSubGroupPornin tests subgroup membership using Pornin's method
// (https://eprint.iacr.org/2022/1164):
// 1 halving on the isogenous Montgomery curve + 1 Legendre symbol.
func (p *PointAffine) isInSubGroupPornin() bool {
	subgroupInitOnce.Do(initSubgroupConstants)

	if isLowOrder(&p.X, &p.Y) {
		return p.IsZero()
	}

	u, w, e := edwardsToPorninMontgomeryScaled(p)

	if !halvePornin(&u, &w, &e) {
		return false
	}

	return u.Legendre() == 1
}

// pracOpsQuartic encodes the PRAC differential addition chain for (p+1)/4.
// Cost: 632 field ops (vs 890 for binary ladder, 29% saving).
var pracOpsQuartic = [592]byte{
	10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
	10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
	10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
	10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
	10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
	10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
	10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
	10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
	10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
	10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
	10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
	10, 10, 3, 10, 3, 0, 3, 10, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0,
	3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0,
	3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0,
	3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0,
	3, 1, 3, 0, 7, 3, 0, 3, 2, 0, 9, 4, 4, 3, 3, 3, 0, 3, 0, 1,
	3, 0, 3, 0, 4, 5, 4, 5, 4, 3, 3, 0, 3, 0, 5, 5, 4, 4, 5, 4,
	4, 3, 3, 0, 3, 3, 0, 3, 0, 3, 3, 3, 0, 5, 3, 3, 3, 0, 3, 0,
	3, 2, 0, 4, 4, 3, 0, 3, 0, 3, 3, 0, 5, 3, 3, 0, 3, 3, 0, 3,
	0, 3, 3, 0, 3, 0, 3, 3, 0, 3, 3, 3, 0, 3, 0, 3, 2, 0, 4, 4,
	3, 0, 3, 0, 3, 3, 3, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 4, 4,
	5, 5, 5, 4, 4, 5, 5, 4, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5,
	5, 4, 5, 4, 5, 4, 5, 5, 5, 4, 4, 4, 4, 5, 4, 5, 4, 5, 5, 4,
	4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4, 4, 4, 4, 5, 5, 4, 5, 5, 4,
	4, 5, 5, 4, 5, 4, 4, 4, 4, 5, 5, 5, 3, 2, 0, 9, 9, 4, 5, 4,
	5, 4, 4, 3, 3, 0, 3, 1, 3, 0, 3, 3, 0, 3, 0, 3, 3, 0, 3, 0,
	3, 0, 3, 3, 0, 3, 3, 3, 0, 3, 0, 3, 0, 3, 0, 4, 3, 3, 0, 3,
	3, 0, 4, 3, 3, 0, 3, 0, 3, 0, 3, 3, 0, 3, 0, 3, 0, 3, 0, 6,
	0, 3, 0, 3, 0, 3, 1, 3, 0, 3, 3, 3, 0, 3, 3, 0, 3, 3, 0, 3,
	1, 3, 0, 3, 3, 0, 4, 4, 4, 4, 3, 10,
}

// pracLucasV evaluates the Lucas V-sequence V_{(p+1)/4}(t) using a PRAC
// differential addition chain. Returns V_n where n is encoded in pracOps.
func pracLucasV(t *fp.Element, pracOps []byte) fp.Element {
	var A, B, C fp.Element
	var T1, T2, T3 fp.Element
	A.Set(t)
	B.Set(t)
	C.Set(&twoFp)

	for _, op := range pracOps {
		switch op {
		case 0: // SWAP
			A, B = B, A
		case 1: // CASE1: 3 ops
			T1.Mul(&A, &B)
			T1.Sub(&T1, &C)
			T2.Mul(&T1, &A)
			T2.Sub(&T2, &B)
			B.Mul(&T1, &B)
			B.Sub(&B, &A)
			A.Set(&T2)
		case 2, 4: // CASE2, CASE4: 2 ops
			T1.Mul(&A, &B)
			T1.Sub(&T1, &C)
			A.Square(&A)
			A.Sub(&A, &twoFp)
			B.Set(&T1)
		case 3: // CASE3: 1 op
			T1.Mul(&A, &B)
			T1.Sub(&T1, &C)
			C.Set(&B)
			B.Set(&T1)
		case 5: // CASE5: 2 ops
			T1.Mul(&A, &C)
			T1.Sub(&T1, &B)
			A.Square(&A)
			A.Sub(&A, &twoFp)
			C.Set(&T1)
		case 6: // CASE6: 4 ops
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
		case 7: // CASE7: 4 ops
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
		case 8: // CASE8: 4 ops
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
		case 9: // CASE9: 2 ops
			T1.Mul(&C, &B)
			T1.Sub(&T1, &A)
			B.Square(&B)
			B.Sub(&B, &twoFp)
			C.Set(&T1)
		case 10: // FINAL: 1 op
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
// V-sequence evaluation. PRAC uses 632 field ops vs 890 for the binary ladder.
func (p *PointAffine) isInSubGroupQuartic() bool {
	subgroupInitOnce.Do(initSubgroupConstants)

	if isLowOrder(&p.X, &p.Y) {
		return p.IsZero()
	}

	u, w := edwardsToPorninMontgomery(p)

	// Alpha construction.
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

	// Projective trace of g = conj(α)/α:
	var a2, b2, T, N fp.Element
	a2.Square(&alpha.A0)
	b2.Square(&alpha.A1)
	T.Sub(&a2, &b2)
	T.Double(&T)
	N.Add(&a2, &b2)

	// Affine trace: t = T/N
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
