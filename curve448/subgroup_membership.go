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

// isInSubGroupQuarticExp tests subgroup membership using a single quartic
// residuosity check in Fp2 (exponentiation method).
//
// Instead of 1 halving (2 square roots) + 1 Legendre symbol, we compute
// the degree-4 Tate pairing t_4(T4, P) where T4 ∈ E[4](Fp2) \ E(Fp) is a
// non-rational 4-torsion point on the Montgomery curve. The pairing value
// ζ ∈ Fp2 satisfies t_4 = 1 iff P is in the prime-order subgroup.
//
// The quartic character χ₄(ζ) = ζ^((p²-1)/4) = 1 is checked by computing
// β = ζ^((p+1)/4) in Fp2 and verifying β ∈ Fp (i.e., β.A1 = 0).
//
// The Miller function value ζ is computed division-free as:
//
//	ℓ₁ = wP - λ·uP - C  (tangent line at T4 evaluated at P)
//	ζ = ℓ₁² · (uP - conj(u_{T2}))
//
// where λ, C are precomputed Fp2 constants and conj(u_{T2}) is the conjugate
// of the 2-torsion u-coordinate.
func (p *PointAffine) isInSubGroupQuarticExp1() bool {
	subgroupInitOnce.Do(initSubgroupConstants)

	if isLowOrder(&p.X, &p.Y) {
		return p.IsZero()
	}

	u, w := edwardsToPorninMontgomery(p)

	// l1 = (wP, 0) - lam * (uP, 0) - C
	// lam * (uP, 0) = (lam.A0*uP, lam.A1*uP)
	var l1 fp2.E2
	l1.A0.Mul(&quarticLam.A0, &u)
	l1.A1.Mul(&quarticLam.A1, &u)
	l1.A0.Sub(&w, &l1.A0)
	l1.A0.Sub(&l1.A0, &quarticC.A0)
	l1.A1.Neg(&l1.A1)
	l1.A1.Sub(&l1.A1, &quarticC.A1)

	// conj_v1 = (uP - 39080, 2*sqrt(39081))
	var conjV1 fp2.E2
	conjV1.A0.Sub(&u, &quarticConjRe)
	conjV1.A1.Set(&quarticConjIm)

	// alpha = l1^2 * conj_v1
	var alpha fp2.E2
	alpha.Square(&l1)
	alpha.Mul(&alpha, &conjV1)

	// beta = alpha^((p+1)/4) in Fp2
	var beta fp2.E2
	beta.ExpBySqrtPp1o4(&alpha)

	// chi_4(alpha) = 1 iff beta ∈ Fp iff beta.A1 = 0
	return beta.A1.IsZero()
}

// isInSubGroupQuarticExp2 tests subgroup membership using the g approach:
// g = conj(α)/α (on the torus), then g^((p+1)/4) with cyclotomic squarings,
// check result == 1.
func (p *PointAffine) isInSubGroupQuarticExp2() bool {
	subgroupInitOnce.Do(initSubgroupConstants)

	if isLowOrder(&p.X, &p.Y) {
		return p.IsZero()
	}

	u, w := edwardsToPorninMontgomery(p)

	// Same alpha construction as Exp1.
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

	// g approach: g = conj(alpha)/alpha (on the norm-1 torus).
	var conjAlpha, alphaInv, g fp2.E2
	conjAlpha.Conjugate(&alpha)
	alphaInv.Inverse(&alpha)
	g.Mul(&conjAlpha, &alphaInv)

	// g^((p+1)/4) with cyclotomic squarings, check == 1.
	var result fp2.E2
	result.CyclotomicExpBySqrtPp1o4(&g)
	return result.IsOne()
}

// isInSubGroupQuarticExp3 tests subgroup membership using the torus/Lucas approach:
// g = conj(α)/α (on the torus), compress to projective trace T/N ∈ Fp,
// then Montgomery ladder for (p+1)/4 = 2^222·(2^224-1) using t→t²−2
// for squarings and t_{m+n} = t_m·t_n − t_{m-n} for additions.
// Check trace == 2 at the end.
func (p *PointAffine) isInSubGroupQuarticExp3() bool {
	subgroupInitOnce.Do(initSubgroupConstants)

	if isLowOrder(&p.X, &p.Y) {
		return p.IsZero()
	}

	u, w := edwardsToPorninMontgomery(p)

	// Same alpha construction.
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
	// t = g + g⁻¹ = 2(a²-b²)/(a²+b²) for α = a+bi.
	// T = 2(a²-b²), N = a²+b² = Norm(α). No inversion needed.
	var a2, b2, T, N fp.Element
	a2.Square(&alpha.A0)
	b2.Square(&alpha.A1)
	T.Sub(&a2, &b2)
	T.Double(&T)
	N.Add(&a2, &b2)

	// Montgomery ladder for exponent e = (p+1)/4 = 2^222·(2^224-1).
	// We maintain (T_n, N_n, T_{n+1}, N_{n+1}) where U_k = T_k/N_k = trace(g^k).
	//
	// Doubling: U_{2k} = U_k² - 2
	//   T_{2k} = T_k² - 2·N_k²,  N_{2k} = N_k²
	//
	// Addition (using difference U_1):
	//   U_{k+1} = U_k·U_1 - U_{k-1}  (not directly usable in Montgomery ladder)
	//
	// Montgomery ladder step for bit b:
	//   if b=0: (U_n, U_{n+1}) → (U_{2n}, U_{2n+1})
	//   if b=1: (U_n, U_{n+1}) → (U_{2n+1}, U_{2n+2})
	// where U_{2n} = U_n²-2, U_{2n+1} = U_n·U_{n+1} - U_1, U_{2n+2} = U_{n+1}²-2
	//
	// Projective Montgomery ladder step:
	//   Double: (T, N) → (T²-2N², N²)
	//   DblAdd: T_{m+n}·N_1 = T_m·T_n·N_1 - ... needs care.
	//
	// For projective, maintain (Tn, Nn, Tn1, Nn1) with known (T1, N1):
	//   U_{2n}:   T' = Tn²-2Nn²,   N' = Nn²
	//   U_{2n+1}: T' = Tn·Tn1 - T1·Nn·Nn1/N1 ... requires division.
	//
	// Simpler: do the inversion once to get affine t = T/N, then run affine ladder.

	var nInv fp.Element
	nInv.Inverse(&N)
	var t fp.Element
	t.Mul(&T, &nInv) // t = trace(g) ∈ Fp

	// Affine Montgomery ladder for e = (p+1)/4.
	// (p+1)/4 in binary: 446 bits.
	// Maintain (u0, u1) = (trace(g^n), trace(g^{n+1})).
	// Bit=0: u0' = u0²-2, u1' = u0·u1 - t
	// Bit=1: u0' = u0·u1 - t, u1' = u1²-2

	// e = (p+1)/4 as big.Int
	var e big.Int
	e.Set(fp.Modulus())
	e.Add(&e, big.NewInt(1))
	e.Rsh(&e, 2)

	u0 := t // trace(g^1) = t
	var u1 fp.Element
	u1.Square(&t)
	u1.Sub(&u1, &twoFp) // trace(g^2) = t²-2

	for i := e.BitLen() - 2; i >= 0; i-- {
		var prod fp.Element
		prod.Mul(&u0, &u1)
		if e.Bit(i) == 0 {
			// u1 = u0*u1 - t, u0 = u0²-2
			u1.Sub(&prod, &t)
			u0.Square(&u0)
			u0.Sub(&u0, &twoFp)
		} else {
			// u0 = u0*u1 - t, u1 = u1²-2
			u0.Sub(&prod, &t)
			u1.Square(&u1)
			u1.Sub(&u1, &twoFp)
		}
	}

	// Check trace(g^e) == 2
	return u0.Equal(&twoFp)
}

// IsInSubGroup tests subgroup membership using the fastest available method.
func (p *PointAffine) IsInSubGroup() bool {
	return p.isInSubGroupQuarticExp1()
}
