package jubjub

import (
	"math/big"
	"sync"

	fp "github.com/yelhousni/divide-and-pair/jubjub/fp"
)

var (
	subgroupInitOnce sync.Once
	aPornin          fp.Element // A = 2(a+d)
	bPornin          fp.Element // B = (a-d)^2
	aPrimePorn       fp.Element // A' = -2*A
	bPrimePorn       fp.Element // B' = A^2 - 4*B
	aMinusD          fp.Element // a - d
	sqrtMinOne       fp.Element // i = sqrt(-1) mod p
	twoI             fp.Element // 2i
	nqr              fp.Element // smallest QNR in Fp (= 5 for BLS12-381 Fr)
	sqrtNqrBp        fp.Element // sqrt(nqr * B')
	subgroupOrder    big.Int

	// Octic constants (Weierstrass model for degree-8 Miller function)
	aW         fp.Element // Weierstrass a coefficient
	three      fp.Element
	octicXT8   fp.Element // 8-torsion point X on Weierstrass
	octicYT8   fp.Element // 8-torsion point Y on Weierstrass
	octicLamT8 fp.Element // tangent slope at T8
	octicXT4   fp.Element // 4-torsion point X on Weierstrass
	octicYT4   fp.Element // 4-torsion point Y on Weierstrass
	octicLamT4 fp.Element // tangent slope at T4
	octicXT2   fp.Element // 2-torsion point X on Weierstrass
	Adiv3      fp.Element // A/3 for Edwards-to-Weierstrass conversion
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

	// sqrt(-1) — exists since p ≡ 1 mod 4
	var negOne fp.Element
	negOne.SetOne()
	negOne.Neg(&negOne)
	sqrtMinOne.Sqrt(&negOne)
	twoI.Double(&sqrtMinOne)

	// For this field (BLS12-381 Fr), p ≡ 1 mod 8, so 2 is a QR.
	// The crrl sqrt(2*up) trick doesn't work. Instead, we use a fixed
	// quadratic non-residue g = 5. When up is NQR, g*up is QR.
	// Precompute sqrt(g * B').
	nqr.SetUint64(5)
	var nqrBp fp.Element
	nqrBp.Mul(&nqr, &bPrimePorn)
	sqrtNqrBp.Sqrt(&nqrBp)

	subgroupOrder.Set(&curveParams.Order)

	// Octic constants: Weierstrass model Y^2 = X^3 + aW*X + bW
	// via X = u + A/3, Y = u*w (where u, w are Pornin Montgomery coords)
	three.SetUint64(3)
	var threeInv fp.Element
	threeInv.Inverse(&three)
	Adiv3.Mul(&aPornin, &threeInv)

	var A2div3 fp.Element
	A2div3.Square(&aPornin)
	A2div3.Mul(&A2div3, &threeInv)
	aW.Sub(&bPornin, &A2div3)

	// 8-torsion point T8, [2]T8=T4, [4]T8=T2 on Weierstrass (from Sage)
	octicXT8.SetString("38074089473419775066521122308184450788044636445201246956909382033745382347397")
	octicYT8.SetString("30713742701501207659057220608756960799104471161253953255407673889302209579815")
	octicLamT8.SetString("48568263891961809311880741918802055033685793886551129072617710864462986515188")
	octicXT4.SetString("11059612379481747039899142612799695948580372366091172512139820602662565702425")
	octicYT4.SetString("4853940988772544534178633755630548083074718905337523302547095884433093026740")
	octicXT2.SetString("30316650416162696399649455282586573940529807768345292798324017494613449779659")

	// Tangent slope at T4
	var XT42, threeXT42, num, twoYT4, twoYT4Inv fp.Element
	XT42.Square(&octicXT4)
	threeXT42.Mul(&three, &XT42)
	num.Add(&threeXT42, &aW)
	twoYT4.Double(&octicYT4)
	twoYT4Inv.Inverse(&twoYT4)
	octicLamT4.Mul(&num, &twoYT4Inv)

}

// isLowOrder checks if an affine point (X, Y) is a low-order point.
// For JubJub (a=-1, cofactor 8), low-order points in E[8] are exactly
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

// edwardsToPorninMontgomeryScaled maps an affine point to scaled Montgomery
// coordinates without the initial inversion. The represented point is on
// Curve(A*e^2, B*e^4) with:
//
//	e = x(1-y)
//	u = (a-d)(1+y)x^2(1-y)
//	w = 2(1-y)
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

// edwardsToPorninMontgomery converts an affine Edwards point (x, y) to
// Pornin's Montgomery (u, w) coordinates, batching the two inversions
// into one using Montgomery's trick.
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

// halvePornin performs one halving step, division-free with scaling factor e.
// When up is NQR, we use the fixed NQR g (= 5 for this field) instead of 2
// (since 2 is QR for p ≡ 1 mod 8).
func halvePornin(u, w *fp.Element, e *fp.Element) (uOut, wOut fp.Element, eOut fp.Element, ok bool) {
	var us, ws fp.Element
	us.Double(u)
	us.Double(&us)
	ws.Double(w)

	var wp fp.Element
	if wp.Sqrt(&us) == nil {
		return fp.Element{}, fp.Element{}, fp.Element{}, false
	}

	var up, tmp, e2 fp.Element
	e2.Square(e)
	tmp.Mul(&e2, &aPrimePorn)
	up.Sub(&us, &tmp)
	tmp.Mul(&wp, &ws)
	up.Sub(&up, &tmp)
	up.Halve()

	eOut.Set(e)

	var sqrtUp fp.Element
	if sqrtUp.Sqrt(&up) != nil {
		// QR case
		wOut.Set(&sqrtUp)
		tmp.Mul(&e2, &aPornin)
		uOut.Square(&wOut)
		uOut.Sub(&uOut, &tmp)
		tmp.Mul(&wOut, &wp)
		uOut.Sub(&uOut, &tmp)
		uOut.Halve()
	} else {
		// NQR case: g*up is QR (g is our fixed NQR)
		var gUp fp.Element
		gUp.Mul(&nqr, &up)
		var tt fp.Element
		tt.Sqrt(&gUp)

		// Update isomorphism: same structure as crrl but with g instead of 2
		wp.Mul(&wp, &tt)
		wOut.Mul(&sqrtNqrBp, &e2)
		eOut.Mul(&eOut, &tt)

		e2.Square(&eOut)
		tmp.Mul(&e2, &aPornin)
		uOut.Square(&wOut)
		uOut.Sub(&uOut, &tmp)
		tmp.Mul(&wOut, &wp)
		uOut.Sub(&uOut, &tmp)
		uOut.Halve()

		wOut.Neg(&wOut)
	}

	return uOut, wOut, eOut, true
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

// isInSubGroupPornin tests subgroup membership using Pornin's method:
// 2 halvings + 1 Legendre symbol (3rd halving check).
func (p *PointAffine) isInSubGroupPornin() bool {
	subgroupInitOnce.Do(initSubgroupConstants)

	if isLowOrder(&p.X, &p.Y) {
		return p.IsZero()
	}

	u, w, e := edwardsToPorninMontgomeryScaled(p)

	for range 2 {
		var uNew, wNew, eNew fp.Element
		var ok bool
		uNew, wNew, eNew, ok = halvePornin(&u, &w, &e)
		if !ok {
			return false
		}
		u, w, e = uNew, wNew, eNew
	}

	return u.Legendre() == 1
}

// edwardsToWeierstrass converts an affine Edwards point (x, y) to
// the Weierstrass model Y^2 = X^3 + aW*X + bW via the Pornin Montgomery
// intermediate: X = u + A/3, Y = u*w.
func edwardsToWeierstrass(p *PointAffine) (X, Y fp.Element) {
	u, w := edwardsToPorninMontgomery(p)
	X.Add(&u, &Adiv3)
	Y.Mul(&u, &w)
	return
}

// edwardsToWeierstrassScaled maps an affine Edwards point to scaled
// Weierstrass coordinates without the initial inversion. If e = x(1-y),
// then the affine Weierstrass coordinates satisfy:
//
//	X = Xs / e^2
//	Y = Ys / e^3
func edwardsToWeierstrassScaled(p *PointAffine) (Xs, Ys, e fp.Element) {
	U, W, e := edwardsToPorninMontgomeryScaled(p)

	var e2 fp.Element
	e2.Square(&e)

	Xs.Mul(&Adiv3, &e2)
	Xs.Add(&Xs, &U)
	Ys.Mul(&U, &W)
	return
}

// isInSubGroupOcticExp tests subgroup membership using the divide-and-pair
// method with 0 halvings + 1 octic residuosity check.
//
// Since d = gcd(8, p-1) = 8 for the JubJub field (p ≡ 1 mod 8), no halvings
// are needed. The degree-8 Miller function f_{8,T₈}(Q) is evaluated on the
// Weierstrass model using a 3-step doubling chain T₈ → T₄ → T₂ → ∞, and
// the octic symbol χ₈(f) = f^((p-1)/8) is checked via exponentiation.
func (p *PointAffine) isInSubGroupOcticExp() bool {
	subgroupInitOnce.Do(initSubgroupConstants)

	if isLowOrder(&p.X, &p.Y) {
		return p.IsZero()
	}

	XQ, YQ, e := edwardsToWeierstrassScaled(p)
	var e2, e3 fp.Element
	e2.Square(&e)
	e3.Mul(&e2, &e)

	// Miller function f_{8,T₈}(Q):
	// R = T₈, f = 1. Three doubling steps for 8 = 1000₂.
	//
	// Step 1 (R=T₈→T₄): use scaled line values. Since X = XQ/e² and Y = YQ/e³,
	// the common denominator contributes an 8th power overall and does not
	// affect the final octic symbol.
	var xd, ell1, v1 fp.Element
	var xt8e2, yt8e3 fp.Element
	xt8e2.Mul(&octicXT8, &e2)
	yt8e3.Mul(&octicYT8, &e3)
	xd.Sub(&XQ, &xt8e2)
	ell1.Mul(&octicLamT8, &xd)
	ell1.Mul(&ell1, &e)
	ell1.Sub(&YQ, &ell1)
	ell1.Sub(&ell1, &yt8e3)

	var xt4e2 fp.Element
	xt4e2.Mul(&octicXT4, &e2)
	v1.Sub(&XQ, &xt4e2)

	// Step 2 (R=T₄→T₂): f = f² · ℓ_{T₄}(Q) / v_{T₂}(Q)
	var xd2, ell2, v2 fp.Element
	var yt4e3, xt2e2 fp.Element
	yt4e3.Mul(&octicYT4, &e3)
	xt2e2.Mul(&octicXT2, &e2)
	xd2.Sub(&XQ, &xt4e2)
	ell2.Mul(&octicLamT4, &xd2)
	ell2.Mul(&ell2, &e)
	ell2.Sub(&YQ, &ell2)
	ell2.Sub(&ell2, &yt4e3)
	v2.Sub(&XQ, &xt2e2)

	// Step 3 (R=T₂→∞): f = f² · (X_Q - X_{T₂})
	// At the 2-torsion point T₂ = (X_{T₂}, 0), the tangent is vertical,
	// so g₃ = X_Q - X_{T₂}. The vertical at ∞ is 1.
	// Note: v2 = X_Q - X_{T₂} = g₃.

	// Full Miller: f = g1⁴ · g2² · g3
	// where g1 = ell1/v1, g2 = ell2/v2, g3 = v2
	// Division-free: multiply by v1⁸ · v2⁸ (both 8th powers):
	//   f · v1⁸ · v2⁸ = ell1⁴ · v1⁴ · ell2² · v2⁶ · v2
	//                  = ell1⁴ · v1⁴ · ell2² · v2⁷
	var ell14, v14, ell22, v27, result fp.Element
	ell14.Square(&ell1)
	ell14.Square(&ell14) // ell1⁴
	v14.Square(&v1)
	v14.Square(&v14)    // v1⁴
	ell22.Square(&ell2) // ell2²
	v27.Square(&v2)
	v27.Mul(&v27, &v2) // v2³
	v27.Square(&v27)
	v27.Mul(&v27, &v2) // v2⁷

	result.Mul(&ell14, &v14)
	result.Mul(&result, &ell22)
	result.Mul(&result, &v27)

	// χ₈(f) = result^((p-1)/8). Check if it equals 1.
	var chi fp.Element
	chi.ExpByOcticExp(result)
	return chi.IsOne()
}

// IsInSubGroup tests subgroup membership using the fastest available method.
func (p *PointAffine) IsInSubGroup() bool {
	return p.isInSubGroupOcticExp()
}
