package fourq

import (
	"math/big"
	"sync"

	fp "github.com/yelhousni/divide-and-pair/fourq/fp"
	fp2 "github.com/yelhousni/divide-and-pair/fourq/fp2"
)

var (
	smtInitOnce sync.Once

	// Weierstrass model constants (over Fp2)
	aW    fp2.E2
	Adiv3 fp2.E2

	// 8-torsion point T8 on Weierstrass
	XT8, YT8, lamT8 fp2.E2
	XT4, YT4, lamT4 fp2.E2
	XT2             fp2.E2

	// Degree-7 Miller loop precomputed constants:
	// S7 (7-torsion), [2]S7, [3]S7 and their slopes
	XS7, YS7 fp2.E2
	lamDbl1  fp2.E2 // tangent slope at S7
	X2S, Y2S fp2.E2 // [2]S7
	lamAdd1  fp2.E2 // chord slope S7, [2]S7
	X3S, Y3S fp2.E2 // [3]S7
	lamDbl2  fp2.E2 // tangent slope at [3]S7
	X6S      fp2.E2 // [6]S7 = -S7 (X only, Y is -YS7)

	three fp2.E2
)

func initSMTConstants() {
	initOnce.Do(initCurveParams)

	amd := new(fp2.E2).Sub(&curveParams.A, &curveParams.D)
	apd := new(fp2.E2).Add(&curveParams.A, &curveParams.D)

	var Amg, Bmg fp2.E2
	Amg.Add(apd, apd)
	Bmg.Square(amd)
	three.A0.SetUint64(3)

	var threeInv fp2.E2
	threeInv.Inverse(&three)
	Adiv3.Mul(&Amg, &threeInv)

	var A2div3 fp2.E2
	A2div3.Square(&Amg)
	A2div3.Mul(&A2div3, &threeInv)
	aW.Sub(&Bmg, &A2div3)

	// 8-torsion T8 (from Sage)
	XT8.A0.SetString("56495978839162271580189242121432892911")
	XT8.A1.SetString("5498965435856981984649943052775033891")
	YT8.A0.SetString("110406726416374112713382695534503225009")
	YT8.A1.SetString("41382232403199127624659984993520567989")
	lamT8.A0.SetString("144442680182921669228377474261162789491")
	lamT8.A1.SetString("78866152311824694793607513268751646905")

	XT4.A0.SetString("170141183460469230329734754113958182802")
	XT4.A1.SetString("71655106159052621705899442625265968763")
	YT4.A0.SetString("80492913427091964959665255396056504603")
	YT4.A1.SetString("170141183460469223319972006104328568185")
	lamT4.A0.SetZero()
	lamT4.A1.SetString("2")

	XT2.A0.SetString("2803905099203851845846")
	XT2.A1.SetString("26830971142363988319888418465352168201")

	// 7-torsion S7 and Miller loop intermediates (from Sage)
	XS7.A0.SetString("170000906615246946031520731709186440780")
	XS7.A1.SetString("94692800017046438403892332000427881085")
	YS7.A0.SetString("113074451165823784128958858701762803038")
	YS7.A1.SetString("4234707370120948597958144423972584168")

	lamDbl1.A0.SetString("149197531863252711357534985748910785211")
	lamDbl1.A1.SetString("106728816833063046549591444245799095908")

	X2S.A0.SetString("54301625611251160727403361627257607395")
	X2S.A1.SetString("98913730175858778379537234706939518181")
	Y2S.A0.SetString("64618348732980520867047493319367212765")
	Y2S.A1.SetString("143042194917927051863741747443651007549")

	lamAdd1.A0.SetString("119059155663838527172282288165420875213")
	lamAdd1.A1.SetString("2641556124042940121879682364519729700")

	X3S.A0.SetString("156066723560704030263055387453284139745")
	X3S.A1.SetString("27481310916039459192906723051714574829")
	Y3S.A0.SetString("61002583797480505059569858362795384815")
	Y3S.A1.SetString("161099302072925642522021123389456975814")

	lamDbl2.A0.SetString("58231718181463545180414387022826039618")
	lamDbl2.A1.SetString("21100793565021504938029727675157445515")

	X6S.Set(&XS7) // [6]S7 = -S7, same X
}

// septicCheck performs the septic residuosity check χ₇(f7) = 1 using
// Norm accumulation: Norm(f7) = ∏ Norm(ℓᵢ)^eᵢ / ∏ Norm(vⱼ)^eⱼ, computed
// entirely in Fp without Fp2 inversions.
//
// f7 = (ℓD1/vD1 · ℓA1/vA1)² · ℓD2/vD2 · gVert
// Norm(f7) = Norm(ℓD1)²·Norm(ℓA1)²·Norm(ℓD2)·Norm(gVert) / (Norm(vD1)²·Norm(vA1)²·Norm(vD2))
func septicCheck(XQ, YQ *fp2.E2) bool {
	// Step 1a: double S7 → [2]S7
	var ellD1, vD1 fp2.E2
	{
		var tmp fp2.E2
		tmp.Sub(XQ, &XS7)
		ellD1.Mul(&lamDbl1, &tmp)
		ellD1.Sub(YQ, &ellD1)
		ellD1.Sub(&ellD1, &YS7)
		vD1.Sub(XQ, &X2S)
	}

	// Step 1b: add [2]S7 + S7 → [3]S7
	var ellA1, vA1 fp2.E2
	{
		var tmp fp2.E2
		tmp.Sub(XQ, &X2S)
		ellA1.Mul(&lamAdd1, &tmp)
		ellA1.Sub(YQ, &ellA1)
		ellA1.Sub(&ellA1, &Y2S)
		vA1.Sub(XQ, &X3S)
	}

	// Step 2a: double [3]S7 → [6]S7
	var ellD2, vD2 fp2.E2
	{
		var tmp fp2.E2
		tmp.Sub(XQ, &X3S)
		ellD2.Mul(&lamDbl2, &tmp)
		ellD2.Sub(YQ, &ellD2)
		ellD2.Sub(&ellD2, &Y3S)
		vD2.Sub(XQ, &X6S)
	}

	// Step 2b: vertical line
	var gVert fp2.E2
	gVert.Sub(XQ, &XS7)

	// Compute Norm of each line/vertical in Fp (2S + 1A each).
	nEllD1 := ellD1.Norm()
	nVD1 := vD1.Norm()
	nEllA1 := ellA1.Norm()
	nVA1 := vA1.Norm()
	nEllD2 := ellD2.Norm()
	nVD2 := vD2.Norm()
	nGVert := gVert.Norm()

	// Numerator: Norm(ℓD1)² · Norm(ℓA1)² · Norm(ℓD2) · Norm(gVert)
	var num fp.Element
	num.Square(&nEllD1) // Norm(ℓD1)²
	var t fp.Element
	t.Square(&nEllA1) // Norm(ℓA1)²
	num.Mul(&num, &t)
	num.Mul(&num, &nEllD2)
	num.Mul(&num, &nGVert)

	// Denominator: Norm(vD1)² · Norm(vA1)² · Norm(vD2)
	var den fp.Element
	den.Square(&nVD1) // Norm(vD1)²
	t.Square(&nVA1)   // Norm(vA1)²
	den.Mul(&den, &t)
	den.Mul(&den, &nVD2)

	// Norm(f7) = num / den
	den.Inverse(&den)
	num.Mul(&num, &den)

	// χ₇(f7) = Norm(f7)^((p-1)/7) ?= 1
	var result fp.Element
	result.ExpBySepticFp(num)
	return result.IsOne()
}

func edwardsToWeierstrass(p *PointAffine) (X, Y fp2.E2) {
	var one fp2.E2
	one.SetOne()
	var onePlusY, oneMinusY, prod, inv, invOMY, invX fp2.E2
	onePlusY.Add(&one, &p.Y)
	oneMinusY.Sub(&one, &p.Y)
	prod.Mul(&p.X, &oneMinusY)
	inv.Inverse(&prod)
	invOMY.Mul(&inv, &p.X)
	invX.Mul(&inv, &oneMinusY)

	amd := new(fp2.E2).Sub(&curveParams.A, &curveParams.D)
	var u, w fp2.E2
	u.Mul(amd, &onePlusY)
	u.Mul(&u, &invOMY)
	w.Add(&invX, &invX)
	X.Add(&u, &Adiv3)
	Y.Mul(&u, &w)
	return
}

// isInSubGroupNaive tests subgroup membership by scalar multiplication by N.
func (p *PointAffine) isInSubGroupNaive() bool {
	initOnce.Do(initCurveParams)
	var res PointAffine
	res.ScalarMultiplication(p, &curveParams.Order)
	return res.IsZero()
}

// ClearCofactor sets p to [392]*p and returns p.
func (p *PointAffine) ClearCofactor() *PointAffine {
	p.ScalarMultiplication(p, big.NewInt(392))
	return p
}

// isInSubGroupTate tests subgroup membership using the divide-and-pair
// method: 0 divisions + octic check (2-part) + septic check (7-part).
//
// Octic: degree-8 Miller function on Weierstrass via 3-step doubling chain,
// then Frobenius trick: z^(p-1) = conj(z)/z followed by 124 squarings
// (since (p+1)/8 = 2^124 for the Mersenne prime p = 2^127-1).
//
// Septic: degree-7 Miller function via precomputed 2-step Miller loop
// (7 = 111₂), then z^((p-1)/7) and Norm check.
func (p *PointAffine) isInSubGroupTate() bool {
	smtInitOnce.Do(initSMTConstants)
	initOnce.Do(initCurveParams)

	if p.IsZero() {
		return true
	}

	XQ, YQ := edwardsToWeierstrass(p)

	// === Octic check (torus/trace approach) ===
	// f8 = ell1^4 * v1^4 * ell2^2 * v2^7 (division-free)
	var xd1, ell1, v1 fp2.E2
	xd1.Sub(&XQ, &XT8)
	ell1.Mul(&lamT8, &xd1)
	ell1.Sub(&YQ, &ell1)
	ell1.Sub(&ell1, &YT8)
	v1.Sub(&XQ, &XT4)

	var xd2, ell2, v2 fp2.E2
	xd2.Sub(&XQ, &XT4)
	ell2.Mul(&lamT4, &xd2)
	ell2.Sub(&YQ, &ell2)
	ell2.Sub(&ell2, &YT4)
	v2.Sub(&XQ, &XT2)

	var ell14, v14, ell22, v27, f8 fp2.E2
	ell14.Square(&ell1)
	ell14.Square(&ell14)
	v14.Square(&v1)
	v14.Square(&v14)
	ell22.Square(&ell2)
	v27.Square(&v2)
	v27.Mul(&v27, &v2)
	v27.Square(&v27)
	v27.Mul(&v27, &v2)

	f8.Mul(&ell14, &v14)
	f8.Mul(&f8, &ell22)
	f8.Mul(&f8, &v27)

	// g = conj(f8)/f8. Trace: t = g + g⁻¹ = 2(a²-b²)/(a²+b²) for f8 = a+bi.
	// Compute projectively: T = 2(a²-b²), N = a²+b² = Norm(f8).
	var a2, b2 fp.Element
	a2.Square(&f8.A0)
	b2.Square(&f8.A1)
	var T, N fp.Element
	T.Sub(&a2, &b2)
	T.Double(&T)    // T = 2(a²-b²)
	N.Add(&a2, &b2) // N = Norm(f8)

	// If f8 = 0 (degenerate Miller value), the point is not in the subgroup.
	if N.IsZero() {
		return false
	}

	// 124 iterations of t → t²−2 in projective form:
	// (T, N) → (T²−2N², N²)
	var T2, N2, twoN2 fp.Element
	for range 124 {
		T2.Square(&T)
		N2.Square(&N)
		twoN2.Double(&N2)
		T.Sub(&T2, &twoN2)
		N.Set(&N2)
	}

	// Check t == 2, i.e., T == 2·N
	var twoN fp.Element
	twoN.Double(&N)
	if !T.Equal(&twoN) {
		return false
	}

	return septicCheck(&XQ, &YQ)
}

// IsInSubGroup tests subgroup membership using the fastest available method.
func (p *PointAffine) IsInSubGroup() bool {
	return p.isInSubGroupTate()
}
