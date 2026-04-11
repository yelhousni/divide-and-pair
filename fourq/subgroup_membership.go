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
	aMinusD fp2.E2
	aW      fp2.E2
	Adiv3   fp2.E2

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

	// Second independent 7-torsion point S7' and its Miller loop intermediates.
	// E[7](Fp2) ≅ Z7 × Z7, so we need two independent points to check both factors.
	XS7p, YS7p fp2.E2
	lamDbl1p   fp2.E2 // tangent slope at S7'
	X2Sp, Y2Sp fp2.E2 // [2]S7'
	lamAdd1p   fp2.E2 // chord slope S7', [2]S7'
	X3Sp, Y3Sp fp2.E2 // [3]S7'
	lamDbl2p   fp2.E2 // tangent slope at [3]S7'
	X6Sp       fp2.E2 // [6]S7' = -S7' (X only)

	three fp2.E2
)

func initSMTConstants() {
	initOnce.Do(initCurveParams)

	aMinusD.Sub(&curveParams.A, &curveParams.D)
	apd := new(fp2.E2).Add(&curveParams.A, &curveParams.D)

	var Amg, Bmg fp2.E2
	Amg.Add(apd, apd)
	Bmg.Square(&aMinusD)
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

	// Second independent 7-torsion point S7' (from Sage, seed=1).
	// S7 and S7' together generate E[7](Fp2) ≅ Z7 × Z7.
	XS7p.A0.SetString("134226817183232424897437757433869201092")
	XS7p.A1.SetString("130939014160170965664116628940157096270")
	YS7p.A0.SetString("44780681017233260662086977535495331704")
	YS7p.A1.SetString("60082122951898920823140075210382386075")

	lamDbl1p.A0.SetString("79876232603258742370078483576637952910")
	lamDbl1p.A1.SetString("93431063624494182359521542759712946209")

	X2Sp.A0.SetString("14985878640945021751153055736915441094")
	X2Sp.A1.SetString("35562414809826163017035549423954142235")
	Y2Sp.A0.SetString("20351043536092277668563817211857088150")
	Y2Sp.A1.SetString("88111568418789077182609138749030999981")

	lamAdd1p.A0.SetString("135655560401356473744585404315093205196")
	lamAdd1p.A1.SetString("158246920428350027040846798030706302090")

	X3Sp.A0.SetString("86614998441207314901553846558717181209")
	X3Sp.A1.SetString("119520286240226947588091143044800503760")
	Y3Sp.A0.SetString("126913250202842280960467879346012678119")
	Y3Sp.A1.SetString("38938712949820441804460416342134306702")

	lamDbl2p.A0.SetString("149509976829250023099150978250294296208")
	lamDbl2p.A1.SetString("150620094415765623422294113057405903169")

	X6Sp.Set(&XS7p) // [6]S7' = -S7', same X
}

// septicCheckWith performs the septic residuosity check χ₇(f_{7,S}(Q)) = 1
// for a given 7-torsion point S and its precomputed Miller loop intermediates.
// Uses Norm accumulation entirely in Fp without Fp2 inversions.
func septicCheckWith(XQ, YQ *fp2.E2,
	xS7, yS7, ld1, x2S, y2S, la1, x3S, y3S, ld2, x6S *fp2.E2) bool {
	var diffXS, diffX2S, diffX3S fp2.E2
	diffXS.Sub(XQ, xS7)
	diffX2S.Sub(XQ, x2S)
	diffX3S.Sub(XQ, x3S)

	// Step 1a: double S → [2]S
	var ellD1, vD1 fp2.E2
	{
		ellD1.Mul(ld1, &diffXS)
		ellD1.Sub(YQ, &ellD1)
		ellD1.Sub(&ellD1, yS7)
		vD1.Set(&diffX2S)
	}

	// Step 1b: add [2]S + S → [3]S
	var ellA1, vA1 fp2.E2
	{
		ellA1.Mul(la1, &diffX2S)
		ellA1.Sub(YQ, &ellA1)
		ellA1.Sub(&ellA1, y2S)
		vA1.Set(&diffX3S)
	}

	// Step 2a: double [3]S → [6]S
	var ellD2, vD2 fp2.E2
	{
		ellD2.Mul(ld2, &diffX3S)
		ellD2.Sub(YQ, &ellD2)
		ellD2.Sub(&ellD2, y3S)
		vD2.Sub(XQ, x6S)
	}

	// Step 2b: vertical line
	gVert := diffXS

	// Compute Norm of each line/vertical in Fp.
	nEllD1 := ellD1.Norm()
	nVD1 := vD1.Norm()
	nEllA1 := ellA1.Norm()
	nVA1 := vA1.Norm()
	nEllD2 := ellD2.Norm()
	nVD2 := vD2.Norm()
	nGVert := gVert.Norm()

	// Numerator: Norm(ℓD1)² · Norm(ℓA1)² · Norm(ℓD2) · Norm(gVert)
	var num fp.Element
	num.Square(&nEllD1)
	var t fp.Element
	t.Square(&nEllA1)
	num.Mul(&num, &t)
	num.Mul(&num, &nEllD2)
	num.Mul(&num, &nGVert)

	// Denominator: Norm(vD1)² · Norm(vA1)² · Norm(vD2)
	var den fp.Element
	den.Square(&nVD1)
	t.Square(&nVA1)
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

// septicCheck performs the septic check with the first 7-torsion point S7.
func septicCheck(XQ, YQ *fp2.E2) bool {
	return septicCheckWith(XQ, YQ,
		&XS7, &YS7, &lamDbl1, &X2S, &Y2S, &lamAdd1, &X3S, &Y3S, &lamDbl2, &X6S)
}

// septicCheckPrime performs the septic check with the second 7-torsion point S7'.
func septicCheckPrime(XQ, YQ *fp2.E2) bool {
	return septicCheckWith(XQ, YQ,
		&XS7p, &YS7p, &lamDbl1p, &X2Sp, &Y2Sp, &lamAdd1p, &X3Sp, &Y3Sp, &lamDbl2p, &X6Sp)
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

	var u, w fp2.E2
	u.Mul(&aMinusD, &onePlusY)
	u.Mul(&u, &invOMY)
	w.Add(&invX, &invX)
	X.Add(&u, &Adiv3)
	Y.Mul(&u, &w)
	return
}

// isInSubGroupNaive tests subgroup membership by scalar multiplication by N.
func (p *PointAffine) isInSubGroupNaive() bool {
	initOnce.Do(initCurveParams)
	var proj PointProj
	proj.FromAffine(p)
	var res PointProj
	res.ScalarMultiplication(&proj, &curveParams.Order)
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

	// === Septic checks ===
	// E[7](Fp2) ≅ Z7 × Z7, so we need two independent 7-torsion points
	// to verify Q ∈ [7]E(Fp2) (i.e., the 7-part of the cofactor is cleared).
	return septicCheck(&XQ, &YQ) && septicCheckPrime(&XQ, &YQ)
}

// IsInSubGroup tests subgroup membership using the fastest available method.
func (p *PointAffine) IsInSubGroup() bool {
	return p.isInSubGroupTate()
}
