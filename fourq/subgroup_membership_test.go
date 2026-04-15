package fourq

import (
	"crypto/rand"
	"math/big"
	"testing"

	"github.com/leanovate/gopter"
	"github.com/leanovate/gopter/prop"
	fp "github.com/yelhousni/divide-and-pair/fourq/fp"
	fp2 "github.com/yelhousni/divide-and-pair/fourq/fp2"
)

type legacySepticData struct {
	xS, yS   fp2.E2
	lamDbl1  fp2.E2
	x2S, y2S fp2.E2
	lamAdd1  fp2.E2
	x3S, y3S fp2.E2
	lamDbl2  fp2.E2
	x6S      fp2.E2
}

var (
	legacyS7 = legacySepticData{
		xS:      fp2.E2{A0: fp.Element{}, A1: fp.Element{}},
		yS:      fp2.E2{A0: fp.Element{}, A1: fp.Element{}},
		lamDbl1: fp2.E2{A0: fp.Element{}, A1: fp.Element{}},
		x2S:     fp2.E2{A0: fp.Element{}, A1: fp.Element{}},
		y2S:     fp2.E2{A0: fp.Element{}, A1: fp.Element{}},
		lamAdd1: fp2.E2{A0: fp.Element{}, A1: fp.Element{}},
		x3S:     fp2.E2{A0: fp.Element{}, A1: fp.Element{}},
		y3S:     fp2.E2{A0: fp.Element{}, A1: fp.Element{}},
		lamDbl2: fp2.E2{A0: fp.Element{}, A1: fp.Element{}},
		x6S:     fp2.E2{A0: fp.Element{}, A1: fp.Element{}},
	}
	legacyS7p = legacySepticData{
		xS:      fp2.E2{A0: fp.Element{}, A1: fp.Element{}},
		yS:      fp2.E2{A0: fp.Element{}, A1: fp.Element{}},
		lamDbl1: fp2.E2{A0: fp.Element{}, A1: fp.Element{}},
		x2S:     fp2.E2{A0: fp.Element{}, A1: fp.Element{}},
		y2S:     fp2.E2{A0: fp.Element{}, A1: fp.Element{}},
		lamAdd1: fp2.E2{A0: fp.Element{}, A1: fp.Element{}},
		x3S:     fp2.E2{A0: fp.Element{}, A1: fp.Element{}},
		y3S:     fp2.E2{A0: fp.Element{}, A1: fp.Element{}},
		lamDbl2: fp2.E2{A0: fp.Element{}, A1: fp.Element{}},
		x6S:     fp2.E2{A0: fp.Element{}, A1: fp.Element{}},
	}
)

func init() {
	legacyS7.xS.A0.SetString("170000906615246946031520731709186440780")
	legacyS7.xS.A1.SetString("94692800017046438403892332000427881085")
	legacyS7.yS.A0.SetString("113074451165823784128958858701762803038")
	legacyS7.yS.A1.SetString("4234707370120948597958144423972584168")
	legacyS7.lamDbl1.A0.SetString("149197531863252711357534985748910785211")
	legacyS7.lamDbl1.A1.SetString("106728816833063046549591444245799095908")
	legacyS7.x2S.A0.SetString("54301625611251160727403361627257607395")
	legacyS7.x2S.A1.SetString("98913730175858778379537234706939518181")
	legacyS7.y2S.A0.SetString("64618348732980520867047493319367212765")
	legacyS7.y2S.A1.SetString("143042194917927051863741747443651007549")
	legacyS7.lamAdd1.A0.SetString("119059155663838527172282288165420875213")
	legacyS7.lamAdd1.A1.SetString("2641556124042940121879682364519729700")
	legacyS7.x3S.A0.SetString("156066723560704030263055387453284139745")
	legacyS7.x3S.A1.SetString("27481310916039459192906723051714574829")
	legacyS7.y3S.A0.SetString("61002583797480505059569858362795384815")
	legacyS7.y3S.A1.SetString("161099302072925642522021123389456975814")
	legacyS7.lamDbl2.A0.SetString("58231718181463545180414387022826039618")
	legacyS7.lamDbl2.A1.SetString("21100793565021504938029727675157445515")
	legacyS7.x6S.Set(&legacyS7.xS)

	legacyS7p.xS.A0.SetString("134226817183232424897437757433869201092")
	legacyS7p.xS.A1.SetString("130939014160170965664116628940157096270")
	legacyS7p.yS.A0.SetString("44780681017233260662086977535495331704")
	legacyS7p.yS.A1.SetString("60082122951898920823140075210382386075")
	legacyS7p.lamDbl1.A0.SetString("79876232603258742370078483576637952910")
	legacyS7p.lamDbl1.A1.SetString("93431063624494182359521542759712946209")
	legacyS7p.x2S.A0.SetString("14985878640945021751153055736915441094")
	legacyS7p.x2S.A1.SetString("35562414809826163017035549423954142235")
	legacyS7p.y2S.A0.SetString("20351043536092277668563817211857088150")
	legacyS7p.y2S.A1.SetString("88111568418789077182609138749030999981")
	legacyS7p.lamAdd1.A0.SetString("135655560401356473744585404315093205196")
	legacyS7p.lamAdd1.A1.SetString("158246920428350027040846798030706302090")
	legacyS7p.x3S.A0.SetString("86614998441207314901553846558717181209")
	legacyS7p.x3S.A1.SetString("119520286240226947588091143044800503760")
	legacyS7p.y3S.A0.SetString("126913250202842280960467879346012678119")
	legacyS7p.y3S.A1.SetString("38938712949820441804460416342134306702")
	legacyS7p.lamDbl2.A0.SetString("149509976829250023099150978250294296208")
	legacyS7p.lamDbl2.A1.SetString("150620094415765623422294113057405903169")
	legacyS7p.x6S.Set(&legacyS7p.xS)
}

func legacySepticCheckWith(XQ, YQ *fp2.E2, data *legacySepticData) bool {
	var diffXS, diffX2S, diffX3S fp2.E2
	diffXS.Sub(XQ, &data.xS)
	diffX2S.Sub(XQ, &data.x2S)
	diffX3S.Sub(XQ, &data.x3S)

	var ellD1, vD1 fp2.E2
	ellD1.Mul(&data.lamDbl1, &diffXS)
	ellD1.Sub(YQ, &ellD1)
	ellD1.Sub(&ellD1, &data.yS)
	vD1.Set(&diffX2S)

	var ellA1, vA1 fp2.E2
	ellA1.Mul(&data.lamAdd1, &diffX2S)
	ellA1.Sub(YQ, &ellA1)
	ellA1.Sub(&ellA1, &data.y2S)
	vA1.Set(&diffX3S)

	var ellD2, vD2 fp2.E2
	ellD2.Mul(&data.lamDbl2, &diffX3S)
	ellD2.Sub(YQ, &ellD2)
	ellD2.Sub(&ellD2, &data.y3S)
	vD2.Sub(XQ, &data.x6S)

	nEllD1 := ellD1.Norm()
	nVD1 := vD1.Norm()
	nEllA1 := ellA1.Norm()
	nVA1 := vA1.Norm()
	nEllD2 := ellD2.Norm()
	nVD2 := vD2.Norm()
	nGVert := diffXS.Norm()

	var num fp.Element
	num.Square(&nEllD1)
	var t fp.Element
	t.Square(&nEllA1)
	num.Mul(&num, &t)
	num.Mul(&num, &nEllD2)
	num.Mul(&num, &nGVert)

	var den fp.Element
	den.Square(&nVD1)
	t.Square(&nVA1)
	den.Mul(&den, &t)
	den.Mul(&den, &nVD2)

	den.Inverse(&den)
	num.Mul(&num, &den)

	var result fp.Element
	result.ExpBySepticFp(num)
	return result.IsOne()
}

func doubleOnWeierstrass(X, Y *fp2.E2) (fp2.E2, fp2.E2, fp2.E2) {
	var lam, x2, threeX2, twoY fp2.E2
	x2.Square(X)
	threeX2.Add(&x2, &x2)
	threeX2.Add(&threeX2, &x2)
	threeX2.Add(&threeX2, &aW)
	twoY.Double(Y)
	lam.Inverse(&twoY)
	lam.Mul(&lam, &threeX2)

	var X2, Y2 fp2.E2
	X2.Square(&lam)
	X2.Sub(&X2, X)
	X2.Sub(&X2, X)
	Y2.Sub(X, &X2)
	Y2.Mul(&lam, &Y2)
	Y2.Sub(&Y2, Y)
	return X2, Y2, lam
}

func TestSubgroupMembership(t *testing.T) {
	t.Parallel()
	parameters := gopter.DefaultTestParameters()
	parameters.MinSuccessfulTests = 20

	properties := gopter.NewProperties(parameters)
	genS := GenBigInt()

	properties.Property("[k]G is in subgroup (naive)", prop.ForAll(
		func(s big.Int) bool {
			params := curveParameters()
			var p PointAffine
			p.ScalarMultiplication(&params.Base, &s)
			return p.isInSubGroupNaive()
		},
		genS,
	))

	properties.Property("[k]G is in subgroup (Tate)", prop.ForAll(
		func(s big.Int) bool {
			params := curveParameters()
			var p PointAffine
			p.ScalarMultiplication(&params.Base, &s)
			return p.isInSubGroupTate()
		},
		genS,
	))

	properties.TestingRun(t, gopter.ConsoleReporter(false))
}

func TestSubgroupAgreement(t *testing.T) {
	t.Parallel()
	parameters := gopter.DefaultTestParameters()
	parameters.MinSuccessfulTests = 20

	properties := gopter.NewProperties(parameters)
	genS := GenBigInt()

	properties.Property("naive and Tate agree on subgroup points", prop.ForAll(
		func(s big.Int) bool {
			params := curveParameters()
			var p PointAffine
			p.ScalarMultiplication(&params.Base, &s)

			naive := p.isInSubGroupNaive()
			tate := p.isInSubGroupTate()

			return naive && tate
		},
		genS,
	))

	properties.TestingRun(t, gopter.ConsoleReporter(false))
}

func TestSepticRejection(t *testing.T) {
	smtInitOnce.Do(initSMTConstants)
	initOnce.Do(initCurveParams)
	params := curveParameters()

	// Construct a 7-torsion point on Edwards from S7' on Weierstrass.
	// Weierstrass: X = u + A/3, Y = u*w => u = X - A/3, w = Y/u
	// Edwards from Montgomery: y = (u_m - 1)/(u_m + 1), x = 2*u_m/(u_m*w)
	// where u_m = u_pornin / (a-d), w_pornin = w
	// Here u_pornin = u (from the Weierstrass inverse), u = X - Adiv3
	// Actually: X = u_pornin + Adiv3, Y = u_pornin * w_pornin
	// So u_pornin = X - Adiv3, w_pornin = Y / u_pornin
	// Then u_montgomery = u_pornin / (a-d)
	// Edwards: y = (u_m - 1)/(u_m + 1), x = u_m * (a-d) * w_pornin ... no, x = 2/w_pornin
	// From the birational: u_pornin = (a-d)*(1+y)/(1-y), w_pornin = 2/x
	// Inverse: (1-y)*u_pornin = (a-d)*(1+y) => u_pornin - u_pornin*y = (a-d) + (a-d)*y
	// => y*(u_pornin + a-d) = u_pornin - (a-d) => y = (u_pornin - (a-d))/(u_pornin + (a-d))
	// x = 2/w_pornin

	amd := new(fp2.E2).Sub(&curveParams.A, &curveParams.D)

	var u_pornin, w_pornin fp2.E2
	u_pornin.Sub(&XS7p, &Adiv3)
	// w_pornin = Y / u_pornin
	var u_inv fp2.E2
	u_inv.Inverse(&u_pornin)
	w_pornin.Mul(&YS7p, &u_inv)

	// y = (u_pornin - (a-d)) / (u_pornin + (a-d))
	var num, den fp2.E2
	num.Sub(&u_pornin, amd)
	den.Add(&u_pornin, amd)
	den.Inverse(&den)
	var yT7 fp2.E2
	yT7.Mul(&num, &den)

	// x = 2 / w_pornin
	var xT7 fp2.E2
	var two fp2.E2
	two.A0.SetUint64(2)
	w_pornin.Inverse(&w_pornin)
	xT7.Mul(&two, &w_pornin)

	var T7 PointAffine
	T7.X.Set(&xT7)
	T7.Y.Set(&yT7)

	if !T7.IsOnCurve() {
		t.Fatal("T7 not on curve")
	}

	// Verify [7]*T7 = O
	var check PointAffine
	check.ScalarMultiplication(&T7, big.NewInt(7))
	if !check.IsZero() {
		t.Fatal("[7]*T7 should be zero")
	}

	// T7 should NOT be in the subgroup
	if T7.isInSubGroupNaive() {
		t.Fatal("7-torsion point should NOT be in subgroup (naive)")
	}
	if T7.isInSubGroupTate() {
		t.Fatal("7-torsion point should NOT be in subgroup (Tate)")
	}

	// Subgroup point + T7 should NOT be in the subgroup
	k, _ := rand.Int(rand.Reader, &params.Order)
	var P, Q PointAffine
	P.ScalarMultiplication(&params.Base, k)
	Q.Add(&P, &T7)

	if Q.isInSubGroupNaive() {
		t.Fatal("subgroup + 7-torsion should NOT be in subgroup (naive)")
	}
	if Q.isInSubGroupTate() {
		t.Fatal("subgroup + 7-torsion should NOT be in subgroup (Tate)")
	}
	t.Log("septic rejection test passed: both S7 and S7' checks needed")
}

func TestSepticMillerConstants(t *testing.T) {
	smtInitOnce.Do(initSMTConstants)

	cases := []struct {
		name         string
		x1, y1, lam1 *fp2.E2
		x2, y2, lam2 *fp2.E2
		x4, y4, lam4 *fp2.E2
	}{
		{"S7", &XS7, &YS7, &lamDbl1, &X2S, &Y2S, &lamDbl2, &X4S, &Y4S, &lamDbl4},
		{"S7'", &XS7p, &YS7p, &lamDbl1p, &X2Sp, &Y2Sp, &lamDbl2p, &X4Sp, &Y4Sp, &lamDbl4p},
	}

	for _, tc := range cases {
		x2, y2, lam1 := doubleOnWeierstrass(tc.x1, tc.y1)
		if !x2.Equal(tc.x2) || !y2.Equal(tc.y2) || !lam1.Equal(tc.lam1) {
			t.Fatalf("%s: first doubling constants mismatch", tc.name)
		}

		x4, y4, lam2 := doubleOnWeierstrass(tc.x2, tc.y2)
		if !x4.Equal(tc.x4) || !y4.Equal(tc.y4) || !lam2.Equal(tc.lam2) {
			t.Fatalf("%s: second doubling constants mismatch", tc.name)
		}

		x8, y8, lam4 := doubleOnWeierstrass(tc.x4, tc.y4)
		if !x8.Equal(tc.x1) || !y8.Equal(tc.y1) || !lam4.Equal(tc.lam4) {
			t.Fatalf("%s: expected [8]S = S in the stored f8 chain", tc.name)
		}
	}
}

func TestSepticF8MatchesLegacyF7(t *testing.T) {
	smtInitOnce.Do(initSMTConstants)
	initOnce.Do(initCurveParams)
	params := curveParameters()

	checkPoint := func(t *testing.T, p *PointAffine) {
		t.Helper()
		XQ, YQ := edwardsToWeierstrass(p)

		if legacySepticCheckWith(&XQ, &YQ, &legacyS7) != septicCheck(&XQ, &YQ) {
			t.Fatal("S7 legacy f7 check disagrees with f8 replacement")
		}
		if legacySepticCheckWith(&XQ, &YQ, &legacyS7p) != septicCheckPrime(&XQ, &YQ) {
			t.Fatal("S7' legacy f7 check disagrees with f8 replacement")
		}
	}

	// Identity and a few subgroup points.
	var id PointAffine
	id.SetInfinity()
	checkPoint(t, &id)

	for range 8 {
		k, err := rand.Int(rand.Reader, &params.Order)
		if err != nil {
			t.Fatal(err)
		}
		var p PointAffine
		p.ScalarMultiplication(&params.Base, k)
		checkPoint(t, &p)
	}

	// Explicit 7-torsion representative and subgroup + torsion points.
	amd := new(fp2.E2).Sub(&curveParams.A, &curveParams.D)
	var uPornin, wPornin fp2.E2
	uPornin.Sub(&XS7p, &Adiv3)
	var uInv fp2.E2
	uInv.Inverse(&uPornin)
	wPornin.Mul(&YS7p, &uInv)

	var num, den fp2.E2
	num.Sub(&uPornin, amd)
	den.Add(&uPornin, amd)
	den.Inverse(&den)
	var yT7 fp2.E2
	yT7.Mul(&num, &den)

	var xT7, two fp2.E2
	two.A0.SetUint64(2)
	wPornin.Inverse(&wPornin)
	xT7.Mul(&two, &wPornin)

	var t7 PointAffine
	t7.X.Set(&xT7)
	t7.Y.Set(&yT7)
	checkPoint(t, &t7)

	for range 4 {
		k, err := rand.Int(rand.Reader, &params.Order)
		if err != nil {
			t.Fatal(err)
		}
		var p, q PointAffine
		p.ScalarMultiplication(&params.Base, k)
		q.Add(&p, &t7)
		checkPoint(t, &q)
	}
}

func TestLowOrderPoints(t *testing.T) {
	// Identity
	var id PointAffine
	id.SetInfinity()
	if !id.isInSubGroupTate() {
		t.Fatal("identity should be in subgroup (Tate)")
	}

	// (0, -1)
	var n PointAffine
	n.SetInfinity()
	n.Y.Neg(&n.Y)
	if n.IsOnCurve() && n.isInSubGroupTate() {
		t.Fatal("(0,-1) should NOT be in subgroup (Tate)")
	}
}

// Benchmarks

func isInSubGroupTateLegacySeptic(p *PointAffine) bool {
	smtInitOnce.Do(initSMTConstants)
	initOnce.Do(initCurveParams)

	if p.IsZero() {
		return true
	}

	XQ, YQ := edwardsToWeierstrass(p)

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

	var a2, b2 fp.Element
	a2.Square(&f8.A0)
	b2.Square(&f8.A1)
	var T, N fp.Element
	T.Sub(&a2, &b2)
	T.Double(&T)
	N.Add(&a2, &b2)
	if N.IsZero() {
		return false
	}

	var T2, N2, twoN2 fp.Element
	for range 124 {
		T2.Square(&T)
		N2.Square(&N)
		twoN2.Double(&N2)
		T.Sub(&T2, &twoN2)
		N.Set(&N2)
	}

	var twoN fp.Element
	twoN.Double(&N)
	if !T.Equal(&twoN) {
		return false
	}

	return legacySepticCheckWith(&XQ, &YQ, &legacyS7) &&
		legacySepticCheckWith(&XQ, &YQ, &legacyS7p)
}

func benchSubgroup(b *testing.B, method func(*PointAffine) bool) {
	params := curveParameters()
	k, _ := rand.Int(rand.Reader, &params.Order)
	var p PointAffine
	p.ScalarMultiplication(&params.Base, k)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		method(&p)
	}
}

func benchSeptic(b *testing.B, method func(*fp2.E2, *fp2.E2) bool) {
	params := curveParameters()
	k, _ := rand.Int(rand.Reader, &params.Order)
	var p PointAffine
	p.ScalarMultiplication(&params.Base, k)
	XQ, YQ := edwardsToWeierstrass(&p)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		method(&XQ, &YQ)
	}
}

func BenchmarkIsInSubGroupNaive(b *testing.B) {
	benchSubgroup(b, (*PointAffine).isInSubGroupNaive)
}

func BenchmarkIsInSubGroupTate(b *testing.B) {
	benchSubgroup(b, (*PointAffine).isInSubGroupTate)
}

func BenchmarkIsInSubGroupTateLegacySeptic(b *testing.B) {
	benchSubgroup(b, isInSubGroupTateLegacySeptic)
}

func BenchmarkIsInSubGroupEndo(b *testing.B) {
	benchSubgroup(b, (*PointAffine).isInSubGroupEndo)
}

func BenchmarkSepticCheckF8(b *testing.B) {
	benchSeptic(b, func(XQ, YQ *fp2.E2) bool {
		return septicCheck(XQ, YQ) && septicCheckPrime(XQ, YQ)
	})
}

func BenchmarkSepticCheckLegacyF7(b *testing.B) {
	benchSeptic(b, func(XQ, YQ *fp2.E2) bool {
		return legacySepticCheckWith(XQ, YQ, &legacyS7) &&
			legacySepticCheckWith(XQ, YQ, &legacyS7p)
	})
}
