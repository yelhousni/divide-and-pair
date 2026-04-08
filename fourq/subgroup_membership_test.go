package fourq

import (
	"crypto/rand"
	"math/big"
	"testing"

	"github.com/leanovate/gopter"
	"github.com/leanovate/gopter/prop"
	fp2 "github.com/yelhousni/divide-and-pair/fourq/fp2"
)

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

func BenchmarkIsInSubGroupNaive(b *testing.B) {
	benchSubgroup(b, (*PointAffine).isInSubGroupNaive)
}

func BenchmarkIsInSubGroupTate(b *testing.B) {
	benchSubgroup(b, (*PointAffine).isInSubGroupTate)
}

func BenchmarkIsInSubGroupEndo(b *testing.B) {
	benchSubgroup(b, (*PointAffine).isInSubGroupEndo)
}
