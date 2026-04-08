package gc256a

import (
	"crypto/rand"
	"math/big"
	"testing"

	"github.com/leanovate/gopter"
	"github.com/leanovate/gopter/prop"
)

func TestCurveParams(t *testing.T) {
	params := curveParameters()

	if !params.Base.IsOnCurve() {
		t.Fatal("base point not on curve")
	}

	var res PointAffine
	res.ScalarMultiplication(&params.Base, &params.Order)
	if !res.IsZero() {
		t.Fatal("ell * Base != O")
	}

	n := new(big.Int).Mul(&params.Order, big.NewInt(4))
	res.ScalarMultiplication(&params.Base, n)
	if !res.IsZero() {
		t.Fatal("4*ell * Base != O")
	}
}

func TestSubgroupMembership(t *testing.T) {
	t.Parallel()
	parameters := gopter.DefaultTestParameters()
	parameters.MinSuccessfulTests = 50

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

	properties.Property("[k]G is in subgroup (Pornin)", prop.ForAll(
		func(s big.Int) bool {
			params := curveParameters()
			var p PointAffine
			p.ScalarMultiplication(&params.Base, &s)
			return p.isInSubGroupPornin()
		},
		genS,
	))

	properties.Property("[k]G is in subgroup (quartic)", prop.ForAll(
		func(s big.Int) bool {
			params := curveParameters()
			var p PointAffine
			p.ScalarMultiplication(&params.Base, &s)
			return p.isInSubGroupQuartic()
		},
		genS,
	))

	properties.TestingRun(t, gopter.ConsoleReporter(false))
}

func TestSubgroupAgreement(t *testing.T) {
	t.Parallel()
	parameters := gopter.DefaultTestParameters()
	parameters.MinSuccessfulTests = 50

	properties := gopter.NewProperties(parameters)
	genS := GenBigInt()

	properties.Property("all methods agree on subgroup points", prop.ForAll(
		func(s big.Int) bool {
			params := curveParameters()
			var p PointAffine
			p.ScalarMultiplication(&params.Base, &s)

			naive := p.isInSubGroupNaive()
			pornin := p.isInSubGroupPornin()
			quartic := p.isInSubGroupQuartic()

			return naive && pornin && quartic
		},
		genS,
	))

	properties.Property("all methods reject order-2 offset", prop.ForAll(
		func(s big.Int) bool {
			params := curveParameters()
			if s.Sign() == 0 {
				return true
			}
			var nPt, p, q PointAffine
			nPt.X.SetZero()
			nPt.Y.SetOne()
			nPt.Y.Neg(&nPt.Y)

			p.ScalarMultiplication(&params.Base, &s)
			q.Add(&p, &nPt)

			return !q.isInSubGroupNaive() && !q.isInSubGroupPornin() && !q.isInSubGroupQuartic()
		},
		genS,
	))

	properties.Property("all methods reject order-4 offset", prop.ForAll(
		func(s big.Int) bool {
			params := curveParameters()
			if s.Sign() == 0 {
				return true
			}
			var order4Pt, p, q PointAffine
			order4Pt.X.SetOne()
			order4Pt.Y.SetZero()

			p.ScalarMultiplication(&params.Base, &s)
			q.Add(&p, &order4Pt)

			return !q.isInSubGroupPornin() && !q.isInSubGroupQuartic()
		},
		genS,
	))

	properties.TestingRun(t, gopter.ConsoleReporter(false))
}

func TestLowOrderPoints(t *testing.T) {
	subgroupInitOnce.Do(initSubgroupConstants)

	var id PointAffine
	id.X.SetZero()
	id.Y.SetOne()
	if !id.isInSubGroupPornin() {
		t.Fatal("identity should be in subgroup (Pornin)")
	}

	var n PointAffine
	n.X.SetZero()
	n.Y.SetOne()
	n.Y.Neg(&n.Y)
	if n.isInSubGroupPornin() {
		t.Fatal("(0,-1) should NOT be in subgroup (Pornin)")
	}

	var t4 PointAffine
	t4.X.SetOne()
	t4.Y.SetZero()
	if t4.IsOnCurve() {
		if t4.isInSubGroupPornin() {
			t.Fatal("(1,0) should NOT be in subgroup (Pornin)")
		}
	}

	var t4neg PointAffine
	t4neg.X.SetOne()
	t4neg.X.Neg(&t4neg.X)
	t4neg.Y.SetZero()
	if t4neg.IsOnCurve() {
		if t4neg.isInSubGroupPornin() {
			t.Fatal("(-1,0) should NOT be in subgroup (Pornin)")
		}
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

func BenchmarkIsInSubGroupPornin(b *testing.B) {
	benchSubgroup(b, (*PointAffine).isInSubGroupPornin)
}

func BenchmarkIsInSubGroupQuartic(b *testing.B) {
	benchSubgroup(b, (*PointAffine).isInSubGroupQuartic)
}
