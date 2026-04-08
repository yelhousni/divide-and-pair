package fourq

import (
	"crypto/rand"
	"testing"

	"github.com/leanovate/gopter"
	"github.com/leanovate/gopter/prop"
	"math/big"
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
