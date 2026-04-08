package jubjub

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

	n := new(big.Int).Mul(&params.Order, big.NewInt(8))
	res.ScalarMultiplication(&params.Base, n)
	if !res.IsZero() {
		t.Fatal("8*ell * Base != O")
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

	properties.Property("[k]G is in subgroup (octic exp)", prop.ForAll(
		func(s big.Int) bool {
			params := curveParameters()
			var p PointAffine
			p.ScalarMultiplication(&params.Base, &s)
			return p.isInSubGroupOcticExp()
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
			octic := p.isInSubGroupOcticExp()

			return naive && pornin && octic
		},
		genS,
	))

	properties.Property("all methods reject non-subgroup points", prop.ForAll(
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

			naive := q.isInSubGroupNaive()
			pornin := q.isInSubGroupPornin()
			octic := q.isInSubGroupOcticExp()

			return !naive && !pornin && !octic
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

	// 4-torsion: (i, 0)
	var t4 PointAffine
	t4.X.Set(&sqrtMinOne)
	t4.Y.SetZero()
	if t4.IsOnCurve() {
		if t4.isInSubGroupPornin() {
			t.Fatal("(i,0) should NOT be in subgroup (Pornin)")
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

func BenchmarkIsInSubGroupOcticExp(b *testing.B) {
	benchSubgroup(b, (*PointAffine).isInSubGroupOcticExp)
}
