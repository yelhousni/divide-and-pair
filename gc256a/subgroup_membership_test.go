package gc256a

import (
	"math/big"
	"math/rand"
	"testing"

	"github.com/yelhousni/divide-and-pair/gc256a/fp"
	"github.com/yelhousni/divide-and-pair/gc256a/fp2"

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

			return naive && pornin && quartic && p.isInSubGroupQuarticGCD()
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

			return !q.isInSubGroupNaive() && !q.isInSubGroupPornin() && !q.isInSubGroupQuartic() && !q.isInSubGroupQuarticGCD()
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

			return !q.isInSubGroupPornin() && !q.isInSubGroupQuartic() && !q.isInSubGroupQuarticGCD()
		},
		genS,
	))

	properties.TestingRun(t, gopter.ConsoleReporter(false))
}

func TestLowOrderPoints(t *testing.T) {
	subgroupInitOnce.Do(initSubgroupConstants)
	if isQuarticResidueTorus(fp2.E2{}) {
		t.Fatal("zero Miller value passed quartic test")
	}

	var id PointAffine
	id.X.SetZero()
	id.Y.SetOne()
	if !id.isInSubGroupPornin() || !id.IsInSubGroup() {
		t.Fatal("identity should be in subgroup")
	}

	var n PointAffine
	n.X.SetZero()
	n.Y.SetOne()
	n.Y.Neg(&n.Y)
	if n.isInSubGroupPornin() || n.IsInSubGroup() {
		t.Fatal("(0,-1) should NOT be in subgroup")
	}

	var t4 PointAffine
	t4.X.SetOne()
	t4.Y.SetZero()
	if t4.IsOnCurve() {
		if t4.isInSubGroupPornin() || t4.IsInSubGroup() {
			t.Fatal("(1,0) should NOT be in subgroup")
		}
	}

	var t4neg PointAffine
	t4neg.X.SetOne()
	t4neg.X.Neg(&t4neg.X)
	t4neg.Y.SetZero()
	if t4neg.IsOnCurve() {
		if t4neg.isInSubGroupPornin() || t4neg.IsInSubGroup() {
			t.Fatal("(-1,0) should NOT be in subgroup")
		}
	}
}

// Benchmarks

func benchSubgroup(b *testing.B, method func(*PointAffine) bool) {
	params := curveParameters()
	// Cycle over many random points so that input-dependent branches (e.g.
	// the extra square root taken in point halving when an intermediate is a
	// non-residue) are averaged over, rather than fixed by a single draw.
	const nPoints = 256
	points := make([]PointAffine, nPoints)
	rng := rand.New(rand.NewSource(1))
	for i := range points {
		k := new(big.Int).Rand(rng, &params.Order)
		points[i].ScalarMultiplication(&params.Base, k)
	}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		method(&points[i%nPoints])
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

func BenchmarkIsInSubGroupQuarticGCD(b *testing.B) {
	benchSubgroup(b, (*PointAffine).isInSubGroupQuarticGCD)
}

func BenchmarkQuarticExtension(b *testing.B) {
	subgroupInitOnce.Do(initSubgroupConstants)
	var inputs [256]fp2.E2
	rng := rand.New(rand.NewSource(1))
	p := fp.Modulus()
	for i := range inputs {
		inputs[i].A0.SetBigInt(new(big.Int).Rand(rng, p))
		inputs[i].A1.SetBigInt(new(big.Int).Rand(rng, p))
		if isQuarticResidueTorus(inputs[i]) != (inputs[i].QuarticSymbol() == 0) {
			b.Fatal("quartic predicates disagree")
		}
	}
	b.Run("TorusPredicate", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			isQuarticResidueTorus(inputs[i%len(inputs)])
		}
	})
	b.Run("GCDIndex", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			inputs[i%len(inputs)].QuarticSymbol()
		}
	})
}
