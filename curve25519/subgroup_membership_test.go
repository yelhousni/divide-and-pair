package curve25519

import (
	"crypto/rand"
	"math/big"
	"testing"

	"github.com/leanovate/gopter"
	"github.com/leanovate/gopter/prop"
	fp "github.com/yelhousni/divide-and-pair/curve25519/fp"
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
			quarticExp := p.isInSubGroupQuarticExp()

			return naive && pornin && quartic && quarticExp
		},
		genS,
	))

	properties.Property("all methods reject non-subgroup points", prop.ForAll(
		func(s big.Int) bool {
			params := curveParameters()
			if s.Sign() == 0 {
				return true // skip zero
			}
			// k*Base + (0,-1) is not in subgroup
			var nPt, p, q PointAffine
			nPt.X.SetZero()
			nPt.Y.SetOne()
			nPt.Y.Neg(&nPt.Y)

			p.ScalarMultiplication(&params.Base, &s)
			q.Add(&p, &nPt)

			naive := q.isInSubGroupNaive()
			pornin := q.isInSubGroupPornin()
			quartic := q.isInSubGroupQuartic()
			quarticExp := q.isInSubGroupQuarticExp()

			return !naive && !pornin && !quartic && !quarticExp
		},
		genS,
	))

	properties.TestingRun(t, gopter.ConsoleReporter(false))
}

func TestQuarticSymbolFromSubgroup(t *testing.T) {
	var one fp.Element
	one.SetOne()
	if one.QuarticSymbol() != 0 {
		t.Fatal("quarticSymbol(1) should be 0 (=1)")
	}

	var a, b, ab fp.Element
	a.SetUint64(3)
	b.SetUint64(7)
	ab.Mul(&a, &b)
	ca := a.QuarticSymbol()
	cb := b.QuarticSymbol()
	cab := ab.QuarticSymbol()
	if cab != (ca+cb)%4 {
		t.Fatalf("quarticSymbol not multiplicative: χ₄(%d)=%d, χ₄(%d)=%d, χ₄(%d·%d)=%d",
			3, ca, 7, cb, 3, 7, cab)
	}
}

func TestLowOrderPoints(t *testing.T) {
	subgroupInitOnce.Do(initSubgroupConstants)

	// Identity (0, 1) - in subgroup
	var id PointAffine
	id.X.SetZero()
	id.Y.SetOne()
	if !id.isInSubGroupPornin() {
		t.Fatal("identity should be in subgroup (Pornin)")
	}
	if !id.isInSubGroupQuartic() {
		t.Fatal("identity should be in subgroup (quartic)")
	}

	// Order-2 point: (0, -1) - NOT in subgroup
	var n PointAffine
	n.X.SetZero()
	n.Y.SetOne()
	n.Y.Neg(&n.Y)
	if n.isInSubGroupPornin() {
		t.Fatal("(0,-1) should NOT be in subgroup (Pornin)")
	}
	if n.isInSubGroupQuartic() {
		t.Fatal("(0,-1) should NOT be in subgroup (quartic)")
	}

	// 4-torsion: (i, 0)
	var t4 PointAffine
	t4.X.Set(&sqrtMinOne)
	t4.Y.SetZero()
	if t4.IsOnCurve() {
		if t4.isInSubGroupPornin() {
			t.Fatal("(i,0) should NOT be in subgroup (Pornin)")
		}
		if t4.isInSubGroupQuartic() {
			t.Fatal("(i,0) should NOT be in subgroup (quartic)")
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

func BenchmarkIsInSubGroupQuarticExp(b *testing.B) {
	benchSubgroup(b, (*PointAffine).isInSubGroupQuarticExp)
}

func BenchmarkIsInSubGroupQuartic(b *testing.B) {
	benchSubgroup(b, (*PointAffine).isInSubGroupQuartic)
}
