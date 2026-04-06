package fourq

import (
	"crypto/rand"
	"math/big"
	"testing"
)

func TestCurveParams(t *testing.T) {
	params := GetFourQCurve()

	if !params.Base.IsOnCurve() {
		t.Fatal("base point not on curve")
	}

	// N * Base should be identity
	var res PointAffine
	res.ScalarMultiplication(&params.Base, &params.Order)
	if !res.IsZero() {
		t.Fatal("N * Base != O")
	}

	// 392 * N * Base = O
	n := new(big.Int).Mul(&params.Order, big.NewInt(392))
	res.ScalarMultiplication(&params.Base, n)
	if !res.IsZero() {
		t.Fatal("392*N * Base != O")
	}
}

func TestIdentity(t *testing.T) {
	var id PointAffine
	id.X.SetZero()
	id.Y.SetOne()
	if !id.IsZero() {
		t.Fatal("(0,1) should be identity")
	}
	if !id.IsOnCurve() {
		t.Fatal("identity should be on curve")
	}
}

func TestAddIdentity(t *testing.T) {
	params := GetFourQCurve()
	var id, res PointAffine
	id.X.SetZero()
	id.Y.SetOne()
	res.Add(&params.Base, &id)
	if !res.Equal(&params.Base) {
		t.Fatal("P + O != P")
	}
}

func TestDoubleConsistency(t *testing.T) {
	params := GetFourQCurve()
	var d1, d2 PointAffine
	d1.Double(&params.Base)
	d2.Add(&params.Base, &params.Base)
	if !d1.Equal(&d2) {
		t.Fatal("Double(P) != P + P")
	}
}

func TestNeg(t *testing.T) {
	params := GetFourQCurve()
	var neg, res PointAffine
	neg.Neg(&params.Base)
	res.Add(&params.Base, &neg)
	if !res.IsZero() {
		t.Fatal("P + (-P) != O")
	}
}

func TestScalarMul(t *testing.T) {
	params := GetFourQCurve()

	// [1]*Base = Base
	var res PointAffine
	res.ScalarMultiplication(&params.Base, big.NewInt(1))
	if !res.Equal(&params.Base) {
		t.Fatal("[1]*Base != Base")
	}

	// [2]*Base = Double(Base)
	var dbl PointAffine
	dbl.Double(&params.Base)
	res.ScalarMultiplication(&params.Base, big.NewInt(2))
	if !res.Equal(&dbl) {
		t.Fatal("[2]*Base != Double(Base)")
	}

	// Random k: [k]*Base should be on curve
	k, _ := rand.Int(rand.Reader, &params.Order)
	res.ScalarMultiplication(&params.Base, k)
	if !res.IsOnCurve() {
		t.Fatal("[k]*Base not on curve")
	}
}

func TestSubgroupNaive(t *testing.T) {
	params := GetFourQCurve()

	if !params.Base.IsInSubGroupNaive() {
		t.Fatal("base point should be in subgroup")
	}

	var id PointAffine
	id.X.SetZero()
	id.Y.SetOne()
	if !id.IsInSubGroupNaive() {
		t.Fatal("identity should be in subgroup")
	}

	// Random k*Base
	k, _ := rand.Int(rand.Reader, &params.Order)
	var p PointAffine
	p.ScalarMultiplication(&params.Base, k)
	if !p.IsInSubGroupNaive() {
		t.Fatal("k*Base should be in subgroup")
	}

	// Non-subgroup point: (0, -1) has order 2
	var n PointAffine
	n.X.SetZero()
	n.Y.SetOne()
	n.Y.Neg(&n.Y)
	if n.IsOnCurve() && n.IsInSubGroupNaive() {
		t.Fatal("(0,-1) should NOT be in subgroup")
	}
}

func TestSubgroupTate(t *testing.T) {
	params := GetFourQCurve()

	// Base point should be in subgroup
	if !params.Base.IsInSubGroupTate() {
		t.Fatal("base point should be in subgroup (Tate)")
	}

	// Identity
	var id PointAffine
	id.X.SetZero()
	id.Y.SetOne()
	if !id.IsInSubGroupTate() {
		t.Fatal("identity should be in subgroup (Tate)")
	}

	// Random subgroup points
	nTests := 20
	for i := 0; i < nTests; i++ {
		k, _ := rand.Int(rand.Reader, &params.Order)
		var p PointAffine
		p.ScalarMultiplication(&params.Base, k)

		naive := p.IsInSubGroupNaive()
		tate := p.IsInSubGroupTate()

		if !naive || !tate {
			t.Fatalf("subgroup point: naive=%v tate=%v", naive, tate)
		}
	}

	// Non-subgroup: (0, -1)
	var n PointAffine
	n.X.SetZero()
	n.Y.SetOne()
	n.Y.Neg(&n.Y)
	if n.IsOnCurve() && n.IsInSubGroupTate() {
		t.Fatal("(0,-1) should NOT be in subgroup (Tate)")
	}

	t.Logf("tested %d subgroup points + non-subgroup: Naive and Tate agree", nTests)
}

func BenchmarkIsInSubGroupTate(b *testing.B) {
	params := GetFourQCurve()
	k, _ := rand.Int(rand.Reader, &params.Order)
	var p PointAffine
	p.ScalarMultiplication(&params.Base, k)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		p.IsInSubGroupTate()
	}
}

func TestClearCofactor(t *testing.T) {
	params := GetFourQCurve()

	// [392]*Base should still be on curve and in subgroup
	var cleared PointAffine
	cleared.Set(&params.Base)
	cleared.ClearCofactor()
	if !cleared.IsOnCurve() {
		t.Fatal("[392]*Base not on curve")
	}
	if !cleared.IsInSubGroupNaive() {
		t.Fatal("[392]*Base should be in subgroup")
	}

	// [392]*identity = identity
	var id PointAffine
	id.X.SetZero()
	id.Y.SetOne()
	id.ClearCofactor()
	if !id.IsZero() {
		t.Fatal("[392]*O != O")
	}

	// Random point on full curve: after cofactor clearing, result is in subgroup
	// Use k*Base (already in subgroup), so [392]*k*Base = k*[392]*Base is still in subgroup
	k, _ := rand.Int(rand.Reader, &params.Order)
	var p PointAffine
	p.ScalarMultiplication(&params.Base, k)
	p.ClearCofactor()
	if !p.IsInSubGroupNaive() {
		t.Fatal("[392]*k*Base should be in subgroup")
	}
}

func BenchmarkClearCofactor(b *testing.B) {
	params := GetFourQCurve()
	k, _ := rand.Int(rand.Reader, &params.Order)
	var p PointAffine
	p.ScalarMultiplication(&params.Base, k)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		var q PointAffine
		q.Set(&p)
		q.ClearCofactor()
	}
}

func BenchmarkIsInSubGroupNaive(b *testing.B) {
	params := GetFourQCurve()
	k, _ := rand.Int(rand.Reader, &params.Order)
	var p PointAffine
	p.ScalarMultiplication(&params.Base, k)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		p.IsInSubGroupNaive()
	}
}

func BenchmarkScalarMul(b *testing.B) {
	params := GetFourQCurve()
	k, _ := rand.Int(rand.Reader, &params.Order)
	var p PointAffine

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		p.ScalarMultiplication(&params.Base, k)
	}
}
