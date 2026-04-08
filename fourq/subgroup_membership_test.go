package fourq

import (
	"crypto/rand"
	"testing"
)

func TestSubgroupNaive(t *testing.T) {
	params := curveParameters()

	if !params.Base.isInSubGroupNaive() {
		t.Fatal("base point should be in subgroup (naive)")
	}

	var id PointAffine
	id.SetInfinity()
	if !id.isInSubGroupNaive() {
		t.Fatal("identity should be in subgroup (naive)")
	}
}

func TestIsInSubGroupTate(t *testing.T) {
	params := curveParameters()

	// Identity
	var id PointAffine
	id.SetInfinity()
	if !id.isInSubGroupTate() {
		t.Fatal("identity should be in subgroup (Tate)")
	}

	// Random subgroup points
	nTests := 20
	for i := 0; i < nTests; i++ {
		k, _ := rand.Int(rand.Reader, &params.Order)
		var p PointAffine
		p.ScalarMultiplication(&params.Base, k)

		naive := p.isInSubGroupNaive()
		tate := p.isInSubGroupTate()

		if !naive || !tate {
			t.Fatalf("subgroup point: naive=%v tate=%v", naive, tate)
		}
	}

	// Non-subgroup: (0, -1)
	var n PointAffine
	n.SetInfinity()
	n.Y.Neg(&n.Y)
	if n.IsOnCurve() && n.isInSubGroupTate() {
		t.Fatal("(0,-1) should NOT be in subgroup (Tate)")
	}

	t.Logf("tested %d subgroup points + non-subgroup: Naive and Tate agree", nTests)
}

// Benchmarks

func BenchmarkIsInSubGroupNaive(b *testing.B) {
	params := curveParameters()
	k, _ := rand.Int(rand.Reader, &params.Order)
	var p PointAffine
	p.ScalarMultiplication(&params.Base, k)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		p.isInSubGroupNaive()
	}
}

func BenchmarkIsInSubGroupTate(b *testing.B) {
	params := curveParameters()
	k, _ := rand.Int(rand.Reader, &params.Order)
	var p PointAffine
	p.ScalarMultiplication(&params.Base, k)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		p.isInSubGroupTate()
	}
}

func BenchmarkIsInSubGroupEndo(b *testing.B) {
	params := curveParameters()
	k, _ := rand.Int(rand.Reader, &params.Order)
	var p PointAffine
	p.ScalarMultiplication(&params.Base, k)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		p.isInSubGroupEndo()
	}
}
