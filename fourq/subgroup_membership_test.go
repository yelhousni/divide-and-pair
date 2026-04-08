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

func TestIsInSubGroupTateExp1Exp2Exp3(t *testing.T) {
	params := curveParameters()

	nTests := 20
	for i := 0; i < nTests; i++ {
		k, _ := rand.Int(rand.Reader, &params.Order)
		var p PointAffine
		p.ScalarMultiplication(&params.Base, k)

		tate := p.isInSubGroupTate()
		exp1 := p.isInSubGroupTateExp1()
		exp2 := p.isInSubGroupTateExp2()
		exp3 := p.isInSubGroupTateExp3()

		if !tate || !exp1 || !exp2 || !exp3 {
			t.Fatalf("subgroup point: tate=%v exp1=%v exp2=%v exp3=%v", tate, exp1, exp2, exp3)
		}
	}

	// Non-subgroup: (0, -1)
	var n PointAffine
	n.X.SetZero()
	n.Y.SetOne()
	n.Y.Neg(&n.Y)

	for i := 0; i < nTests; i++ {
		k, _ := rand.Int(rand.Reader, &params.Order)
		if k.Sign() == 0 {
			continue
		}
		var p, q PointAffine
		p.ScalarMultiplication(&params.Base, k)
		q.Add(&p, &n)

		tate := q.isInSubGroupTate()
		exp1 := q.isInSubGroupTateExp1()
		exp2 := q.isInSubGroupTateExp2()
		exp3 := q.isInSubGroupTateExp3()

		if tate || exp1 || exp2 || exp3 {
			t.Fatalf("non-subgroup: tate=%v exp1=%v exp2=%v exp3=%v", tate, exp1, exp2, exp3)
		}
	}
	t.Logf("tested %d subgroup + %d non-subgroup: all methods agree", nTests, nTests)
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

func BenchmarkIsInSubGroupTateExp1(b *testing.B) {
	params := curveParameters()
	k, _ := rand.Int(rand.Reader, &params.Order)
	var p PointAffine
	p.ScalarMultiplication(&params.Base, k)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		p.isInSubGroupTateExp1()
	}
}

func BenchmarkIsInSubGroupTateExp2(b *testing.B) {
	params := curveParameters()
	k, _ := rand.Int(rand.Reader, &params.Order)
	var p PointAffine
	p.ScalarMultiplication(&params.Base, k)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		p.isInSubGroupTateExp2()
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

func BenchmarkIsInSubGroupTateExp3(b *testing.B) {
	params := curveParameters()
	k, _ := rand.Int(rand.Reader, &params.Order)
	var p PointAffine
	p.ScalarMultiplication(&params.Base, k)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		p.isInSubGroupTateExp3()
	}
}
