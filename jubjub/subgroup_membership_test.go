package jubjub

import (
	"crypto/rand"
	"math/big"
	"testing"
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

func TestSubgroupNaive(t *testing.T) {
	params := curveParameters()

	if !params.Base.isInSubGroupNaive() {
		t.Fatal("base point should be in subgroup")
	}

	var id PointAffine
	id.X.SetZero()
	id.Y.SetOne()
	if !id.isInSubGroupNaive() {
		t.Fatal("identity should be in subgroup")
	}

	k, _ := rand.Int(rand.Reader, &params.Order)
	var p PointAffine
	p.ScalarMultiplication(&params.Base, k)
	if !p.isInSubGroupNaive() {
		t.Fatal("k*Base should be in subgroup")
	}
}

func TestSubgroupPornin(t *testing.T) {
	params := curveParameters()

	if !params.Base.isInSubGroupPornin() {
		t.Fatal("base point should be in subgroup (Pornin)")
	}
}

func TestSubgroupOcticExp(t *testing.T) {
	params := curveParameters()

	if !params.Base.isInSubGroupOcticExp() {
		t.Fatal("base point should be in subgroup (OcticExp)")
	}
}

func TestSubgroupAgreement(t *testing.T) {
	params := curveParameters()

	nTests := 50
	for i := 0; i < nTests; i++ {
		k, _ := rand.Int(rand.Reader, &params.Order)
		var p PointAffine
		p.ScalarMultiplication(&params.Base, k)

		naive := p.isInSubGroupNaive()
		pornin := p.isInSubGroupPornin()
		octic := p.isInSubGroupOcticExp()

		if !naive || !pornin || !octic {
			t.Fatalf("subgroup point k*Base: naive=%v pornin=%v octic=%v", naive, pornin, octic)
		}
	}
	t.Logf("tested %d subgroup points: all 3 methods agree (all true)", nTests)

	var nPt PointAffine
	nPt.X.SetZero()
	nPt.Y.SetOne()
	nPt.Y.Neg(&nPt.Y)

	nNonSub := 0
	for i := 0; i < nTests; i++ {
		k, _ := rand.Int(rand.Reader, &params.Order)
		if k.Sign() == 0 {
			continue
		}
		var p, q PointAffine
		p.ScalarMultiplication(&params.Base, k)
		q.Add(&p, &nPt)

		naive := q.isInSubGroupNaive()
		pornin := q.isInSubGroupPornin()
		octic := q.isInSubGroupOcticExp()

		if naive || pornin || octic {
			t.Fatalf("non-subgroup point k*Base+N: naive=%v pornin=%v octic=%v", naive, pornin, octic)
		}
		nNonSub++
	}
	t.Logf("tested %d non-subgroup points: all 3 methods agree (all false)", nNonSub)
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
