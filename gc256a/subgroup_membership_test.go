package gc256a

import (
	"crypto/rand"
	"math/big"
	"testing"
)

func TestCurveParams(t *testing.T) {
	params := GetEdwardsCurve()

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

func TestSubgroupNaive(t *testing.T) {
	params := GetEdwardsCurve()

	if !params.Base.IsInSubGroupNaive() {
		t.Fatal("base point should be in subgroup")
	}

	var id PointAffine
	id.X.SetZero()
	id.Y.SetOne()
	if !id.IsInSubGroupNaive() {
		t.Fatal("identity should be in subgroup")
	}

	k, _ := rand.Int(rand.Reader, &params.Order)
	var p PointAffine
	p.ScalarMultiplication(&params.Base, k)
	if !p.IsInSubGroupNaive() {
		t.Fatal("k*Base should be in subgroup")
	}
}

func TestSubgroupPornin(t *testing.T) {
	params := GetEdwardsCurve()

	if !params.Base.IsInSubGroupPornin() {
		t.Fatal("base point should be in subgroup (Pornin)")
	}
}

func TestSubgroupAgreement(t *testing.T) {
	params := GetEdwardsCurve()

	nTests := 50
	for i := 0; i < nTests; i++ {
		k, _ := rand.Int(rand.Reader, &params.Order)
		var p PointAffine
		p.ScalarMultiplication(&params.Base, k)

		naive := p.IsInSubGroupNaive()
		pornin := p.IsInSubGroupPornin()

		if !naive || !pornin {
			t.Fatalf("subgroup point k*Base: naive=%v pornin=%v", naive, pornin)
		}
	}
	t.Logf("tested %d subgroup points: both methods agree (all true)", nTests)

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

		naive := q.IsInSubGroupNaive()
		pornin := q.IsInSubGroupPornin()

		if naive || pornin {
			t.Fatalf("non-subgroup point k*Base+N: naive=%v pornin=%v", naive, pornin)
		}
		nNonSub++
	}
	t.Logf("tested %d non-subgroup points: both methods agree (all false)", nNonSub)
}

func TestLowOrderPoints(t *testing.T) {
	subgroupInitOnce.Do(initSubgroupConstants)

	var id PointAffine
	id.X.SetZero()
	id.Y.SetOne()
	if !id.IsInSubGroupPornin() {
		t.Fatal("identity should be in subgroup (Pornin)")
	}

	var n PointAffine
	n.X.SetZero()
	n.Y.SetOne()
	n.Y.Neg(&n.Y)
	if n.IsInSubGroupPornin() {
		t.Fatal("(0,-1) should NOT be in subgroup (Pornin)")
	}

	var t4 PointAffine
	t4.X.SetOne()
	t4.Y.SetZero()
	if t4.IsOnCurve() {
		if t4.IsInSubGroupPornin() {
			t.Fatal("(1,0) should NOT be in subgroup (Pornin)")
		}
	}

	var t4neg PointAffine
	t4neg.X.SetOne()
	t4neg.X.Neg(&t4neg.X)
	t4neg.Y.SetZero()
	if t4neg.IsOnCurve() {
		if t4neg.IsInSubGroupPornin() {
			t.Fatal("(-1,0) should NOT be in subgroup (Pornin)")
		}
	}
}

// Benchmarks

func BenchmarkIsInSubGroupNaive(b *testing.B) {
	params := GetEdwardsCurve()
	k, _ := rand.Int(rand.Reader, &params.Order)
	var p PointAffine
	p.ScalarMultiplication(&params.Base, k)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		p.IsInSubGroupNaive()
	}
}

func BenchmarkIsInSubGroupPornin(b *testing.B) {
	params := GetEdwardsCurve()
	k, _ := rand.Int(rand.Reader, &params.Order)
	var p PointAffine
	p.ScalarMultiplication(&params.Base, k)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		p.IsInSubGroupPornin()
	}
}
