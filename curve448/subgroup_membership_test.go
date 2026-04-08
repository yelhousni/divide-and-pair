package curve448

import (
	"crypto/rand"
	"math/big"
	"testing"
)

func TestCurveParams(t *testing.T) {
	params := curveParameters()

	// Base point should be on the curve
	if !params.Base.IsOnCurve() {
		t.Fatal("base point not on curve")
	}

	// ell * Base should be identity
	var res PointAffine
	res.ScalarMultiplication(&params.Base, &params.Order)
	if !res.IsZero() {
		t.Fatal("ell * Base != O")
	}

	// 4 * ell * Base = O (just checking cofactor)
	n := new(big.Int).Mul(&params.Order, big.NewInt(4))
	res.ScalarMultiplication(&params.Base, n)
	if !res.IsZero() {
		t.Fatal("4*ell * Base != O")
	}
}

func TestSubgroupNaive(t *testing.T) {
	params := curveParameters()

	// Base point is in subgroup
	if !params.Base.isInSubGroupNaive() {
		t.Fatal("base point should be in subgroup")
	}

	// Identity is in subgroup
	var id PointAffine
	id.X.SetZero()
	id.Y.SetOne()
	if !id.isInSubGroupNaive() {
		t.Fatal("identity should be in subgroup")
	}

	// Random subgroup point: k * Base
	k, _ := rand.Int(rand.Reader, &params.Order)
	var p PointAffine
	p.ScalarMultiplication(&params.Base, k)
	if !p.isInSubGroupNaive() {
		t.Fatal("k*Base should be in subgroup")
	}
}

func TestSubgroupPornin(t *testing.T) {
	params := curveParameters()

	// Base point
	if !params.Base.isInSubGroupPornin() {
		t.Fatal("base point should be in subgroup (Pornin)")
	}
}

func TestSubgroupAgreement(t *testing.T) {
	params := curveParameters()

	// Test subgroup points: k * Base
	nTests := 50
	for i := 0; i < nTests; i++ {
		k, _ := rand.Int(rand.Reader, &params.Order)
		var p PointAffine
		p.ScalarMultiplication(&params.Base, k)

		naive := p.isInSubGroupNaive()
		pornin := p.isInSubGroupPornin()

		if !naive || !pornin {
			t.Fatalf("subgroup point k*Base: naive=%v pornin=%v", naive, pornin)
		}
	}
	t.Logf("tested %d subgroup points: both methods agree (all true)", nTests)

	// Test non-subgroup points: add a low-order component
	// N = (0, -1) has order 2. k*Base + N is NOT in subgroup
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

func TestSubgroupQuarticAgreement(t *testing.T) {
	params := curveParameters()

	nTests := 50
	for i := 0; i < nTests; i++ {
		k, _ := rand.Int(rand.Reader, &params.Order)
		var p PointAffine
		p.ScalarMultiplication(&params.Base, k)

		pornin := p.isInSubGroupPornin()
		quartic := p.isInSubGroupQuartic()

		if !pornin || !quartic {
			t.Fatalf("subgroup point: pornin=%v quartic=%v", pornin, quartic)
		}
	}
	t.Logf("tested %d subgroup points: both methods agree (all true)", nTests)

	// Non-subgroup: order-2 offset
	var order2Pt PointAffine
	order2Pt.X.SetZero()
	order2Pt.Y.SetOne()
	order2Pt.Y.Neg(&order2Pt.Y)

	nNonSub := 0
	for i := 0; i < nTests; i++ {
		k, _ := rand.Int(rand.Reader, &params.Order)
		if k.Sign() == 0 {
			continue
		}
		var p, q PointAffine
		p.ScalarMultiplication(&params.Base, k)
		q.Add(&p, &order2Pt)

		pornin := q.isInSubGroupPornin()
		quartic := q.isInSubGroupQuartic()

		if pornin || quartic {
			t.Fatalf("non-subgroup (order-2): pornin=%v quartic=%v", pornin, quartic)
		}
		nNonSub++
	}
	t.Logf("tested %d non-subgroup points (order-2): all false", nNonSub)

	// Non-subgroup: order-4 offset
	var order4Pt PointAffine
	order4Pt.X.SetOne()
	order4Pt.Y.SetZero()

	nNonSub = 0
	for i := 0; i < nTests; i++ {
		k, _ := rand.Int(rand.Reader, &params.Order)
		if k.Sign() == 0 {
			continue
		}
		var p, q PointAffine
		p.ScalarMultiplication(&params.Base, k)
		q.Add(&p, &order4Pt)

		pornin := q.isInSubGroupPornin()
		quartic := q.isInSubGroupQuartic()

		if pornin || quartic {
			t.Fatalf("non-subgroup (order-4): pornin=%v quartic=%v", pornin, quartic)
		}
		nNonSub++
	}
	t.Logf("tested %d non-subgroup points (order-4): all false", nNonSub)
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

func BenchmarkIsInSubGroupPornin(b *testing.B) {
	params := curveParameters()
	k, _ := rand.Int(rand.Reader, &params.Order)
	var p PointAffine
	p.ScalarMultiplication(&params.Base, k)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		p.isInSubGroupPornin()
	}
}

func BenchmarkIsInSubGroupQuartic(b *testing.B) {
	params := curveParameters()
	k, _ := rand.Int(rand.Reader, &params.Order)
	var p PointAffine
	p.ScalarMultiplication(&params.Base, k)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		p.isInSubGroupQuartic()
	}
}
