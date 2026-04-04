package curve25519

import (
	"crypto/rand"
	"math/big"
	"testing"

	fp "github.com/yelhousni/divide-and-pair/curve25519/fp"
)

func TestCurveParams(t *testing.T) {
	params := GetEdwardsCurve()

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

	// 8 * ell * Base = O (just checking cofactor)
	n := new(big.Int).Mul(&params.Order, big.NewInt(8))
	res.ScalarMultiplication(&params.Base, n)
	if !res.IsZero() {
		t.Fatal("8*ell * Base != O")
	}
}

func TestSubgroupNaive(t *testing.T) {
	params := GetEdwardsCurve()

	// Base point is in subgroup
	if !params.Base.IsInSubGroupNaive() {
		t.Fatal("base point should be in subgroup")
	}

	// Identity is in subgroup
	var id PointAffine
	id.X.SetZero()
	id.Y.SetOne()
	if !id.IsInSubGroupNaive() {
		t.Fatal("identity should be in subgroup")
	}

	// Random subgroup point: k * Base
	k, _ := rand.Int(rand.Reader, &params.Order)
	var p PointAffine
	p.ScalarMultiplication(&params.Base, k)
	if !p.IsInSubGroupNaive() {
		t.Fatal("k*Base should be in subgroup")
	}
}

func TestSubgroupPornin(t *testing.T) {
	params := GetEdwardsCurve()

	// Base point
	if !params.Base.IsInSubGroupPornin() {
		t.Fatal("base point should be in subgroup (Pornin)")
	}
}

func TestSubgroupQuartic(t *testing.T) {
	params := GetEdwardsCurve()

	// Base point
	if !params.Base.IsInSubGroupQuartic() {
		t.Fatal("base point should be in subgroup (quartic)")
	}
}

func TestSubgroupAgreement(t *testing.T) {
	params := GetEdwardsCurve()

	// Test subgroup points: k * Base (always in subgroup since Base has order ell)
	nTests := 50
	for i := 0; i < nTests; i++ {
		k, _ := rand.Int(rand.Reader, &params.Order)
		var p PointAffine
		p.ScalarMultiplication(&params.Base, k)

		naive := p.IsInSubGroupNaive()
		pornin := p.IsInSubGroupPornin()
		quartic := p.IsInSubGroupQuartic()

		if !naive || !pornin || !quartic {
			t.Fatalf("subgroup point k*Base: naive=%v pornin=%v quartic=%v", naive, pornin, quartic)
		}
	}
	t.Logf("tested %d subgroup points: all 3 methods agree (all true)", nTests)

	// Test non-subgroup points: add a low-order component
	// N = (0, -1) has order 2. k*Base + N is NOT in subgroup (unless k*Base = -N, negligible)
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
		quartic := q.IsInSubGroupQuartic()

		if naive || pornin || quartic {
			t.Fatalf("non-subgroup point k*Base+N: naive=%v pornin=%v quartic=%v", naive, pornin, quartic)
		}
		nNonSub++
	}
	t.Logf("tested %d non-subgroup points: all 3 methods agree (all false)", nNonSub)
}

func TestQuarticSymbolFromSubgroup(t *testing.T) {
	// quarticSymbol(1) = 1
	var one fp.Element
	one.SetOne()
	if one.QuarticSymbol() != 0 {
		t.Fatal("quarticSymbol(1) should be 0 (=1)")
	}

	// quarticSymbol is multiplicative: quarticSymbol(a*b) = quarticSymbol(a) + quarticSymbol(b) mod 4
	var a, b, ab fp.Element
	a.SetUint64(3)
	b.SetUint64(7)
	ab.Mul(&a, &b)
	ca := a.QuarticSymbol()
	cb := b.QuarticSymbol()
	cab := ab.QuarticSymbol()
	if cab != (ca+cb)%4 {
		t.Fatalf("quarticSymbol not multiplicative: quarticSymbol(%d)=%d, quarticSymbol(%d)=%d, quarticSymbol(%d*%d)=%d",
			3, ca, 7, cb, 3, 7, cab)
	}
}

func TestLowOrderPoints(t *testing.T) {
	subgroupInitOnce.Do(initSubgroupConstants)

	// Identity (0, 1) - in subgroup
	var id PointAffine
	id.X.SetZero()
	id.Y.SetOne()
	if !id.IsInSubGroupPornin() {
		t.Fatal("identity should be in subgroup (Pornin)")
	}
	if !id.IsInSubGroupQuartic() {
		t.Fatal("identity should be in subgroup (quartic)")
	}

	// Order-2 point: (0, -1) - NOT in subgroup
	var n PointAffine
	n.X.SetZero()
	n.Y.SetOne()
	n.Y.Neg(&n.Y)
	if n.IsInSubGroupPornin() {
		t.Fatal("(0,-1) should NOT be in subgroup (Pornin)")
	}
	if n.IsInSubGroupQuartic() {
		t.Fatal("(0,-1) should NOT be in subgroup (quartic)")
	}

	// 4-torsion: (i, 0)
	var t4 PointAffine
	t4.X.Set(&sqrtMinOne)
	t4.Y.SetZero()
	if t4.IsOnCurve() {
		// (i, 0): check a*i^2 + 0 = 1 + d*i^2*0 => a*(-1) = 1 => -a = 1 => a = -1. Yes!
		if t4.IsInSubGroupPornin() {
			t.Fatal("(i,0) should NOT be in subgroup (Pornin)")
		}
		if t4.IsInSubGroupQuartic() {
			t.Fatal("(i,0) should NOT be in subgroup (quartic)")
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

func BenchmarkIsInSubGroupQuarticExp(b *testing.B) {
	params := GetEdwardsCurve()
	k, _ := rand.Int(rand.Reader, &params.Order)
	var p PointAffine
	p.ScalarMultiplication(&params.Base, k)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		p.IsInSubGroupQuarticExp()
	}
}

func BenchmarkIsInSubGroupQuartic(b *testing.B) {
	params := GetEdwardsCurve()
	k, _ := rand.Int(rand.Reader, &params.Order)
	var p PointAffine
	p.ScalarMultiplication(&params.Base, k)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		p.IsInSubGroupQuartic()
	}
}
