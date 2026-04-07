package fourq

import (
	"crypto/rand"
	"math/big"
	"testing"
)

func TestCurveParams(t *testing.T) {
	params := curveParameters()
	g := Generators()
	a, d := CurveCoefficients()
	order := Order()

	if !g.Equal(&params.Base) {
		t.Fatal("Generators() returned unexpected base point")
	}
	if !a.Equal(&params.A) || !d.Equal(&params.D) {
		t.Fatal("CurveCoefficients() returned unexpected coefficients")
	}
	if order.Cmp(&params.Order) != 0 {
		t.Fatal("Order() returned unexpected subgroup order")
	}

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
	id.SetInfinity()
	if !id.IsZero() {
		t.Fatal("(0,1) should be identity")
	}
	if !id.IsOnCurve() {
		t.Fatal("identity should be on curve")
	}
}

func TestAddIdentity(t *testing.T) {
	params := curveParameters()
	var id, res PointAffine
	id.SetInfinity()
	res.Add(&params.Base, &id)
	if !res.Equal(&params.Base) {
		t.Fatal("P + O != P")
	}
}

func TestDoubleConsistency(t *testing.T) {
	params := curveParameters()
	var d1, d2 PointAffine
	d1.Double(&params.Base)
	d2.Add(&params.Base, &params.Base)
	if !d1.Equal(&d2) {
		t.Fatal("Double(P) != P + P")
	}
}

func TestNeg(t *testing.T) {
	params := curveParameters()
	var neg, res PointAffine
	neg.Neg(&params.Base)
	res.Add(&params.Base, &neg)
	if !res.IsZero() {
		t.Fatal("P + (-P) != O")
	}
}

func TestScalarMul(t *testing.T) {
	params := curveParameters()

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
	params := curveParameters()

	if !params.Base.isInSubGroupNaive() {
		t.Fatal("base point should be in subgroup")
	}

	var id PointAffine
	id.SetInfinity()
	if !id.isInSubGroupNaive() {
		t.Fatal("identity should be in subgroup")
	}

	// Random k*Base
	k, _ := rand.Int(rand.Reader, &params.Order)
	var p PointAffine
	p.ScalarMultiplication(&params.Base, k)
	if !p.isInSubGroupNaive() {
		t.Fatal("k*Base should be in subgroup")
	}

	// Non-subgroup point: (0, -1) has order 2
	var n PointAffine
	n.SetInfinity()
	n.Y.Neg(&n.Y)
	if n.IsOnCurve() && n.isInSubGroupNaive() {
		t.Fatal("(0,-1) should NOT be in subgroup")
	}
}

func TestIsInSubGroupTate(t *testing.T) {
	params := curveParameters()

	// Base point should be in subgroup
	if !params.Base.isInSubGroupTate() {
		t.Fatal("base point should be in subgroup (Tate)")
	}
	if !params.Base.IsInSubGroup() {
		t.Fatal("base point should be in subgroup")
	}

	// Identity
	var id PointAffine
	id.SetInfinity()
	if !id.isInSubGroupTate() {
		t.Fatal("identity should be in subgroup (Tate)")
	}
	if !id.IsInSubGroup() {
		t.Fatal("identity should be in subgroup")
	}

	// Random subgroup points
	nTests := 20
	for i := 0; i < nTests; i++ {
		k, _ := rand.Int(rand.Reader, &params.Order)
		var p PointAffine
		p.ScalarMultiplication(&params.Base, k)

		naive := p.isInSubGroupNaive()
		tate := p.isInSubGroupTate()
		defaultCheck := p.IsInSubGroup()

		if !naive || !tate || !defaultCheck {
			t.Fatalf("subgroup point: naive=%v tate=%v default=%v", naive, tate, defaultCheck)
		}
	}

	// Non-subgroup: (0, -1)
	var n PointAffine
	n.SetInfinity()
	n.Y.Neg(&n.Y)
	if n.IsOnCurve() && n.isInSubGroupTate() {
		t.Fatal("(0,-1) should NOT be in subgroup (Tate)")
	}
	if n.IsOnCurve() && n.IsInSubGroup() {
		t.Fatal("(0,-1) should NOT be in subgroup")
	}

	t.Logf("tested %d subgroup points + non-subgroup: Naive and Tate agree", nTests)
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

func TestClearCofactor(t *testing.T) {
	params := curveParameters()

	// [392]*Base should still be on curve and in subgroup
	var cleared PointAffine
	cleared.Set(&params.Base)
	cleared.ClearCofactor()
	if !cleared.IsOnCurve() {
		t.Fatal("[392]*Base not on curve")
	}
	if !cleared.isInSubGroupNaive() {
		t.Fatal("[392]*Base should be in subgroup")
	}

	// [392]*identity = identity
	var id PointAffine
	id.SetInfinity()
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
	if !p.isInSubGroupNaive() {
		t.Fatal("[392]*k*Base should be in subgroup")
	}
}

func BenchmarkClearCofactor(b *testing.B) {
	params := curveParameters()
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

func TestIsInSubGroupEndo(t *testing.T) {
	params := curveParameters()

	// Base point should be in subgroup
	if !params.Base.isInSubGroupEndo() {
		t.Fatal("base point should be in subgroup (Endo)")
	}

	// Identity
	var id PointAffine
	id.SetInfinity()
	if !id.isInSubGroupEndo() {
		t.Fatal("identity should be in subgroup (Endo)")
	}

	// Random subgroup points: endo test should agree with naive
	nTests := 10
	for i := 0; i < nTests; i++ {
		k, _ := rand.Int(rand.Reader, &params.Order)
		var p PointAffine
		p.ScalarMultiplication(&params.Base, k)

		if !p.isInSubGroupEndo() {
			t.Fatalf("subgroup point not detected by endo test (i=%d)", i)
		}
	}

	// Non-subgroup point: (0, -1) has order 2
	var n PointAffine
	n.SetInfinity()
	n.Y.Neg(&n.Y)
	if n.IsOnCurve() && n.isInSubGroupEndo() {
		t.Fatal("(0,-1) should NOT be in subgroup (Endo)")
	}

	t.Logf("tested %d subgroup + 1 non-subgroup: all correct", nTests)
}

func TestEndomorphismEigenvalues(t *testing.T) {
	params := curveParameters()

	// Verify ψ(G) == [λ_ψ]G and φ(G) == [λ_φ]G on the base point
	var psiG, phiG, expected PointAffine

	psiG.Psi(&params.Base)
	expected.ScalarMultiplication(&params.Base, &lambdaPsi)
	if !psiG.Equal(&expected) {
		t.Fatal("ψ(G) ≠ [λ_ψ]G")
	}

	phiG.Phi(&params.Base)
	expected.ScalarMultiplication(&params.Base, &lambdaPhi)
	if !phiG.Equal(&expected) {
		t.Fatal("φ(G) ≠ [λ_φ]G")
	}

	t.Log("eigenvalue verification passed")
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

func BenchmarkScalarMul(b *testing.B) {
	params := curveParameters()
	k, _ := rand.Int(rand.Reader, &params.Order)
	var p PointAffine

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		p.ScalarMultiplication(&params.Base, k)
	}
}
