package fp2

import (
	"math/big"
	"testing"

	fp "github.com/yelhousni/divide-and-pair/curve448/fp"
)

func TestE2Mul(t *testing.T) {
	var a, b, c, d E2
	a.SetRandom()
	b.SetRandom()

	// Commutativity
	c.Mul(&a, &b)
	d.Mul(&b, &a)
	if !c.Equal(&d) {
		t.Fatal("Mul not commutative")
	}

	// Associativity
	var e E2
	e.SetRandom()
	c.Mul(&a, &b)
	c.Mul(&c, &e)
	d.Mul(&b, &e)
	d.Mul(&a, &d)
	if !c.Equal(&d) {
		t.Fatal("Mul not associative")
	}
}

func TestE2Square(t *testing.T) {
	var a, sq1, sq2 E2
	a.SetRandom()
	sq1.Square(&a)
	sq2.Mul(&a, &a)
	if !sq1.Equal(&sq2) {
		t.Fatal("Square != Mul(a,a)")
	}
}

func TestE2NormConjugate(t *testing.T) {
	var a, conj E2
	a.SetRandom()
	conj.Conjugate(&a)

	norm := a.Norm()
	var prod E2
	prod.Mul(&a, &conj)
	if !prod.A0.Equal(&norm) || !prod.A1.IsZero() {
		t.Fatal("Norm != a * conj(a)")
	}
}

func TestE2Inverse(t *testing.T) {
	var a, inv, prod E2
	a.SetRandom()
	inv.Inverse(&a)
	prod.Mul(&a, &inv)
	if !prod.IsOne() {
		t.Fatal("a * a^(-1) != 1")
	}
}

func TestE2ExpBySqrtPp1o4(t *testing.T) {
	// Verify ExpBySqrtPp1o4 against generic Exp with the same exponent.
	var a, fast, slow E2
	a.SetRandom()

	fast.ExpBySqrtPp1o4(&a)

	// (p+1)/4 as big.Int
	p := fp.Modulus()
	exp := new(big.Int).Add(p, big.NewInt(1))
	exp.Rsh(exp, 2) // (p+1)/4

	slow.Exp(&a, exp)

	if !fast.Equal(&slow) {
		t.Fatal("ExpBySqrtPp1o4 != generic Exp")
	}
}

func TestQuarticSymbol(t *testing.T) {
	var one E2
	one.SetOne()

	// chi_4(1) = 0 (1 is a quartic residue)
	if one.QuarticSymbol() != 0 {
		t.Fatalf("chi_4(1) = %d, want 0", one.QuarticSymbol())
	}

	// chi_4(z^4) = 0 for random z
	var z, z4 E2
	z.SetRandom()
	z4.Square(&z)
	z4.Square(&z4)
	if z4.QuarticSymbol() != 0 {
		t.Fatalf("chi_4(z^4) = %d, want 0", z4.QuarticSymbol())
	}

	// Multiplicativity: chi_4(a*b) = chi_4(a) + chi_4(b) mod 4
	var a, b, ab E2
	for i := 0; i < 20; i++ {
		a.SetRandom()
		b.SetRandom()
		ab.Mul(&a, &b)

		ca := a.QuarticSymbol()
		cb := b.QuarticSymbol()
		cab := ab.QuarticSymbol()

		if cab != (ca+cb)%4 {
			t.Fatalf("multiplicativity failed: chi4(a)=%d, chi4(b)=%d, chi4(ab)=%d, want %d",
				ca, cb, cab, (ca+cb)%4)
		}
	}
}

func BenchmarkE2Square(b *testing.B) {
	var z E2
	z.SetRandom()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		z.Square(&z)
	}
}

func BenchmarkE2Mul(b *testing.B) {
	var x, y E2
	x.SetRandom()
	y.SetRandom()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		x.Mul(&x, &y)
	}
}

func BenchmarkExpBySqrtPp1o4(b *testing.B) {
	var z E2
	z.SetRandom()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		z.ExpBySqrtPp1o4(&z)
	}
}
