package fp2

import (
	"testing"
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
