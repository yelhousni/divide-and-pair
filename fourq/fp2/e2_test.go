package fp

import (
	"math/big"
	"testing"
)

func TestE2MulOne(t *testing.T) {
	var a, one, res E2
	a.A0.SetUint64(3)
	a.A1.SetUint64(7)
	one.SetOne()
	res.Mul(&a, &one)
	if !res.Equal(&a) {
		t.Fatal("a * 1 != a")
	}
}

func TestE2MulInverse(t *testing.T) {
	var a, inv, res E2
	a.A0.SetUint64(3)
	a.A1.SetUint64(7)
	inv.Inverse(&a)
	res.Mul(&a, &inv)
	if !res.IsOne() {
		t.Fatalf("a * a^-1 != 1, got %s", res.String())
	}
}

func TestE2Square(t *testing.T) {
	var a, sq, mul E2
	a.A0.SetUint64(5)
	a.A1.SetUint64(11)
	sq.Square(&a)
	mul.Mul(&a, &a)
	if !sq.Equal(&mul) {
		t.Fatal("a^2 != a*a")
	}
}

func TestE2Conjugate(t *testing.T) {
	// (a+bi)(a-bi) = a^2 + b^2 (real)
	var a, conj, prod E2
	a.A0.SetUint64(3)
	a.A1.SetUint64(7)
	conj.Conjugate(&a)
	prod.Mul(&a, &conj)
	if !prod.A1.IsZero() {
		t.Fatal("a * conj(a) should be real")
	}
	norm := a.Norm()
	if !prod.A0.Equal(&norm) {
		t.Fatal("a * conj(a) should equal norm(a)")
	}
}

func TestE2Exp(t *testing.T) {
	var a, res E2
	a.A0.SetUint64(2)
	a.A1.SetUint64(3)
	// a^0 = 1
	res.Exp(&a, big.NewInt(0))
	if !res.IsOne() {
		t.Fatal("a^0 != 1")
	}
	// a^1 = a
	res.Exp(&a, big.NewInt(1))
	if !res.Equal(&a) {
		t.Fatal("a^1 != a")
	}
	// a^2 = a*a
	var sq E2
	sq.Square(&a)
	res.Exp(&a, big.NewInt(2))
	if !res.Equal(&sq) {
		t.Fatal("a^2 != a*a")
	}
}

func TestE2Random(t *testing.T) {
	for i := 0; i < 20; i++ {
		var a, inv, prod E2
		a.SetRandom()
		if a.IsZero() {
			continue
		}
		inv.Inverse(&a)
		prod.Mul(&a, &inv)
		if !prod.IsOne() {
			t.Fatalf("random a * a^-1 != 1")
		}
	}
}
