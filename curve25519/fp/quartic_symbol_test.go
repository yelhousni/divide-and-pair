package fp

import (
	"crypto/rand"
	"testing"
)

func TestQuarticSymbolKnown(t *testing.T) {
	// quarticSymbol(1) = 1
	var one Element
	one.SetOne()
	if c := one.QuarticSymbol(); c != 0 {
		t.Fatalf("quarticSymbol(1) = %d, want 0", c)
	}

	// quarticSymbol(0) = 0 (convention)
	var zero Element
	if c := zero.QuarticSymbol(); c != 0 {
		t.Fatalf("quarticSymbol(0) = %d, want 0", c)
	}

	// quarticSymbol(-1) = (-1)^((q-1)/4); for q=2^255-19, (q-1)/4 is odd, so quarticSymbol(-1) = -1 = i^2
	var minusOne Element
	minusOne.SetOne()
	minusOne.Neg(&minusOne)
	if c := minusOne.QuarticSymbol(); c != 2 {
		t.Fatalf("quarticSymbol(-1) = %d, want 2 (= -1)", c)
	}
}

func TestQuarticSymbolMultiplicative(t *testing.T) {
	// quarticSymbol(a*b) = quarticSymbol(a) + quarticSymbol(b) mod 4
	for i := 0; i < 100; i++ {
		var a, b, ab Element
		a.SetRandom()
		b.SetRandom()
		ab.Mul(&a, &b)

		ca := a.QuarticSymbol()
		cb := b.QuarticSymbol()
		cab := ab.QuarticSymbol()

		if cab != (ca+cb)%4 {
			t.Fatalf("quarticSymbol not multiplicative: quarticSymbol(a)=%d, quarticSymbol(b)=%d, quarticSymbol(a*b)=%d, expected %d",
				ca, cb, cab, (ca+cb)%4)
		}
	}
}

func TestQuarticSymbolVsLegendre(t *testing.T) {
	// quarticSymbol(x)^2 = Legendre(x): (quarticSymbol mod 2) should give 0 for QR, 1 for QNR
	for i := 0; i < 100; i++ {
		var x Element
		x.SetRandom()

		c := x.QuarticSymbol()
		l := x.Legendre()

		// quarticSymbol^2: 0->1, 1->-1, 2->1, 3->-1
		var expectedLeg int
		switch c {
		case 0, 2:
			expectedLeg = 1
		case 1, 3:
			expectedLeg = -1
		}

		if l != expectedLeg {
			t.Fatalf("quarticSymbol=%d implies Legendre=%d, but got Legendre=%d", c, expectedLeg, l)
		}
	}
}

func TestQuarticSymbolVsExp(t *testing.T) {
	// Compare QuarticSymbol (Weilert Euclidean GCD) against quarticSymbolExp (addition chain)
	for v := uint64(1); v <= 100; v++ {
		var x Element
		x.SetUint64(v)
		gcd := x.QuarticSymbol()
		exp := x.QuarticSymbolExp()
		if gcd != exp {
			t.Fatalf("quarticSymbol(%d): GCD=%d, Exp=%d", v, gcd, exp)
		}
	}

	for i := 0; i < 1000; i++ {
		var x Element
		x.SetRandom()
		gcd := x.QuarticSymbol()
		exp := x.QuarticSymbolExp()
		if gcd != exp {
			t.Fatalf("mismatch at i=%d: QuarticSymbol=%d, quarticSymbolExp=%d", i, gcd, exp)
		}
	}
	t.Log("1100/1100 QuarticSymbol vs quarticSymbolExp match")
}

func TestQuarticSymbolConsistency(t *testing.T) {
	for i := 0; i < 200; i++ {
		var x Element
		x.SetRandom()
		c := x.QuarticSymbol()
		l := x.Legendre()
		var expectedLeg int
		switch c {
		case 0, 2:
			expectedLeg = 1
		case 1, 3:
			expectedLeg = -1
		}
		if l != expectedLeg {
			t.Fatalf("quarticSymbol=%d implies Legendre=%d, but got Legendre=%d", c, expectedLeg, l)
		}
	}
}

func BenchmarkQuarticSymbol(b *testing.B) {
	var x Element
	x.SetRandom()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		x.QuarticSymbol()
	}
}

func BenchmarkQuarticSymbolExp(b *testing.B) {
	var x Element
	x.SetRandom()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		x.QuarticSymbolExp()
	}
}

func BenchmarkLegendre(b *testing.B) {
	var x Element
	x.SetRandom()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		x.Legendre()
	}
}

func BenchmarkSqrt(b *testing.B) {
	var x Element
	x.SetRandom()
	// Make x a QR
	x.Square(&x)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		var r Element
		r.Sqrt(&x)
	}
}

// Generate random elements for testing
func init() {
	// Ensure we have a working random source
	_ = rand.Reader
}
