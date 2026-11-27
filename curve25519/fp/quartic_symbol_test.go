package fp

import (
	"math/big"
	"math/rand"
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
	// Compare QuarticSymbol (Gaussian binary GCD) against quarticSymbolExp (addition chain)
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

// Exercise dyadic valuations around word boundaries and both sides of the
// rounded quotients used to lift the raw Montgomery representative.
func TestQuarticSymbolReductionBoundaries(t *testing.T) {
	modulus := Modulus()
	check := func(raw *big.Int) {
		var x Element
		x.SetBigInt(raw)
		// Feed these exact bits into the reduction, without the SetBigInt
		// Montgomery conversion. The exponentiation uses the same element.
		x = Element(x.Bits())
		if got, want := x.QuarticSymbol(), x.QuarticSymbolExp(); got != want {
			t.Fatalf("raw residue %x: got %d, want %d", raw, got, want)
		}
	}
	one := big.NewInt(1)
	imagUnit := sqrtMinusOneFp.BigInt(new(big.Int))
	for _, k := range []uint{0, 1, 59, 60, 61, 62, 63, 64, 65, 125, 126, 127, 128, 129, 190, 191, 192, 193, 253, 254} {
		v := new(big.Int).Lsh(one, k)
		// For k=128 one of these lifts is pi + 2^128 i. Its first
		// subtraction from pi shifts a full two words at once.
		imaginary := new(big.Int).Mul(v, imagUnit)
		check(imaginary)
		check(new(big.Int).Neg(imaginary))
		for _, delta := range []int64{-1, 0, 1} {
			raw := new(big.Int).Add(v, big.NewInt(delta))
			check(raw)
			check(new(big.Int).Sub(modulus, raw))
		}
	}
	for _, coord := range [][2]uint64{{piRe0, piRe1}, {piIm0, piIm1}} {
		c := new(big.Int).Lsh(new(big.Int).SetUint64(coord[1]), 64)
		c.Add(c, new(big.Int).SetUint64(coord[0]))
		// The first and last half-integer quotient boundaries.
		for _, twiceQ := range []*big.Int{one, new(big.Int).Sub(new(big.Int).Lsh(c, 1), one)} {
			boundary := new(big.Int).Lsh(twiceQ, 254)
			boundary.Quo(boundary, c)
			check(boundary)
			check(new(big.Int).Add(boundary, one))
		}
	}
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

func benchQuarticSymbol(b *testing.B, symbol func(*Element) uint8) {
	// Use the same corpus across runs and algorithms. Reusing one input
	// lets the branch predictor learn its entire GCD path.
	rng := rand.New(rand.NewSource(1))
	var inputs [256]Element
	for i := range inputs {
		var buf [Bytes]byte
		_, _ = rng.Read(buf[:])
		inputs[i].SetBytes(buf[:])
	}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		symbol(&inputs[i%len(inputs)])
	}
}

func BenchmarkQuarticSymbol(b *testing.B)    { benchQuarticSymbol(b, (*Element).QuarticSymbol) }
func BenchmarkQuarticSymbolExp(b *testing.B) { benchQuarticSymbol(b, (*Element).QuarticSymbolExp) }

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
