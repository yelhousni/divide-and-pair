package fp

import (
	"filippo.io/edwards25519/field"
)

// ToFilippo converts a Montgomery-form Element to a filippo field.Element
// (5×51-bit unsaturated representation).
func (z *Element) ToFilippo() field.Element {
	b := z.Bytes()
	// Bytes() returns big-endian; filippo expects little-endian
	var le [32]byte
	for i := range 32 {
		le[i] = b[31-i]
	}
	var fe field.Element
	fe.SetBytes(le[:])
	return fe
}

// SetFilippo sets z from a filippo field.Element.
func (z *Element) SetFilippo(fe *field.Element) *Element {
	le := fe.Bytes()
	// Bytes() returns little-endian; SetBytes expects big-endian
	var be [Bytes]byte
	for i := range 32 {
		be[i] = le[31-i]
	}
	z.SetBytes(be[:])
	return z
}

// SqrtFilippo computes the square root using filippo's field arithmetic
// (5×51-bit unsaturated representation with dedicated squaring).
// Returns nil if x is not a QR.
func (z *Element) SqrtFilippo(x *Element) *Element {
	fe := x.ToFilippo()

	// For p = 2^255-19 (p ≡ 5 mod 8):
	// candidate = x^((p+3)/8) = x * x^((p-5)/8) = x * Pow22523(x)
	var p58, candidate, check field.Element
	p58.Pow22523(&fe)             // x^((p-5)/8)
	candidate.Multiply(&fe, &p58) // x^((p+3)/8)

	// Verify: candidate² == x or candidate² == -x
	check.Square(&candidate)
	if check.Equal(&fe) == 1 {
		z.SetFilippo(&candidate)
		return z
	}

	// Try candidate * sqrt(-1)
	var sqrtM1 field.Element
	sqrtM1.SetBytes([]byte{
		0xb0, 0xa0, 0x0e, 0x4a, 0x27, 0x1b, 0xee, 0xc4,
		0x78, 0xe4, 0x2f, 0xad, 0x06, 0x18, 0x43, 0x2f,
		0xa7, 0xd7, 0xfb, 0x3d, 0x99, 0x00, 0x4d, 0x2b,
		0x0b, 0xdf, 0xc1, 0x4f, 0x80, 0x24, 0x83, 0x2b,
	})
	candidate.Multiply(&candidate, &sqrtM1)
	check.Square(&candidate)
	if check.Equal(&fe) == 1 {
		z.SetFilippo(&candidate)
		return z
	}

	return nil
}

// LegendreFilippo computes the Legendre symbol using filippo's field.
func (z *Element) LegendreFilippo() int {
	fe := z.ToFilippo()

	// Legendre = x^((p-1)/2). Since p = 2^255-19:
	// (p-1)/2 = 2^254 - 10
	// We compute via: x^((p-5)/8)^4 * x^3
	// Since x^((p-1)/2) = (x^((p-5)/8))^4 * x^3
	// Proof: 4*(p-5)/8 + 3 = (p-5)/2 + 3 = (p-5+6)/2 = (p+1)/2
	// Hmm that gives (p+1)/2 not (p-1)/2. Let me redo.
	// x^((p-1)/2) = x^(2^254 - 10)
	// Use: x^((p-5)/8) = Pow22523(x). Exponent = (p-5)/8 = 2^252 - 3.
	// x^(2^254 - 10) = x^(4*(2^252-3) + 2) = (x^((p-5)/8))^4 * x^2
	var p58, t field.Element
	p58.Pow22523(&fe)   // x^((p-5)/8) = x^(2^252-3)
	t.Square(&p58)      // x^(2*(2^252-3))
	t.Square(&t)        // x^(4*(2^252-3)) = x^(2^254-12)
	t.Multiply(&t, &fe) // x^(2^254-11)
	t.Multiply(&t, &fe) // x^(2^254-10) = x^((p-1)/2)

	var one field.Element
	one.One()
	if t.Equal(&one) == 1 {
		return 1
	}

	var zero field.Element
	zero.Zero()
	if t.Equal(&zero) == 1 {
		return 0
	}
	return -1
}

// QuarticSymbolExpFilippo computes the quartic symbol using filippo's field.
// χ₄(z) = z·(z^((p-5)/8))². Returns 0,1,2,3.
func (z *Element) QuarticSymbolExpFilippo() uint8 {
	if z.IsZero() {
		return 0
	}

	fe := z.ToFilippo()

	var p58, t, result field.Element
	p58.Pow22523(&fe)
	t.Square(&p58)
	result.Multiply(&fe, &t)

	var one field.Element
	one.One()
	if result.Equal(&one) == 1 {
		return 0
	}

	sqrtM1 := sqrtMinusOneFp.ToFilippo()
	if result.Equal(&sqrtM1) == 1 {
		return 1
	}

	var minusOne field.Element
	minusOne.One()
	minusOne.Negate(&minusOne)
	if result.Equal(&minusOne) == 1 {
		return 2
	}

	return 3
}
