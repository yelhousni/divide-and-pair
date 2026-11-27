package fp

import (
	"math/bits"

	"github.com/yelhousni/divide-and-pair/internal/quartic"
)

const (
	// pi = A - Bi is primary, with norm p = 2^255 - 19.
	piRe0 uint64 = 0x3feab578735893c3
	piRe1 uint64 = 0x33a5cbdded73544f
	piIm0 uint64 = 0x43c900683eb6254a
	piIm1 uint64 = 0xad7eb9766c0b7b36
)

var sqrtMinusOneFp, minusOneFp Element

func init() {
	sqrtMinusOneFp.SetString("19681161376707505956807079304988542015446066515923890162744021073123829784752")
	minusOneFp.SetOne()
	minusOneFp.Neg(&minusOneFp)
}

// roundProduct returns round(z*c/2^255). Both inputs are unsigned, z < p
// and c < 2^128. Only the two high words of the rounded quotient are needed.
func roundProduct(z *Element, c0, c1 uint64) (uint64, uint64) {
	var product [6]uint64
	for j, c := range [2]uint64{c0, c1} {
		var carry uint64
		for i, v := range z {
			hi, lo := bits.Mul64(v, c)
			lo, k := bits.Add64(lo, product[i+j], 0)
			hi += k
			lo, k = bits.Add64(lo, carry, 0)
			product[i+j] = lo
			carry = hi + k
		}
		product[j+4] = carry
	}
	lo := product[3]>>63 | product[4]<<1
	hi := product[4]>>63 | product[5]<<1
	lo, carry := bits.Add64(lo, (product[3]>>62)&1, 0)
	return lo, hi + carry
}

// quarticInt is a two's complement Gaussian coordinate with a guard word for
// subtraction and division by 1+i.
// Normalized coordinates have magnitude < 2^128.
type quarticInt struct {
	lo, mid uint64
	hi      int64
}

func (a quarticInt) add(b quarticInt) quarticInt {
	lo, c := bits.Add64(a.lo, b.lo, 0)
	mid, c := bits.Add64(a.mid, b.mid, c)
	return quarticInt{lo, mid, a.hi + b.hi + int64(c)}
}

func (a quarticInt) sub(b quarticInt) quarticInt {
	lo, c := bits.Sub64(a.lo, b.lo, 0)
	mid, c := bits.Sub64(a.mid, b.mid, c)
	return quarticInt{lo, mid, a.hi - b.hi - int64(c)}
}

func (a quarticInt) neg() quarticInt { return (quarticInt{}).sub(a) }

func (a quarticInt) shift(k int) quarticInt {
	if k >= 64 {
		if k >= 128 {
			return quarticInt{uint64(a.hi >> (k - 128)), uint64(a.hi >> 63), a.hi >> 63}
		}
		a = quarticInt{a.mid, uint64(a.hi), a.hi >> 63}
		k -= 64
	}
	return quarticInt{a.lo>>k | a.mid<<(64-k), a.mid>>k | uint64(a.hi)<<(64-k), a.hi >> k}
}

// lowProduct computes a product modulo 2^128.
func lowProduct(a0, a1, b0, b1 uint64) quarticInt {
	hi, lo := bits.Mul64(a0, b0)
	return quarticInt{lo, hi + a0*b1 + a1*b0, 0}
}

type gaussian struct{ r, i quarticInt }

// Reduce the Montgomery representative modulo pi. Its factor R = 2^256
// is a fourth power. Replacing p by 2^255 in each rounded quotient adds
// < 19*2^128/2^255 to its rounding error. Consequently each remainder
// coordinate is bounded by (A+B)*(1/2 + 19*2^-127) < 2^127.
func quarticReduce(z *Element) gaussian {
	q0, q1 := roundProduct(z, piRe0, piRe1)
	t0, t1 := roundProduct(z, piIm0, piIm1)
	// Subtract (q + ti)(A - Bi) from z, keeping the low 128 bits.
	r := quarticInt{z[0], z[1], 0}.sub(lowProduct(q0, q1, piRe0, piRe1))
	r = r.sub(lowProduct(t0, t1, piIm0, piIm1))
	i := lowProduct(q0, q1, piIm0, piIm1).sub(lowProduct(t0, t1, piRe0, piRe1))
	r.hi = int64(r.mid) >> 63
	i.hi = int64(i.mid) >> 63
	return gaussian{r, i}
}

// Strip 2^k = (-i)^k (1+i)^(2k), then at most one more 1+i.
// Rotate the remaining odd element by i^rotation to make it primary.
func (g gaussian) primary(f gaussian) (gaussian, uint64) {
	k := bits.TrailingZeros64(g.r.lo | g.i.lo)
	if k == 64 {
		k += bits.TrailingZeros64(g.r.mid | g.i.mid)
	}
	r, i := g.r.shift(k), g.i.shift(k)
	unit, ramified := quartic.Supp(f.r.lo, f.i.lo)
	e := uint64(k) * (2*ramified - unit)
	if r.lo&i.lo&1 != 0 {
		r, i = r.add(i).shift(1), i.sub(r).shift(1)
		e += ramified
	}
	if i.lo&1 != 0 {
		r, i = i.neg(), r
		e -= unit
	}
	if (r.lo+i.lo)&2 != 0 {
		r, i = r.neg(), i.neg()
		e -= 2 * unit
	}
	return gaussian{r, i}, e
}

// Primary inputs are both 1 modulo (1+i)^3. The 1-norm orders them
// without squaring either coordinate. For the larger input g,
// removing those three factors after subtraction gives
// |g'| <= (1+sqrt(2))*|g|/sqrt(8) < |g| in complex absolute value.
func (g gaussian) size() quarticInt {
	r, i := g.r, g.i
	if r.hi < 0 {
		r = r.neg()
	}
	if i.hi < 0 {
		i = i.neg()
	}
	return r.add(i)
}

// QuarticSymbol computes z^((p-1)/4) as an exponent of i (0, 1, 2, or 3)
// using a variable-time Gaussian binary GCD. It returns 0 on zero by convention.
func (z *Element) QuarticSymbol() uint8 {
	if z.IsZero() {
		return 0
	}
	f := gaussian{quarticInt{piRe0, piRe1, 0}, quarticInt{piIm0, piIm1, 0}.neg()}
	g, e := quarticReduce(z).primary(f)
	// The running character is i^e times (g/f)_4. The operands are coprime
	// since pi has prime norm and z is nonzero modulo pi. Equal primary
	// operands can therefore only be 1. Their 1-norms fit in 128 bits
	// since the initial complex magnitudes are at most sqrt(p).
	sf, sg := f.size(), g.size()
	for {
		if sf.mid|sg.mid == 0 && (sf.lo|sg.lo) < 1<<61 {
			return quartic.Word(int64(f.r.lo), int64(f.i.lo), int64(g.r.lo), int64(g.i.lo), e)
		}
		if sf.mid == 0 && sf.lo == 1 || sg.mid == 0 && sg.lo == 1 {
			return uint8(e & 3)
		}
		if sf.mid > sg.mid || sf.mid == sg.mid && sf.lo > sg.lo {
			// Biquadratic reciprocity contributes -1 iff both real parts are 3 mod 4.
			e += f.r.lo & g.r.lo & 2
			f, g = g, f
			sf = sg
		}
		g = gaussian{g.r.sub(f.r), g.i.sub(f.i)}
		var update uint64
		g, update = g.primary(f)
		e += update
		sg = g.size()
	}
}

// QuarticSymbolExp computes the quartic symbol via addition chain exponentiation.
// χ₄(z) = z·(z^((p-5)/8))². Returns 0,1,2,3.
func (z *Element) QuarticSymbolExp() uint8 {
	if z.IsZero() {
		return 0
	}
	var t Element
	t.ExpBySqrtPm5o8(*z)
	var result Element
	result.Square(&t)
	result.Mul(&result, z)
	if result.IsOne() {
		return 0
	}
	if result.Equal(&sqrtMinusOneFp) {
		return 1
	}
	if result.Equal(&minusOneFp) {
		return 2
	}
	return 3
}
