package octic

import (
	"math"
	"math/bits"
)

// Joye's I,K,L are the exponents for 1+t, t, and 1+sqrt(2).
// All returned exponents are used modulo 8. Wrapping unsigned arithmetic
// preserves enough low bits for the exact power-of-two divisions.
type octicLaws struct{ I, K, L, d, f, n uint64 }

func LawsOf(v [4]uint64) octicLaws {
	x, y, z, t := v[0], v[1], v[2], v[3]
	a := x*x - z*z + 2*y*t
	b := t*t - y*y + 2*x*z
	d := x*y + x*t - y*z + z*t
	f := x*y - x*t + y*z + z*t
	dd := d * d
	return octicLaws{
		I: (a - 1 + b + 6*d + b*b + 2*b*d + (dd * dd >> 1)) >> 3,
		K: (a - 1 + 4*b + 2*b*d + 2*dd) >> 2,
		L: (d - 3*b - b*d - 2*dd) >> 1,
		d: d, f: f, n: a*a + b*b,
	}
}

func (a octicLaws) Reciprocity(b octicLaws) uint64 {
	return ((a.n - 1) * (b.n - 1) >> 4) + ((b.d*a.f - a.d*b.f) >> 2)
}

// Word is a signed two's-complement integer. Arithmetic wraps modulo
// 2^128. The lift and descent below bound the signed results we reconstruct.
type Word struct {
	Lo uint64
	Hi int64
}

func (x Word) add(y Word) Word {
	lo, carry := bits.Add64(x.Lo, y.Lo, 0)
	return Word{lo, x.Hi + y.Hi + int64(carry)}
}

func (x Word) Sub(y Word) Word {
	lo, borrow := bits.Sub64(x.Lo, y.Lo, 0)
	return Word{lo, x.Hi - y.Hi - int64(borrow)}
}

func (x Word) mul(y Word) Word {
	hi, lo := bits.Mul64(x.Lo, y.Lo)
	return Word{lo, int64(hi + x.Lo*uint64(y.Hi) + uint64(x.Hi)*y.Lo)}
}

func (x Word) half() Word {
	return Word{x.Lo>>1 | uint64(x.Hi)<<63, x.Hi >> 1}
}

func (x Word) approx() float64 {
	return float64(x.Hi+int64(x.Lo>>63))*0x1p64 + float64(int64(x.Lo))
}

// Vector holds coefficients in ascending powers of a primitive eighth root.
type Vector [4]Word

// b - a*q modulo (t^4+1, 2^128).
func (b Vector) Remainder(a, q Vector) Vector {
	return Vector{
		b[0].Sub(a[0].mul(q[0])).add(a[1].mul(q[3])).add(a[2].mul(q[2])).add(a[3].mul(q[1])),
		b[1].Sub(a[0].mul(q[1])).Sub(a[1].mul(q[0])).add(a[2].mul(q[3])).add(a[3].mul(q[2])),
		b[2].Sub(a[0].mul(q[2])).Sub(a[1].mul(q[1])).Sub(a[2].mul(q[0])).add(a[3].mul(q[3])),
		b[3].Sub(a[0].mul(q[3])).Sub(a[1].mul(q[2])).Sub(a[2].mul(q[1])).Sub(a[3].mul(q[0])),
	}
}

func (a Vector) approx() [4]float64 {
	return [4]float64{a[0].approx(), a[1].approx(), a[2].approx(), a[3].approx()}
}

func (a Vector) narrow() (r [4]int64, ok bool) {
	for i, v := range a {
		n := int64(v.Lo)
		if v.Hi != n>>63 || n < -(1<<48) || n >= 1<<48 {
			return r, false
		}
		r[i] = n
	}
	return r, true
}

func (a Vector) low() [4]uint64 {
	return [4]uint64{a[0].Lo, a[1].Lo, a[2].Lo, a[3].Lo}
}

func (a Vector) divRamified() Vector {
	a[0] = a[0].add(a[1]).Sub(a[2]).add(a[3]).half()
	for i := 1; i < 4; i++ {
		a[i] = a[i].Sub(a[i-1])
	}
	return a
}

// Multiplying by epsilon^-1 gives a shorter vector than multiplying by
// epsilon when x*(y-w)+z*(y+w) is positive. Both choices give primary
// parity, including when rounding changes the sign of this estimate.
func octicPreferInverse(x, y, z, w float64) bool {
	return x*(y-w)+z*(y+w) >= 0
}

func (a *Vector) makePrimary() (k, l uint64) {
	var zero Word
	for a[0].Lo&1 == 0 || a[2].Lo&1 != 0 {
		*a = Vector{a[1], a[2], a[3], zero.Sub(a[0])}
		k++
	}
	if a[1].Lo&1 != 0 {
		x, y, z, w := a[0], a[1], a[2], a[3]
		s := Vector{y.Sub(w), x.add(z), y.add(w), z.Sub(x)}
		if octicPreferInverse(x.approx(), y.approx(), z.approx(), w.approx()) {
			for i := range a {
				a[i] = s[i].Sub(a[i])
			}
			l = 1
		} else {
			for i := range a {
				a[i] = s[i].add(a[i])
			}
			l = ^uint64(0)
		}
	}
	if (a[0].Lo+a[1].Lo+a[2].Lo+a[3].Lo)&2 != 0 {
		for i := range a {
			a[i] = zero.Sub(a[i])
		}
		k += 4
	}
	return
}

// Approximate b/a in the power basis via the relative Gaussian norm.
// Rounding only selects an integral multiple to subtract. Reciprocity uses
// the exact integer remainder even if rounding picks a neighboring quotient.
func octicQuotient(b, a [4]float64, limit float64) (q [4]int64, ok bool) {
	x, y, z, w := a[0], a[1], a[2], a[3]
	A, B := x*x-z*z+2*y*w, w*w-y*y+2*x*z
	inv := 1 / (A*A + B*B)
	u, v, s, t := (x*A+z*B)*inv, (-y*A-w*B)*inv, (z*A-x*B)*inv, (y*B-w*A)*inv
	x, y, z, w = b[0], b[1], b[2], b[3]
	values := [4]float64{
		x*u - y*t - z*s - w*v,
		x*v + y*u - z*t - w*s,
		x*s + y*v + z*u - w*t,
		x*t + y*s + z*v + w*u,
	}
	for i, v := range values {
		n := math.Round(v)
		if !(math.Abs(n) < limit) {
			return q, false
		}
		q[i] = int64(n)
	}
	return q, true
}

const octicStepLimit = 512

// Descent computes the octic index for a nonzero residue a modulo a primary
// prime b. Both inputs must have coordinate magnitudes at most 2^88.
// A false result signals that a size or iteration bound was reached.
func Descent(a, b Vector) (int, bool) {
	den := LawsOf(b.low())
	var j uint64
	for step := 0; step < octicStepLimit; step++ {
		// |a_i| <= 2^88 before normalization. A primary adjustment costs
		// at most a factor 3, and quotient coefficients are below 2^30.
		// Each four-product remainder therefore fits a signed 128-bit word.
		for _, v := range a {
			if v.Hi < -(1<<24) || v.Hi >= 1<<24 {
				return 0, false
			}
		}
		if (a[0].Lo^a[1].Lo^a[2].Lo^a[3].Lo)&1 == 0 {
			a = a.divRamified()
			j += den.I
			continue
		}
		k, l := a.makePrimary()
		v := a.low()
		j += k*den.K + l*den.L
		if a == (Vector{{1, 0}}) {
			return int(j & 7), true
		}
		small, ok := octicQuotient(b.approx(), a.approx(), 1<<30)
		if !ok {
			return 0, false
		}
		var q Vector
		for i, n := range small {
			q[i] = Word{uint64(n), n >> 63}
		}
		r := b.Remainder(a, q)
		// Adjacent operands remain coprime. A zero remainder certifies a
		// primary unit; its residual exponents are a0-1 and a3 modulo 8.
		if r == (Vector{}) {
			return int((j + (v[0]-1)*den.K + v[3]*den.L) & 7), true
		}
		num := LawsOf(v)
		j += num.Reciprocity(den)
		den = num
		if aa, ok := r.narrow(); ok {
			if bb, ok := a.narrow(); ok {
				return octicDescent64(aa, bb, den, j, step+1)
			}
		}
		a, b = r, a
	}
	return 0, false
}

// octicWidth bounds the largest absolute coordinate by 2^width. A negative
// coordinate contributes the bits of its complement, including -1 with width 0.
func octicWidth(x, y, z, w int64) int {
	return bits.Len64(uint64((x ^ (x >> 63)) | (y ^ (y >> 63)) | (z ^ (z >> 63)) | (w ^ (w >> 63))))
}

// Enter once all coordinates have magnitude at most 2^48. Before each
// product, the normalized divisor sets a quotient limit. Each product fits
// in 60 magnitude bits, and their sum with b fits in int64.
func octicDescent64(a, b [4]int64, den octicLaws, j uint64, step int) (int, bool) {
	for ; step < octicStepLimit; step++ {
		x, y, z, w := a[0], a[1], a[2], a[3]
		if octicWidth(x, y, z, w) > 56 {
			return 0, false
		}
		if (x^y^z^w)&1 == 0 {
			x = (x + y - z + w) >> 1
			y -= x
			z -= y
			w -= z
			a = [4]int64{x, y, z, w}
			j += den.I
			continue
		}
		for x&1 == 0 || z&1 != 0 {
			x, y, z, w = y, z, w, -x
			j += den.K
		}
		if y&1 != 0 {
			if octicPreferInverse(float64(x), float64(y), float64(z), float64(w)) {
				x, y, z, w = y-w-x, x+z-y, y+w-z, z-x-w
				j += den.L
			} else {
				x, y, z, w = y-w+x, x+z+y, y+w+z, z-x+w
				j -= den.L
			}
		}
		if (x+y+z+w)&2 != 0 {
			x, y, z, w = -x, -y, -z, -w
			j += 4 * den.K
		}
		if x == 1 && y == 0 && z == 0 && w == 0 {
			return int(j & 7), true
		}
		width := octicWidth(x, y, z, w)
		limit := int64(1) << min(30, 60-width)
		u, v, s, t := b[0], b[1], b[2], b[3]
		q, ok := octicQuotient(
			[4]float64{float64(u), float64(v), float64(s), float64(t)},
			[4]float64{float64(x), float64(y), float64(z), float64(w)}, float64(limit))
		if !ok {
			return 0, false
		}
		h, i, k, l := q[0], q[1], q[2], q[3]
		a = [4]int64{
			u - x*h + y*l + z*k + w*i,
			v - x*i - y*h + z*l + w*k,
			s - x*k - y*i - z*h + w*l,
			t - x*l - y*k - z*i - w*h,
		}
		if a == ([4]int64{}) {
			return int((j + (uint64(x)-1)*den.K + uint64(w)*den.L) & 7), true
		}
		b = [4]int64{x, y, z, w}
		num := LawsOf([4]uint64{uint64(x), uint64(y), uint64(z), uint64(w)})
		j += num.Reciprocity(den)
		den = num
	}
	return 0, false
}
