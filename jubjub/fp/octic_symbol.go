package fp

import (
	"math"
	"math/bits"
)

// Joye's I,K,L are the exponents for 1+t, t, and 1+sqrt(2).
// All returned exponents are used modulo 8. Wrapping unsigned arithmetic
// preserves enough low bits for the exact power-of-two divisions.
type octicLaws struct{ I, K, L, d, f, n uint64 }

func octicLawsOf(v [4]uint64) octicLaws {
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

func (a octicLaws) reciprocity(b octicLaws) uint64 {
	return ((a.n - 1) * (b.n - 1) >> 4) + ((b.d*a.f - a.d*b.f) >> 2)
}

// octicWord is a signed two's-complement integer. Arithmetic wraps modulo
// 2^128. The lift and descent below bound the signed results we reconstruct.
type octicWord struct {
	lo uint64
	hi int64
}

func (x octicWord) add(y octicWord) octicWord {
	lo, carry := bits.Add64(x.lo, y.lo, 0)
	return octicWord{lo, x.hi + y.hi + int64(carry)}
}

func (x octicWord) sub(y octicWord) octicWord {
	lo, borrow := bits.Sub64(x.lo, y.lo, 0)
	return octicWord{lo, x.hi - y.hi - int64(borrow)}
}

func (x octicWord) mul(y octicWord) octicWord {
	hi, lo := bits.Mul64(x.lo, y.lo)
	return octicWord{lo, int64(hi + x.lo*uint64(y.hi) + uint64(x.hi)*y.lo)}
}

func (x octicWord) half() octicWord {
	return octicWord{x.lo>>1 | uint64(x.hi)<<63, x.hi >> 1}
}

func (x octicWord) approx() float64 {
	return float64(x.hi+int64(x.lo>>63))*0x1p64 + float64(int64(x.lo))
}

type octicVector [4]octicWord

// Prime above p with t mapping to 7^((p-1)/8).
var octicMod = octicVector{
	{6451692399033376397, 0}, {11914338290826914742, 0},
	{1094064715653685168, 0}, {^uint64(11306383988181076382 - 1), -1},
}

// Magnitudes of round(2^256 * adj(pi)[i] / p), little endian.
// The signs are (-,+,-,-), and pi * adj(pi) = p.
var octicRecip = [4][3]uint64{
	{0x77ea6bb3e629c5d4, 0x57eaddeb671d2e65, 0x850d212844f1dedf},
	{0xbaea9c1d6ca0c4e9, 0x82c817811719611d, 0xf5b4ad4fb1a4715a},
	{0x765508a038de26d9, 0x75e960d7ad66b7d0, 0x169004276828668d},
	{0x53a12b0cdfc4f794, 0x0e4a680a305fee95, 0xe92b0976144bbeca},
}

// Low 128 bits of round(x*c/2^256).
func octicLiftQuotient(x *Element, c [3]uint64) octicWord {
	var p [7]uint64
	for j, v := range c {
		var carry uint64
		for i, u := range x {
			hi, lo := bits.Mul64(u, v)
			lo, k := bits.Add64(lo, p[i+j], 0)
			hi += k
			lo, k = bits.Add64(lo, carry, 0)
			p[i+j] = lo
			carry = hi + k
		}
		p[j+4] = carry
	}
	lo, carry := bits.Add64(p[4], p[3]>>63, 0)
	return octicWord{lo, int64(p[5] + carry)}
}

// b - a*q modulo (t^4+1, 2^128).
func (b octicVector) remainder(a, q octicVector) octicVector {
	return octicVector{
		b[0].sub(a[0].mul(q[0])).add(a[1].mul(q[3])).add(a[2].mul(q[2])).add(a[3].mul(q[1])),
		b[1].sub(a[0].mul(q[1])).sub(a[1].mul(q[0])).add(a[2].mul(q[3])).add(a[3].mul(q[2])),
		b[2].sub(a[0].mul(q[2])).sub(a[1].mul(q[1])).sub(a[2].mul(q[0])).add(a[3].mul(q[3])),
		b[3].sub(a[0].mul(q[3])).sub(a[1].mul(q[2])).sub(a[2].mul(q[1])).sub(a[3].mul(q[0])),
	}
}

func octicLift(x *Element) octicVector {
	var q octicVector
	for i, c := range octicRecip {
		q[i] = octicLiftQuotient(x, c)
		if i != 1 {
			q[i] = (octicWord{}).sub(q[i])
		}
	}
	// A raw Montgomery representative differs from x by the eighth power
	// 2^256. Each quotient coefficient has error <= 1/2 + p/2^257 < 3/4.
	// Thus every remainder coefficient is < 3/4 * sum |pi_i| < 2^65.
	b := octicVector{{x[0], int64(x[1])}}
	return b.remainder(octicMod, q)
}

func (a octicVector) approx() [4]float64 {
	return [4]float64{a[0].approx(), a[1].approx(), a[2].approx(), a[3].approx()}
}

func (a octicVector) narrow() (r [4]int64, ok bool) {
	for i, v := range a {
		n := int64(v.lo)
		if v.hi != n>>63 || n < -(1<<48) || n >= 1<<48 {
			return r, false
		}
		r[i] = n
	}
	return r, true
}

func (a octicVector) low() [4]uint64 {
	return [4]uint64{a[0].lo, a[1].lo, a[2].lo, a[3].lo}
}

func (a octicVector) divRamified() octicVector {
	a[0] = a[0].add(a[1]).sub(a[2]).add(a[3]).half()
	for i := 1; i < 4; i++ {
		a[i] = a[i].sub(a[i-1])
	}
	return a
}

// Multiplying by epsilon^-1 gives a shorter vector than multiplying by
// epsilon when x*(y-w)+z*(y+w) is positive. Both choices give primary
// parity, including when rounding changes the sign of this estimate.
func octicPreferInverse(x, y, z, w float64) bool {
	return x*(y-w)+z*(y+w) >= 0
}

func (a *octicVector) makePrimary() (k, l uint64) {
	var zero octicWord
	for a[0].lo&1 == 0 || a[2].lo&1 != 0 {
		*a = octicVector{a[1], a[2], a[3], zero.sub(a[0])}
		k++
	}
	if a[1].lo&1 != 0 {
		x, y, z, w := a[0], a[1], a[2], a[3]
		s := octicVector{y.sub(w), x.add(z), y.add(w), z.sub(x)}
		if octicPreferInverse(x.approx(), y.approx(), z.approx(), w.approx()) {
			for i := range a {
				a[i] = s[i].sub(a[i])
			}
			l = 1
		} else {
			for i := range a {
				a[i] = s[i].add(a[i])
			}
			l = ^uint64(0)
		}
	}
	if (a[0].lo+a[1].lo+a[2].lo+a[3].lo)&2 != 0 {
		for i := range a {
			a[i] = zero.sub(a[i])
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

var octicRoots = func() (r [8]Element) {
	r[0].SetOne()
	var g Element
	g.SetUint64(7)
	r[1].ExpByOcticExp(g)
	for i := 2; i < len(r); i++ {
		r[i].Mul(&r[i-1], &r[1])
	}
	return
}()

// OcticSymbol returns j in [0,7] such that x^((p-1)/8) = t^j in Fp,
// with t = 7^((p-1)/8). It returns -1 for zero. This variable-time method
// uses Algorithm 1 of Marc Joye (https://marcjoye.github.io/papers/Joy22octic.pdf).
func (x *Element) OcticSymbol() int {
	if x.IsZero() {
		return -1
	}
	if j, ok := octicDescent(octicLift(x)); ok {
		return j
	}
	var r Element
	r.ExpByOcticExp(*x)
	for i := range octicRoots {
		if r == octicRoots[i] {
			return i
		}
	}
	panic("invalid octic character")
}

const octicStepLimit = 512

// a is the nonzero input modulo the prime. The limit also covers a quotient
// estimate which fails to decrease the norm. The caller can always finish
// with the field's fixed-width addition chain.
func octicDescent(a octicVector) (int, bool) {
	b := octicMod
	den := octicLawsOf(b.low())
	var j uint64
	for step := 0; step < octicStepLimit; step++ {
		// |a_i| <= 2^88 before normalization. A primary adjustment costs
		// at most a factor 3, and quotient coefficients are below 2^30.
		// Each four-product remainder therefore fits a signed 128-bit word.
		for _, v := range a {
			if v.hi < -(1<<24) || v.hi >= 1<<24 {
				return 0, false
			}
		}
		if (a[0].lo^a[1].lo^a[2].lo^a[3].lo)&1 == 0 {
			a = a.divRamified()
			j += den.I
			continue
		}
		k, l := a.makePrimary()
		v := a.low()
		j += k*den.K + l*den.L
		if a == (octicVector{{1, 0}}) {
			return int(j & 7), true
		}
		small, ok := octicQuotient(b.approx(), a.approx(), 1<<30)
		if !ok {
			return 0, false
		}
		var q octicVector
		for i, n := range small {
			q[i] = octicWord{uint64(n), n >> 63}
		}
		r := b.remainder(a, q)
		// Adjacent operands remain coprime. A zero remainder certifies a
		// primary unit; its residual exponents are a0-1 and a3 modulo 8.
		if r == (octicVector{}) {
			return int((j + (v[0]-1)*den.K + v[3]*den.L) & 7), true
		}
		num := octicLawsOf(v)
		j += num.reciprocity(den)
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
		num := octicLawsOf([4]uint64{uint64(x), uint64(y), uint64(z), uint64(w)})
		j += num.reciprocity(den)
		den = num
	}
	return 0, false
}
