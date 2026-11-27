package fp

import "math/bits"

// Caranay–Scheidler, An Efficient Seventh Power Residue Symbol Algorithm,
// Algorithms 5.3, 7.1–7.2 and 8.1 (doi:10.1142/S1793042110003770).
// Represent sum a[i]*zeta^i with a[0]=0. Subtracting the same coefficient
// from all seven entries preserves the element since Phi_7(zeta)=0.
type septicInteger [7]int64

func (a septicInteger) canonical() septicInteger {
	c := a[0]
	for i := range a {
		a[i] -= c
	}
	return a
}

func (a septicInteger) mul(b septicInteger) septicInteger {
	// Cyclic convolution, followed by reduction by Phi_7.
	return (septicInteger{
		a[1]*b[6] + a[2]*b[5] + a[3]*b[4] + a[4]*b[3] + a[5]*b[2] + a[6]*b[1],
		a[2]*b[6] + a[3]*b[5] + a[4]*b[4] + a[5]*b[3] + a[6]*b[2],
		a[1]*b[1] + a[3]*b[6] + a[4]*b[5] + a[5]*b[4] + a[6]*b[3],
		a[1]*b[2] + a[2]*b[1] + a[4]*b[6] + a[5]*b[5] + a[6]*b[4],
		a[1]*b[3] + a[2]*b[2] + a[3]*b[1] + a[5]*b[6] + a[6]*b[5],
		a[1]*b[4] + a[2]*b[3] + a[3]*b[2] + a[4]*b[1] + a[6]*b[6],
		a[1]*b[5] + a[2]*b[4] + a[3]*b[3] + a[4]*b[2] + a[5]*b[1],
	}).canonical()
}

func (a septicInteger) conjugate(k int) (c septicInteger) {
	for i := 1; i < 7; i++ {
		c[(i*k)%7] = a[i]
	}
	return
}

func (a septicInteger) rotate(k int) (c septicInteger) {
	for i := range a {
		c[(i+k)%7] = a[i]
	}
	return c.canonical()
}

func septicMod7(x int64) int64 {
	x %= 7
	return x + (x >> 63 & 7)
}

var septicInv = [7]int64{0, 1, 4, 5, 2, 3, 6}

// pi generates the prime (p, zeta-root), for the root in septicRoots.
// The signed reciprocal coefficients are round(2^128 * (p/pi)[i] / p).
var septicPrime = septicInteger{0, 71954, 2030085, 837349, 1616928, 2162954, 354058}
var septicRecip = [...]struct {
	lo, hi   uint64
	negative bool
}{
	{0x2d0d5043cddf2592, 0xfdb5af6935, true},
	{0x7800bac0f191df40, 0x47dc41dc080, false},
	{0x500f0a703afa91da, 0x25b93965cb3, false},
	{0x43974665647c5cd8, 0x9a82d8de10, true},
	{0xcbc76168819d0d9c, 0x529995e36d0, false},
	{0xf174ae3591caab1a, 0x1deeafc639d, true},
}

// Low word of round(x*c/2^128). Higher quotient words disappear in the
// lift because its subtraction is carried out modulo 2^64.
func septicLiftWord(x [2]uint64, lo, hi uint64) uint64 {
	h00, _ := bits.Mul64(x[0], lo)
	h01, l01 := bits.Mul64(x[0], hi)
	h10, l10 := bits.Mul64(x[1], lo)
	mid, c1 := bits.Add64(h00, l01, 0)
	mid, c2 := bits.Add64(mid, l10, 0)
	return x[1]*hi + h01 + h10 + c1 + c2 + (mid >> 63)
}

func septicLift(x [2]uint64) septicInteger {
	var q septicInteger
	for i, c := range septicRecip {
		v := int64(septicLiftWord(x, c.lo, c.hi))
		if c.negative {
			v = -v
		}
		q[i+1] = v
	}
	r := q.mul(septicPrime)
	// Each quotient error is <= 1/2+p/2^129 < 3/4. Canonicalizing
	// the cyclic product gives |r_i| <= 2*(3/4)*sum |pi_i| < 2^24
	// for the exact remainder. Its signed low words therefore suffice.
	for i := 1; i < 7; i++ {
		r[i] = -int64(x[0]) - r[i]
	}
	return r
}

// eta1=zeta+zeta^-1 and eta2=zeta^2+zeta^-2. Only their products
// modulo 49 enter the complementary laws.
var septicUnits = func() (units [7][7]septicInteger) {
	eta1 := septicInteger{0, 1, 0, 0, 0, 0, 1}
	eta2 := septicInteger{0, 0, 1, 0, 0, 1, 0}
	units[0][0] = septicInteger{0, -1, -1, -1, -1, -1, -1}
	for m := range units {
		if m > 0 {
			units[m][0] = units[m-1][0].mul(eta1)
		}
		for n := 1; n < 7; n++ {
			units[m][n] = units[m][n-1].mul(eta2)
		}
	}
	return
}()

type septicUnit struct{ k, m, n int64 }
type septicPrimary struct {
	unit  septicUnit // u such that u*a is primary
	laws  septicUnit // indices of zeta, eta1, eta2 modulo u*a
	omega int64      // index of 1-zeta modulo u*a
}

func (u septicUnit) index(laws septicUnit) int64 {
	return u.k*laws.k + u.m*laws.m + u.n*laws.n
}

func (a septicInteger) sum() int64     { return a[1] + a[2] + a[3] + a[4] + a[5] + a[6] }
func (a septicInteger) moment1() int64 { return a[1] + 2*a[2] + 3*a[3] + 4*a[4] + 5*a[5] + 6*a[6] }

// primary computes Algorithm 5.3 and Theorem 6.3 on small residues alone.
// Its input is coprime to 1-zeta.
func (a septicInteger) primary() septicPrimary {
	for i := 1; i < 7; i++ {
		a[i] %= 49
	}
	inv := septicInv[septicMod7(a.sum())]
	k := septicMod7(-a.moment1() * inv)
	a = a.rotate(int(k))
	h := septicMod7((a[1] - 3*a[2] + 2*a[3] + 2*a[4] - 3*a[5] + a[6]) * inv)
	e := septicMod7((a[1] + 2*a[2] - 3*a[3] - 3*a[4] + 2*a[5] + a[6]) * inv)
	m, n := septicMod7(h+e-3*h*h), septicMod7(3*h-2*e-h*h)
	a = a.mul(septicUnits[m][n])
	for i := 1; i < 7; i++ {
		a[i] %= 49
	}
	inv = septicInv[septicMod7(a.sum())]
	d := septicMod7((a[1] + a[2] - a[3] + a[4] - a[5] - a[6]) * inv)
	f := septicMod7((a[1] - 3*a[2] - 2*a[3] + 2*a[4] + 3*a[5] - a[6]) * inv)

	// a*conj(a) = U + V*eta1 + W*eta2, with V,W divisible by 7
	// by primarity. Its real cubic norm modulo 49 keeps only constant
	// and linear terms; both real units have trace -1.
	c := a.mul(a.conjugate(6))
	U, V, W := -c[3]%49, (c[1]-c[3])%49, (c[2]-c[3])%49
	norm := (U*U*(U-V-W)%49 + 49) % 49
	zeta := (norm - 1) / 7
	return septicPrimary{
		unit:  septicUnit{k, m, n},
		laws:  septicUnit{zeta, septicMod7(f - 2*d), septicMod7(3 * (d - f))},
		omega: septicMod7(4*(zeta+1) - (a.moment1()/7)*inv),
	}
}

func (a septicInteger) width() int {
	var mask uint64
	for _, v := range a {
		mask |= uint64(v ^ (v >> 63))
	}
	return bits.Len64(mask)
}

// Round all coordinates with a common power-of-two scale. This is used only
// to estimate the inverse for division, never for the Euclidean operands.
func (a septicInteger) scaled() septicInteger {
	if shift := a.width() - 26; shift > 0 {
		half := int64(1) << (shift - 1)
		for i := 1; i < 7; i++ {
			a[i] = (a[i] + half) >> shift
		}
	}
	return a
}

// Balance the three complex absolute values using the real units. This
// keeps the integer approximation of the inverse well conditioned. The
// squared Euclidean length is Tr(a*conj(a)); choose a unit only if it
// decreases that length. Account for every chosen unit in the symbol.
func (a septicInteger) balanced(laws septicUnit) (septicInteger, int64) {
	var correction int64
	indices := [4]int64{0, laws.m, laws.n, -laws.m - laws.n}
	for range 16 {
		// Use at most 25 magnitude bits for the trace comparisons.
		s := a
		if shift := a.width() - 25; shift > 0 {
			for i := 1; i < 7; i++ {
				s[i] >>= shift
			}
		}
		sum := s.sum()
		var tr [4]int64
		for j := range tr {
			var dot int64
			for i := 1; i < 7; i++ {
				dot += s[i] * s[(i+j)%7]
			}
			tr[j] = 7*dot - sum*sum
		}
		energy := tr[0]
		best, choice, inverse := energy, 0, false
		for j := 1; j <= 3; j++ {
			j2 := min(2*j, 7-2*j)
			// eta_j^2=2+eta_2j and eta_j^-2=1-eta_j-2*eta_2j.
			if v := 2*energy + 2*tr[j2]; v < best {
				best, choice, inverse = v, j, false
			}
			if v := energy - 2*tr[j] - 4*tr[j2]; v < best {
				best, choice, inverse = v, j, true
			}
		}
		if choice == 0 {
			break
		}
		j := choice
		var next septicInteger
		for i := range next {
			next[i] = a[(i+j)%7] + a[(i+7-j)%7]
			if inverse {
				next[i] += a[(i+2*j)%7] + a[(i+7-2*j)%7]
			}
		}
		// eta_j^-1=eta_j+eta_2j. Remove the multiplier from the
		// numerator's symbol to keep the running index invariant.
		if inverse {
			correction += indices[j]
		} else {
			correction -= indices[j]
		}
		a = next.canonical()
		if a.width() > 28 {
			break
		}
	}
	return a, correction
}

// Lenstra's quotient selection on a fixed-point approximation to b/a.
// With operand magnitudes <= 2^28 and scaled magnitudes <= 2^26, each
// canonical convolution is bounded by 12 products and fits in int64.
func septicRemainder(b, a septicInteger) (septicInteger, bool) {
	bar := a.conjugate(6)
	real := a.mul(bar).scaled()
	h := real.conjugate(2).mul(real.conjugate(4)).scaled()
	adj := bar.mul(h).scaled()
	// Project a*adj onto Q using trace/6. The scaled product of the five
	// conjugates approximates a scalar multiple of 1/a. The common scale
	// cancels in b*adj / den.
	den := -a.mul(adj).sum() / 6
	if den <= 0 {
		return septicInteger{}, false
	}
	num := b.mul(adj)
	var q, frac septicInteger
	const scale = int64(1) << 16
	shift := max(bits.Len64(uint64(den))-45, 0)
	for i, v := range num {
		n, r := v/den, v%den
		if r < 0 {
			n--
			r += den
		}
		if n < -(1<<28) || n >= 1<<28 {
			return septicInteger{}, false
		}
		q[i] = n
		frac[i] = ((r >> shift) * scale) / (den >> shift)
	}
	// Adding one to the smallest fractional coordinate explores the seven
	// translates in Algorithm 7.1. Trace(z*conj(z)) = 7 sum z_i²-(sum z_i)².
	for range 7 {
		var sum, squares int64
		smallest := 0
		for i, v := range frac {
			sum += v
			squares += v * v
			if v < frac[smallest] {
				smallest = i
			}
		}
		if 7*squares-sum*sum <= 4*scale*scale {
			break
		}
		frac[smallest] += scale
		q[smallest]--
	}
	// Any integral q preserves the residue class. Approximation can affect
	// progress, while the remainder and all symbol updates are exact.
	r := a.mul(q.canonical())
	for i := 1; i < 7; i++ {
		r[i] = b[i] - r[i]
	}
	return r, true
}

var septicPrimePrimary = septicPrime.primary()

// Descend on the small unnormalized operands. If A=u*a and B=v*b are
// primary, reciprocity replaces (a/B) by (b/A) with corrections (v/A)
// and (u/B)^-1. Thus the unit corrections need only arithmetic modulo 49.
func septicDescent(a, b septicInteger, den septicPrimary) (int, bool) {
	var index int64
	for range 256 {
		if a.width() > 28 {
			return 0, false
		}
		sum := a.sum()
		if sum%7 == 0 {
			// For a=(1-zeta)*r, sum(a)=7*r[6] and
			// r[i]-r[i-1]=a[i]-r[6] in our canonical basis.
			last := sum / 7
			for i := 1; i < 7; i++ {
				a[i] += a[i-1] - last
			}
			index += den.omega
			continue
		}
		var adjust int64
		a, adjust = a.balanced(den.laws)
		if a.width() > 28 {
			return 0, false
		}
		index += adjust
		num := a.primary()
		index -= num.unit.index(den.laws)
		r, ok := septicRemainder(b, a)
		if !ok {
			return 0, false
		}
		if r == (septicInteger{}) {
			// Coprimality makes a a unit. Its primary associate is a
			// seventh power up to sign, so its residual character is one.
			return int(septicMod7(index)), true
		}
		index += den.unit.index(num.laws)
		a, b, den = r, a, num
	}
	return 0, false
}

var septicRoots = func() (r [7]Element) {
	r[0].SetOne()
	r[1] = Element{0xd055abe62be39447, 0x5d25498c8b778e76}
	for i := 2; i < 7; i++ {
		r[i].Mul(&r[i-1], &r[1])
	}
	return
}()

// SepticSymbol returns the index j in [0,6] such that x^((p-1)/7)=root^j,
// where root=9^((p-1)/7). Zero returns -1.
// The variable-time Euclidean descent uses bounded int64 coefficients. If a
// size or step bound is reached, the fixed-width addition chain finishes it.
func (x *Element) SepticSymbol() int {
	if x.IsZero() {
		return -1
	}
	if j, ok := septicDescent(septicLift(x.Bits()), septicPrime, septicPrimePrimary); ok {
		return j
	}
	var y Element
	y.ExpBySepticFp(*x)
	for i, r := range septicRoots {
		if y == r {
			return i
		}
	}
	panic("invalid septic character")
}
