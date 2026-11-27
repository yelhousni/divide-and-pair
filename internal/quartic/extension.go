package quartic

import "math/bits"

// Seven data limbs cover Curve448. The eighth is a sign and carry guard.
// Only the active low n limbs participate as the operands shrink.
type coord [8]uint64
type gaussian struct{ r, i coord }

func neg(z *coord, n int) {
	var borrow uint64
	for j := 0; j < n; j++ {
		z[j], borrow = bits.Sub64(0, z[j], borrow)
	}
}

func subtract(z *coord, x *coord, n int) {
	var borrow uint64
	for j := 0; j < n; j++ {
		z[j], borrow = bits.Sub64(z[j], x[j], borrow)
	}
}

func shift(z *coord, k, n int) {
	words, s := k/64, uint(k%64)
	sign := uint64(int64(z[n-1]) >> 63)
	for j := 0; j < n-words; j++ {
		next := sign
		if j+words+1 < n {
			next = z[j+words+1]
		}
		z[j] = z[j+words]>>s | next<<(64-s)
	}
	for j := n - words; j < n; j++ {
		z[j] = sign
	}
}

func (g *gaussian) primary(f *gaussian, n int) uint64 {
	k := 0
	for j := 0; j < n; j++ {
		count := bits.TrailingZeros64(g.r[j] | g.i[j])
		k += count
		if count < 64 {
			break
		}
	}
	if k != 0 {
		shift(&g.r, k, n)
		shift(&g.i, k, n)
	}
	unit, ramified := Supp(f.r[0], f.i[0])
	e := uint64(k) * (2*ramified - unit)
	if g.r[0]&g.i[0]&1 != 0 {
		var c, b uint64
		for j := 0; j < n; j++ {
			r, i := g.r[j], g.i[j]
			g.r[j], c = bits.Add64(r, i, c)
			g.i[j], b = bits.Sub64(i, r, b)
		}
		shift(&g.r, 1, n)
		shift(&g.i, 1, n)
		e += ramified
	}
	if g.i[0]&1 != 0 {
		g.r, g.i = g.i, g.r
		neg(&g.r, n)
		e -= unit
	}
	if (g.r[0]+g.i[0])&2 != 0 {
		neg(&g.r, n)
		neg(&g.i, n)
		e -= 2 * unit
	}
	return e
}

func (g *gaussian) size(n int) (s coord) {
	// Absolute values via conditional two's-complement negation.
	rs, is := uint64(int64(g.r[n-1])>>63), uint64(int64(g.i[n-1])>>63)
	rc, ic := rs&1, is&1
	var carry uint64
	for j := 0; j < n; j++ {
		r, c := bits.Add64(g.r[j]^rs, 0, rc)
		rc = c
		i, c := bits.Add64(g.i[j]^is, 0, ic)
		ic = c
		s[j], carry = bits.Add64(r, i, carry)
	}
	return
}

func greater(a, b *coord, n int) bool {
	for j := n - 1; j >= 0; j-- {
		if a[j] != b[j] {
			return a[j] > b[j]
		}
	}
	return false
}

// Extension returns the quartic character index relative to i in Fp[i],
// or -1 for zero. The input prime is congruent to 3 modulo 4. Components and p are little endian, on at most
// seven limbs. A common nonzero base-field scale has trivial quartic character.
func Extension(re, im, p []uint64) int {
	n := len(p) + 1
	var a, b gaussian
	copy(a.r[:], p)
	neg(&a.r, n) // -p is primary and generates (p) in Z[i].
	copy(b.r[:], re)
	copy(b.i[:], im)
	if b == (gaussian{}) {
		return -1
	}
	f, g := &a, &b
	e := g.primary(f, n)
	sf, sg := f.size(n), g.size(n)
	one := coord{1}
	for {
		// Leave three guard bits for the difference, the 1+i step and
		// any growth in the 1-norm as the complex magnitude decreases.
		for n > 1 && sf[n-1]|sg[n-1] == 0 && sf[n-2]|sg[n-2] < 1<<61 {
			n--
		}
		if n == 1 {
			return int(Word(int64(f.r[0]), int64(f.i[0]), int64(g.r[0]), int64(g.i[0]), e))
		}
		if sf == one || sg == one {
			return int(e & 3)
		}
		// For the operand with larger 1-norm, |f| <= sqrt(2)*|g|.
		// Both are primary, so (1+i)^3 divides g-f. Removing it gives
		// a new complex magnitude <= (1+sqrt(2))*|g|/sqrt(8) < |g|.
		// Coprimality rules out equal nonunit operands.
		if greater(&sf, &sg, n) {
			e += f.r[0] & g.r[0] & 2
			f, g = g, f
			sf = sg
		}
		subtract(&g.r, &f.r, n)
		subtract(&g.i, &f.i, n)
		e += g.primary(f, n)
		sg = g.size(n)
	}
}
