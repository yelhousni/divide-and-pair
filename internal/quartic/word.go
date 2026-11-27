package quartic

import "math/bits"

// chi(i) and chi(1+i) as exponents of i for a primary denominator c+di.
// The numerators in the supplementary laws are respectively even and
// divisible by four, so their low bits suffice.
func Supp(c, d uint64) (unit, ramified uint64) {
	return (1 - c) >> 1, (c - d - d*d - 1) >> 2
}

// Once both 1-norms are < 2^61, the same reduction fits in signed words,
// including intermediate sums. Complex magnitudes decrease, so all
// subsequent 1-norms are at most sqrt(2) times the starting maximum.
func Word(fr, fi, gr, gi int64, e uint64) uint8 {
	for {
		if (fi == 0 && fr == 1) || (gi == 0 && gr == 1) {
			return uint8(e & 3)
		}
		if absWord(fr)+absWord(fi) > absWord(gr)+absWord(gi) {
			e += uint64(fr & gr & 2)
			fr, fi, gr, gi = gr, gi, fr, fi
		}
		gr -= fr
		gi -= fi
		k := bits.TrailingZeros64(uint64(gr | gi))
		gr, gi = gr>>k, gi>>k
		unit, ramified := Supp(uint64(fr), uint64(fi))
		e += uint64(k) * (2*ramified - unit)
		if gr&gi&1 != 0 {
			gr, gi = (gr+gi)>>1, (gi-gr)>>1
			e += ramified
		}
		if gi&1 != 0 {
			gr, gi = -gi, gr
			e -= unit
		}
		if (gr+gi)&2 != 0 {
			gr, gi = -gr, -gi
			e -= 2 * unit
		}
	}
}

func absWord(v int64) int64 {
	if v < 0 {
		return -v
	}
	return v
}
