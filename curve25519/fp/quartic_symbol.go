package fp

import (
	"math/big"
	"math/bits"
	"sync"
)

const (
	piRe0 uint64 = 0x3feab578735893c3
	piRe1 uint64 = 0x33a5cbdded73544f
	piIm0 uint64 = 0x43c900683eb6254a
	piIm1 uint64 = 0xad7eb9766c0b7b36
)

type signed128 struct {
	lo, hi uint64
	neg    bool
}

func (s signed128) isZero() bool   { return s.lo == 0 && s.hi == 0 }
func (s signed128) parity() uint64 { return s.lo & 1 }

func (s signed128) low8() int64 {
	v := int64(s.lo & 7)
	if s.neg && v != 0 {
		v = 8 - v
	}
	return v
}

func (s signed128) low64() int64 {
	v := int64(s.lo & 63)
	if s.neg && v != 0 {
		v = 64 - v
	}
	return v
}

func neg128(s signed128) signed128 {
	if s.isZero() {
		return s
	}
	return signed128{s.lo, s.hi, !s.neg}
}

func add128(a, b signed128) signed128 {
	if a.neg == b.neg {
		lo, c := bits.Add64(a.lo, b.lo, 0)
		hi, _ := bits.Add64(a.hi, b.hi, c)
		return signed128{lo, hi, a.neg}
	}
	if a.hi > b.hi || (a.hi == b.hi && a.lo >= b.lo) {
		lo, bw := bits.Sub64(a.lo, b.lo, 0)
		hi, _ := bits.Sub64(a.hi, b.hi, bw)
		return signed128{lo, hi, a.neg}
	}
	lo, bw := bits.Sub64(b.lo, a.lo, 0)
	hi, _ := bits.Sub64(b.hi, a.hi, bw)
	return signed128{lo, hi, b.neg}
}

func sub128(a, b signed128) signed128 { return add128(a, neg128(b)) }

func rsh128(s signed128) signed128 {
	return signed128{(s.lo >> 1) | (s.hi << 63), s.hi >> 1, s.neg}
}

func mulSmall128(s signed128, k int64) signed128 {
	neg := s.neg
	if k < 0 {
		neg = !neg
		k = -k
	}
	uk := uint64(k)
	hi, lo := bits.Mul64(s.lo, uk)
	hi += s.hi * uk
	return signed128{lo, hi, neg}
}

func s128ToFloat(s signed128) float64 {
	f := float64(s.hi)*0x1p64 + float64(s.lo)
	if s.neg {
		f = -f
	}
	return f
}

func gaussQuotientSmall(aRe, aIm, bRe, bIm signed128) (int64, int64) {
	ar := s128ToFloat(aRe)
	ai := s128ToFloat(aIm)
	br := s128ToFloat(bRe)
	bi := s128ToFloat(bIm)
	nB := br*br + bi*bi
	if nB == 0 {
		return 0, 0
	}
	numRe := ar*br + ai*bi
	numIm := ai*br - ar*bi
	return roundFloat(numRe / nB), roundFloat(numIm / nB)
}

func roundFloat(x float64) int64 {
	if x >= 0 {
		return int64(x + 0.5)
	}
	return -int64(-x + 0.5)
}

func gaussRem128(aRe, aIm, bRe, bIm signed128) (signed128, signed128) {
	qr, qi := gaussQuotientSmall(aRe, aIm, bRe, bIm)
	qbRe := sub128(mulSmall128(bRe, qr), mulSmall128(bIm, qi))
	qbIm := add128(mulSmall128(bIm, qr), mulSmall128(bRe, qi))
	return sub128(aRe, qbRe), sub128(aIm, qbIm)
}

func divBy1PlusI128(e, f *signed128) {
	sum := add128(*e, *f)
	diff := sub128(*f, *e)
	*e = rsh128(sum)
	*f = rsh128(diff)
}

func makePrimary128(e, f *signed128) int64 {
	for n := int64(0); n < 4; n++ {
		if f.parity() == 0 {
			sum := add128(*e, *f)
			m4 := int64(sum.lo & 3)
			if sum.neg && m4 != 0 {
				m4 = 4 - m4
			}
			if m4 == 1 {
				return n
			}
		}
		*e, *f = signed128{f.lo, f.hi, !f.neg}, *e
	}
	return 0
}

// Pre-allocated big.Int constants for Phase 1 (avoid allocations per call).
var (
	biPiRe    big.Int
	biPiIm    big.Int // negative (conj)
	biPiImAbs big.Int // |piIm| (positive) — avoids per-call Neg
	biNormB   big.Int // = q
	biOne     = big.NewInt(1)

	// sqrtMinusOneFp is sqrt(-1) mod q, used by QuarticSymbolExp.
	sqrtMinusOneFp Element
	minusOneFp     Element
)

func init() {
	biPiRe.SetUint64(piRe1)
	biPiRe.Lsh(&biPiRe, 64)
	biPiRe.Or(&biPiRe, new(big.Int).SetUint64(piRe0))
	biPiIm.SetUint64(piIm1)
	biPiIm.Lsh(&biPiIm, 64)
	biPiIm.Or(&biPiIm, new(big.Int).SetUint64(piIm0))
	biPiImAbs.Set(&biPiIm)
	biPiIm.Neg(&biPiIm)
	var t1, t2 big.Int
	t1.Mul(&biPiRe, &biPiRe)
	t2.Mul(&biPiIm, &biPiIm)
	biNormB.Add(&t1, &t2)

	sqrtMinusOneFp.SetString("19681161376707505956807079304988542015446066515923890162744021073123829784752")
	minusOneFp.SetOne()
	minusOneFp.Neg(&minusOneFp)
}

// phase1Scratch holds pre-allocated big.Int temporaries for Phase 1.
// Buffers are pre-sized to avoid allocation during Mul.
type phase1Scratch struct {
	zBI, numRe, numIm, qRe, qIm, t1, t2, e, f, s, d big.Int
}

var phase1Pool = sync.Pool{
	New: func() any {
		sc := new(phase1Scratch)
		// Pre-size buffers to 6 words (384 bits) to avoid growth during Mul.
		for _, p := range []*big.Int{&sc.zBI, &sc.numRe, &sc.numIm, &sc.qRe, &sc.qIm, &sc.t1, &sc.t2, &sc.e, &sc.f, &sc.s, &sc.d} {
			p.SetBits(make([]big.Word, 6))
			p.SetUint64(0)
		}
		return sc
	},
}

// quarticSymbolWeilert implements the Weilert GCD algorithm for the quartic symbol.
// The fallback parameter is called if the GCD exceeds 200 iterations.
func quarticSymbolWeilert(z *Element, fallback func(*Element) uint8) uint8 {
	if z.IsZero() {
		return 0
	}

	bRe := signed128{piRe0, piRe1, false}
	bIm := signed128{piIm0, piIm1, true}
	result := int64(0)
	var aRe, aIm signed128

	// Phase 1: one big Euclidean step (pooled big.Int to avoid allocation).
	sc := phase1Pool.Get().(*phase1Scratch)
	z.BigInt(&sc.zBI)

	sc.numRe.Mul(&sc.zBI, &biPiRe)
	sc.numIm.Mul(&sc.zBI, &biPiImAbs)
	roundDivBI(&sc.qRe, &sc.numRe, &biNormB, &sc.t1, &sc.t2)
	roundDivBI(&sc.qIm, &sc.numIm, &biNormB, &sc.t1, &sc.t2)

	sc.t1.Mul(&sc.qRe, &biPiRe)
	sc.t2.Mul(&sc.qIm, &biPiIm)
	sc.e.Sub(&sc.t1, &sc.t2)
	sc.e.Sub(&sc.zBI, &sc.e)
	sc.t1.Mul(&sc.qRe, &biPiIm)
	sc.t2.Mul(&sc.qIm, &biPiRe)
	sc.f.Add(&sc.t1, &sc.t2)
	sc.f.Neg(&sc.f)

	if sc.e.Sign() == 0 && sc.f.Sign() == 0 {
		return 0
	}

	// Remove (1+i) factors using only parity checks (no big.Int.Mod).
	m := int64(0)
	for biParity(&sc.e) == biParity(&sc.f) {
		sc.s.Add(&sc.e, &sc.f)
		sc.d.Sub(&sc.f, &sc.e)
		sc.e.Rsh(&sc.s, 1)
		sc.f.Rsh(&sc.d, 1)
		m++
	}

	// Make primary: f even, e+f ≡ 1 mod 4
	nc := int64(0)
	for n := int64(0); n < 4; n++ {
		if biParity(&sc.f) == 0 {
			sc.s.Add(&sc.e, &sc.f)
			if biMod4(&sc.s) == 1 {
				nc = n
				break
			}
		}
		sc.s.Set(&sc.e)
		sc.e.Neg(&sc.f)
		sc.f.Set(&sc.s)
		nc = n + 1
	}
	np := (4 - nc) % 4

	bl := int64(piRe0 & 63)
	dl := (-(int64(piIm0 & 63))) & 63
	el := biMod8(&sc.e)
	result = (m*sup1PlusILow(bl, dl) + int64(np)*supILow(bl&7) + 2*reciprocityLow(el, bl&7)) & 3

	bigToS128(&aRe, &sc.e)
	bigToS128(&aIm, &sc.f)
	phase1Pool.Put(sc) // return to pool
	aRe, bRe = bRe, aRe
	aIm, bIm = bIm, aIm

	// Phase 2: pure uint64 Euclidean GCD.
	for iter := 0; ; iter++ {
		if iter > 200 {
			return fallback(z)
		}
		if bIm.isZero() && bRe.hi == 0 && bRe.lo <= 1 {
			break
		}
		if bRe.isZero() && bIm.hi == 0 && bIm.lo <= 1 {
			break
		}
		if aRe.isZero() && aIm.isZero() {
			return 0
		}
		isAUnit := (aIm.isZero() && aRe.hi == 0 && aRe.lo <= 1) ||
			(aRe.isZero() && aIm.hi == 0 && aIm.lo <= 1)
		if isAUnit {
			k := unitPow128(aRe, aIm)
			result = (result + int64(k)*supILow(bRe.low8())) & 3
			break
		}
		eRe, eIm := gaussRem128(aRe, aIm, bRe, bIm)
		if eRe.isZero() && eIm.isZero() {
			return 0
		}
		m := int64(0)
		for eRe.parity() == eIm.parity() {
			divBy1PlusI128(&eRe, &eIm)
			m++
		}
		nc := makePrimary128(&eRe, &eIm)
		np := (4 - nc) % 4
		s1pi := sup1PlusILow(bRe.low64(), bIm.low64())
		sI := supILow(bRe.low8())
		rec := reciprocityLow(eRe.low8(), bRe.low8())
		result = (result + m*s1pi + int64(np)*sI + 2*rec) & 3
		aRe, aIm = bRe, bIm
		bRe, bIm = eRe, eIm
	}

	return uint8(result & 3)
}

// QuarticSymbol computes χ₄(z) = z^((q-1)/4) mod q using Weilert's algorithm.
// Returns 0 (χ₄=1), 1 (χ₄=i), 2 (χ₄=-1), or 3 (χ₄=-i).
func (z *Element) QuarticSymbol() uint8 {
	return quarticSymbolWeilert(z, (*Element).QuarticSymbolExp)
}

func unitPow128(re, im signed128) int64 {
	if im.isZero() {
		if !re.neg {
			return 0
		}
		return 2
	}
	if re.isZero() {
		if !im.neg {
			return 1
		}
		return 3
	}
	return 0
}

func supILow(c int64) int64 { c = ((c % 8) + 8) % 8; return ((1 - c) / 2) & 3 }
func sup1PlusILow(c, d int64) int64 {
	c = ((c % 64) + 64) % 64
	d = ((d % 64) + 64) % 64
	return ((c - d - d*d - 1) / 4) & 3
}
func reciprocityLow(a, c int64) int64 {
	a = ((a % 8) + 8) % 8
	c = ((c % 8) + 8) % 8
	return (((a - 1) / 2) & 1) * (((c - 1) / 2) & 1)
}

// --- Zero-alloc big.Int helpers ---

// biParity returns z mod 2 (0 or 1), without allocation.
func biParity(z *big.Int) uint {
	w := z.Bits()
	if len(w) == 0 {
		return 0
	}
	return uint(w[0]) & 1
}

// biMod4 returns z mod 4 in [0,3], without allocation.
func biMod4(z *big.Int) int64 {
	w := z.Bits()
	if len(w) == 0 {
		return 0
	}
	v := int64(w[0] & 3)
	if z.Sign() < 0 && v != 0 {
		v = 4 - v
	}
	return v
}

// biMod8 returns z mod 8 in [0,7], without allocation.
func biMod8(z *big.Int) int64 {
	w := z.Bits()
	if len(w) == 0 {
		return 0
	}
	v := int64(w[0] & 7)
	if z.Sign() < 0 && v != 0 {
		v = 8 - v
	}
	return v
}

// roundDivBI computes z = round(a/b) for b > 0, using pre-allocated temps.
func roundDivBI(z, a, b *big.Int, q, r *big.Int) {
	q.QuoRem(a, b, r)
	r.Abs(r)
	r.Lsh(r, 1)
	if r.Cmp(b) > 0 {
		if a.Sign() >= 0 {
			z.Add(q, biOne)
		} else {
			z.Sub(q, biOne)
		}
	} else {
		z.Set(q)
	}
}

func bigToS128(s *signed128, x *big.Int) {
	s.neg = x.Sign() < 0
	w := x.Bits() // no allocation — returns internal slice
	s.lo, s.hi = 0, 0
	if len(w) > 0 {
		s.lo = uint64(w[0])
	}
	if len(w) > 1 {
		s.hi = uint64(w[1])
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
