package fp2

import (
	"math/big"
	"sync"
)

// QuarticSymbolGCD computes χ₄(z) = z^((p²-1)/4) for z ∈ Fp2 using a
// Euclidean GCD algorithm over Z[i].
//
// Since p ≡ 3 mod 4, p is a Gaussian prime in Z[i] (it remains inert).
// The quartic residue symbol [z/p]_4 is computed by reducing z modulo p
// in Z[i] via iterated Gaussian Euclidean division, tracking the symbol
// via quartic reciprocity and supplementary laws.
//
// Returns 0 (χ₄=1), 1 (χ₄=i), 2 (χ₄=-1), or 3 (χ₄=-i).
func (z *E2) QuarticSymbolGCD() uint8 {
	if z.IsZero() {
		return 0
	}

	// Convert z to big.Int Gaussian integer (a + b*i)
	sc := gcdPool.Get().(*gcdScratch)
	defer gcdPool.Put(sc)

	z.A0.BigInt(&sc.aRe)
	z.A1.BigInt(&sc.aIm)

	// b = (-p, 0): the primary Gaussian prime (p ≡ 3 mod 4, so -p ≡ 1 mod 4 is primary)
	sc.bRe.Neg(&biP)
	sc.bIm.SetUint64(0)

	result := int64(0)

	// p is primary: p ≡ 3 mod 4, so p + 0i has p ≡ 3 mod (1+i)^3.
	// For the quartic symbol [a/b]_4, both a and b need to be primary.
	// We make a primary first: extract (1+i) factors and i factors.

	// Extract (1+i) factors from a: while a.re ≡ a.im mod 2, divide by (1+i)
	m := int64(0)
	for biParity(&sc.aRe) == biParity(&sc.aIm) {
		biDivBy1PlusI(&sc.aRe, &sc.aIm, &sc.t1)
		m++
		if sc.aRe.Sign() == 0 && sc.aIm.Sign() == 0 {
			return 0 // z divides a power of (1+i), hence divides p
		}
	}

	// Make a primary: a ≡ 1 mod (1+i)^3, i.e., a.im even and a.re+a.im ≡ 1 mod 4
	np := biMakePrimary(&sc.aRe, &sc.aIm, &sc.t1)

	// Initial symbol contribution from (1+i)^m and i^np extracted from a
	bl := biMod64(&sc.bRe)
	dl := int64(0) // bIm = 0
	el := biMod8(&sc.aRe)
	result = (m*biSup1PlusI(bl, dl) + int64(np)*biSupI(bl&7)) & 3

	// Euclidean GCD loop
	for iter := 0; iter < 2000; iter++ {
		if sc.bRe.Sign() == 0 && sc.bIm.Sign() == 0 {
			break
		}
		// Check if b is a unit
		isUnit := false
		if sc.bIm.Sign() == 0 {
			if sc.bRe.CmpAbs(biOne) <= 0 {
				isUnit = true
			}
		} else if sc.bRe.Sign() == 0 {
			if sc.bIm.CmpAbs(biOne) <= 0 {
				isUnit = true
			}
		}
		if isUnit {
			k := biUnitPow(&sc.bRe, &sc.bIm)
			result = (result + int64(k)*biSupI(biMod8(&sc.aRe))) & 3
			break
		}

		// Gaussian remainder: e = a mod b
		biGaussRem(&sc.eRe, &sc.eIm, &sc.aRe, &sc.aIm, &sc.bRe, &sc.bIm, &sc.t1, &sc.t2, &sc.t3, &sc.t4)

		if sc.eRe.Sign() == 0 && sc.eIm.Sign() == 0 {
			return 0 // gcd found, not coprime
		}

		// Extract (1+i) factors
		m = 0
		for biParity(&sc.eRe) == biParity(&sc.eIm) {
			biDivBy1PlusI(&sc.eRe, &sc.eIm, &sc.t1)
			m++
		}

		// Make primary
		np = biMakePrimary(&sc.eRe, &sc.eIm, &sc.t1)

		// Update symbol
		bl = biMod64(&sc.bRe)
		dl = biMod64(&sc.bIm)
		el = biMod8(&sc.eRe)
		s1pi := biSup1PlusI(bl, dl)
		sI := biSupI(bl & 7)
		rec := biReciprocity(el, bl&7)
		result = (result + m*s1pi + int64(np)*sI + 2*rec) & 3

		// Swap: a, b = b, e
		sc.aRe.Set(&sc.bRe)
		sc.aIm.Set(&sc.bIm)
		sc.bRe.Set(&sc.eRe)
		sc.bIm.Set(&sc.eIm)
	}

	return uint8(result & 3)
}

// gcdScratch holds pre-allocated big.Int temporaries.
type gcdScratch struct {
	aRe, aIm, bRe, bIm, eRe, eIm big.Int
	t1, t2, t3, t4                 big.Int
}

var gcdPool = sync.Pool{
	New: func() any {
		return new(gcdScratch)
	},
}

var (
	biP   big.Int
	biOne = big.NewInt(1)
)

func init() {
	// p = 2^448 - 2^224 - 1
	biP.SetBit(&biP, 448, 1)
	var t big.Int
	t.SetBit(&t, 224, 1)
	biP.Sub(&biP, &t)
	biP.Sub(&biP, big.NewInt(1))
}

// biParity returns z mod 2.
func biParity(z *big.Int) uint {
	w := z.Bits()
	if len(w) == 0 {
		return 0
	}
	return uint(w[0]) & 1
}

// biMod4 returns z mod 4 in [0,3].
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

// biMod8 returns z mod 8 in [0,7].
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

// biMod64 returns z mod 64 in [0,63].
func biMod64(z *big.Int) int64 {
	w := z.Bits()
	if len(w) == 0 {
		return 0
	}
	v := int64(w[0] & 63)
	if z.Sign() < 0 && v != 0 {
		v = 64 - v
	}
	return v
}

// biDivBy1PlusI divides a+bi by (1+i): (a+b)/2 + (b-a)/2 * i
func biDivBy1PlusI(a, b, t *big.Int) {
	t.Add(a, b)
	b.Sub(b, a)
	a.Rsh(t, 1)
	b.Rsh(b, 1)
}

// biMakePrimary adjusts (a, b) to be primary (b even, a+b ≡ 1 mod 4)
// by multiplying by a power of i. Returns the power n ∈ [0,3].
func biMakePrimary(a, b, t *big.Int) int64 {
	for n := int64(0); n < 4; n++ {
		if biParity(b) == 0 {
			t.Add(a, b)
			if biMod4(t) == 1 {
				return n
			}
		}
		// Multiply by i: (a,b) -> (-b, a)
		t.Set(a)
		a.Neg(b)
		b.Set(t)
	}
	return 0
}

// biGaussRem computes e = a - round(a/b)*b in Z[i].
func biGaussRem(eRe, eIm, aRe, aIm, bRe, bIm, t1, t2, t3, t4 *big.Int) {
	// q = round((a * conj(b)) / N(b))
	// a * conj(b) = (aRe*bRe + aIm*bIm) + (aIm*bRe - aRe*bIm)*i
	// N(b) = bRe^2 + bIm^2
	t1.Mul(aRe, bRe)
	t2.Mul(aIm, bIm)
	t3.Add(t1, t2) // numRe = aRe*bRe + aIm*bIm

	t1.Mul(aIm, bRe)
	t2.Mul(aRe, bIm)
	t4.Sub(t1, t2) // numIm = aIm*bRe - aRe*bIm

	t1.Mul(bRe, bRe)
	t2.Mul(bIm, bIm)
	t1.Add(t1, t2) // norm = bRe^2 + bIm^2

	// qRe = round(numRe / norm), qIm = round(numIm / norm)
	biRoundDiv(eRe, t3, t1, t2) // eRe = round(numRe/norm) temporarily
	biRoundDiv(eIm, t4, t1, t2) // eIm = round(numIm/norm) temporarily

	// e = a - q*b
	// q*b = (qRe*bRe - qIm*bIm) + (qRe*bIm + qIm*bRe)*i
	qRe := new(big.Int).Set(eRe)
	qIm := new(big.Int).Set(eIm)

	t1.Mul(qRe, bRe)
	t2.Mul(qIm, bIm)
	t3.Sub(t1, t2) // qb_re

	t1.Mul(qRe, bIm)
	t2.Mul(qIm, bRe)
	t4.Add(t1, t2) // qb_im

	eRe.Sub(aRe, t3)
	eIm.Sub(aIm, t4)
}

// biRoundDiv computes z = round(a/b). Uses t as temporary.
func biRoundDiv(z, a, b, t *big.Int) {
	z.QuoRem(a, b, t)
	t.Abs(t)
	t.Lsh(t, 1)
	if t.Cmp(b) > 0 {
		if a.Sign() >= 0 {
			z.Add(z, biOne)
		} else {
			z.Sub(z, biOne)
		}
	}
}

// biUnitPow returns n such that z = i^n for a Gaussian unit z.
func biUnitPow(re, im *big.Int) int64 {
	if im.Sign() == 0 {
		if re.Sign() > 0 {
			return 0
		}
		return 2
	}
	if re.Sign() == 0 {
		if im.Sign() > 0 {
			return 1
		}
		return 3
	}
	return 0
}

func biSupI(c int64) int64 {
	c = ((c % 8) + 8) % 8
	return ((1 - c) / 2) & 3
}

func biSup1PlusI(c, d int64) int64 {
	c = ((c % 64) + 64) % 64
	d = ((d % 64) + 64) % 64
	return ((c - d - d*d - 1) / 4) & 3
}

func biReciprocity(a, c int64) int64 {
	a = ((a % 8) + 8) % 8
	c = ((c % 8) + 8) % 8
	return (((a - 1) / 2) & 1) * (((c - 1) / 2) & 1)
}
