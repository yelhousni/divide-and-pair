package fp2

import (
	"github.com/yelhousni/divide-and-pair/internal/octic"
)

// pi=1+2^63*(t-t^3) has absolute norm (1-2^127)^2=p^2.
// In the residue field t^2=i and sqrt(2)=t-t^3=-2^64, so
// t=-2^63*(1+i). The quotient is Fp[i] and pi is primary.
var octicPrime = octic.Vector{{Lo: 1}, {Lo: 1 << 63}, {}, {Lo: 1 << 63, Hi: -1}}

func octicLift(x *E2) octic.Vector {
	// The raw limbs share a nonzero Fp Montgomery factor. Every such factor
	// is invisible to chi_8 since 8 divides p+1. Replace each high limb's
	// weight 2^64 by -sqrt(2), also using i*sqrt(2)=t+t^3.
	a, b := x.A0, x.A1
	ah, bh := octic.Word{Lo: a[1]}, octic.Word{Lo: b[1]}
	return octic.Vector{{Lo: a[0]}, (octic.Word{}).Sub(ah).Sub(bh), {Lo: b[0]}, ah.Sub(bh)}
}

var octicRoots = func() (roots [8]E2) {
	roots[0].SetOne()
	roots[1].A0.SetUint64(1 << 63)
	roots[1].A0.Neg(&roots[1].A0)
	roots[1].A1 = roots[1].A0
	for i := 2; i < 8; i++ {
		roots[i].Mul(&roots[i-1], &roots[1])
	}
	return
}()

// OcticSymbol returns j with x^((p^2-1)/8)=t^j, t=-2^63*(1+i),
// or -1 for zero. Joye's Euclidean descent is variable-time.
func (x *E2) OcticSymbol() int {
	if x.IsZero() {
		return -1
	}
	if j, ok := octic.Descent(octicLift(x), octicPrime); ok {
		return j
	}
	return x.octicSymbolExp()
}

// For nonzero x, y=x^((p+1)/8) satisfies conjugate(y)=chi_8(x)*y.
// Recover the full index without an inversion.
func (x *E2) octicSymbolExp() int {
	y := *x
	for range 124 {
		y.Square(&y)
	}
	var conjugate, product E2
	conjugate.Conjugate(&y)
	for i := range octicRoots {
		product.Mul(&y, &octicRoots[i])
		if product == conjugate {
			return i
		}
	}
	panic("invalid octic character")
}
