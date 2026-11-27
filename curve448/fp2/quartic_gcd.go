package fp2

import (
	"github.com/yelhousni/divide-and-pair/curve448/fp"
	"github.com/yelhousni/divide-and-pair/internal/quartic"
)

var quarticModulus = func() [fp.Limbs]uint64 {
	var minusOne fp.Element
	minusOne.SetOne()
	minusOne.Neg(&minusOne)
	p := minusOne.Bits()
	p[0]++ // p is odd, so p-1 has an even low word.
	return p
}()

// QuarticSymbol returns j such that x^((p²-1)/4)=i^j, or -1 for zero.
// The Gaussian binary GCD is variable-time. Its input may retain the common
// Montgomery factor because every element of Fp* is a fourth power in Fp².
func (x *E2) QuarticSymbol() int {
	return quartic.Extension(x.A0[:], x.A1[:], quarticModulus[:])
}
