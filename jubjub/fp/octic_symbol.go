package fp

import (
	"math/bits"

	"github.com/yelhousni/divide-and-pair/internal/octic"
)

// Prime above p with t mapping to 7^((p-1)/8).
var octicMod = octic.Vector{
	{Lo: 6451692399033376397}, {Lo: 11914338290826914742},
	{Lo: 1094064715653685168}, {Lo: ^uint64(11306383988181076382 - 1), Hi: -1},
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
func octicLiftQuotient(x *Element, c [3]uint64) octic.Word {
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
	return octic.Word{Lo: lo, Hi: int64(p[5] + carry)}
}

func octicLift(x *Element) octic.Vector {
	var q octic.Vector
	for i, c := range octicRecip {
		q[i] = octicLiftQuotient(x, c)
		if i != 1 {
			q[i] = (octic.Word{}).Sub(q[i])
		}
	}
	// A raw Montgomery representative differs from x by the eighth power
	// 2^256. Each quotient coefficient has error <= 1/2 + p/2^257 < 3/4.
	// Thus every remainder coefficient is < 3/4 * sum |pi_i| < 2^65.
	b := octic.Vector{{Lo: x[0], Hi: int64(x[1])}}
	return b.Remainder(octicMod, q)
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
	if j, ok := octic.Descent(octicLift(x), octicMod); ok {
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
