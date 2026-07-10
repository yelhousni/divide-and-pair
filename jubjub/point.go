package jubjub

import (
	"math/big"
	"math/bits"

	fp "github.com/yelhousni/divide-and-pair/jubjub/fp"
)

// PointAffine point on a twisted Edwards curve in affine coordinates
type PointAffine struct {
	X, Y fp.Element
}

// PointProj point in projective coordinates (X:Y:Z) with x=X/Z, y=Y/Z
type PointProj struct {
	X, Y, Z fp.Element
}

// Set sets p to p1 and returns p
func (p *PointAffine) Set(p1 *PointAffine) *PointAffine {
	p.X.Set(&p1.X)
	p.Y.Set(&p1.Y)
	return p
}

// Equal returns true if p and p1 are equal
func (p *PointAffine) Equal(p1 *PointAffine) bool {
	return p.X.Equal(&p1.X) && p.Y.Equal(&p1.Y)
}

// IsZero returns true if p is the identity (0, 1)
func (p *PointAffine) IsZero() bool {
	return p.X.IsZero() && p.Y.IsOne()
}

// IsOnCurve checks if the point is on the curve: Ax^2 + y^2 = 1 + Dx^2*y^2
func (p *PointAffine) IsOnCurve() bool {
	initOnce.Do(initCurveParams)

	var lhs, rhs, tmp fp.Element

	tmp.Mul(&p.Y, &p.Y)
	lhs.Mul(&p.X, &p.X)
	mulByA(&lhs)
	lhs.Add(&lhs, &tmp) // A*x^2 + y^2

	tmp.Mul(&p.X, &p.X).
		Mul(&tmp, &p.Y).
		Mul(&tmp, &p.Y).
		Mul(&tmp, &curveParams.D)
	rhs.SetOne().Add(&rhs, &tmp) // 1 + D*x^2*y^2

	return lhs.Equal(&rhs)
}

// Neg negates a point
func (p *PointAffine) Neg(p1 *PointAffine) *PointAffine {
	p.X.Neg(&p1.X)
	p.Y.Set(&p1.Y)
	return p
}

// Add adds two affine points using the unified addition formula for twisted Edwards curves
func (p *PointAffine) Add(p1, p2 *PointAffine) *PointAffine {
	initOnce.Do(initCurveParams)

	// x3 = (x1*y2 + y1*x2) / (1 + d*x1*x2*y1*y2)
	// y3 = (y1*y2 - a*x1*x2) / (1 - d*x1*x2*y1*y2)
	var x1y2, y1x2, y1y2, x1x2, dxy fp.Element

	x1y2.Mul(&p1.X, &p2.Y)
	y1x2.Mul(&p1.Y, &p2.X)
	y1y2.Mul(&p1.Y, &p2.Y)
	x1x2.Mul(&p1.X, &p2.X)
	dxy.Mul(&curveParams.D, &x1x2).Mul(&dxy, &y1y2)

	var numX, numY, denX, denY fp.Element
	numX.Add(&x1y2, &y1x2)
	ax1x2 := x1x2
	mulByA(&ax1x2)
	numY.Sub(&y1y2, &ax1x2) // y1y2 - a*x1x2 = y1y2 + x1x2 (since a=-1)

	denX.SetOne().Add(&denX, &dxy)
	denY.SetOne().Sub(&denY, &dxy)

	denX.Inverse(&denX)
	denY.Inverse(&denY)

	p.X.Mul(&numX, &denX)
	p.Y.Mul(&numY, &denY)
	return p
}

// Double doubles a point
func (p *PointAffine) Double(p1 *PointAffine) *PointAffine {
	return p.Add(p1, p1)
}

// ScalarMultiplication scalar multiplication of a point
func (p *PointAffine) ScalarMultiplication(p1 *PointAffine, scalar *big.Int) *PointAffine {
	var proj, res PointProj
	proj.FromAffine(p1)
	res.ScalarMultiplication(&proj, scalar)
	res.ToAffine(p)
	return p
}

// --- Projective coordinates ---

// Set sets p to p1
func (p *PointProj) Set(p1 *PointProj) *PointProj {
	p.X.Set(&p1.X)
	p.Y.Set(&p1.Y)
	p.Z.Set(&p1.Z)
	return p
}

// FromAffine converts affine to projective
func (p *PointProj) FromAffine(a *PointAffine) *PointProj {
	p.X.Set(&a.X)
	p.Y.Set(&a.Y)
	p.Z.SetOne()
	return p
}

// ToAffine converts projective to affine
func (p *PointProj) ToAffine(a *PointAffine) *PointAffine {
	var zInv fp.Element
	zInv.Inverse(&p.Z)
	a.X.Mul(&p.X, &zInv)
	a.Y.Mul(&p.Y, &zInv)
	return a
}

// IsZero returns true if p is the identity (0:1:1)
func (p *PointProj) IsZero() bool {
	return p.X.IsZero() && p.Y.Equal(&p.Z)
}

func (p *PointProj) setInfinity() *PointProj {
	p.X.SetZero()
	p.Y.SetOne()
	p.Z.SetOne()
	return p
}

func (p *PointProj) Neg(p1 *PointProj) *PointProj {
	p.X.Neg(&p1.X)
	p.Y.Set(&p1.Y)
	p.Z.Set(&p1.Z)
	return p
}

func (p *PointProj) Add(p1, p2 *PointProj) *PointProj {
	initOnce.Do(initCurveParams)
	var A, B, C, D, E, F, G, H, I fp.Element
	A.Mul(&p1.Z, &p2.Z)
	B.Square(&A)
	C.Mul(&p1.X, &p2.X)
	D.Mul(&p1.Y, &p2.Y)
	E.Mul(&curveParams.D, &C)
	E.Mul(&E, &D)
	F.Sub(&B, &E)
	G.Add(&B, &E)
	H.Add(&p1.X, &p1.Y)
	I.Add(&p2.X, &p2.Y)
	p.X.Mul(&H, &I).Sub(&p.X, &C).Sub(&p.X, &D).Mul(&p.X, &A).Mul(&p.X, &F)
	var aC fp.Element
	aC.Set(&C)
	mulByA(&aC)
	p.Y.Sub(&D, &aC).Mul(&p.Y, &A).Mul(&p.Y, &G)
	p.Z.Mul(&F, &G)
	return p
}

func (p *PointProj) Double(p1 *PointProj) *PointProj {
	var B, C, D, E, F, H, J fp.Element
	B.Add(&p1.X, &p1.Y).Square(&B)
	C.Square(&p1.X)
	D.Square(&p1.Y)
	E.Set(&C)
	mulByA(&E)
	F.Add(&E, &D)
	H.Square(&p1.Z)
	J.Sub(&F, &H).Sub(&J, &H)
	p.X.Sub(&B, &C).Sub(&p.X, &D).Mul(&p.X, &J)
	p.Y.Sub(&E, &D).Mul(&p.Y, &F)
	p.Z.Mul(&F, &J)
	return p
}

func (p *PointProj) ScalarMultiplication(p1 *PointProj, scalar *big.Int) *PointProj {
	var _scalar big.Int
	_scalar.Set(scalar)
	p.Set(p1)
	if _scalar.Sign() == -1 {
		_scalar.Neg(&_scalar)
		p.Neg(p)
	}
	var res PointProj
	res.setInfinity()
	const wordSize = bits.UintSize
	sWords := _scalar.Bits()
	for i := len(sWords) - 1; i >= 0; i-- {
		ithWord := sWords[i]
		for k := 0; k < wordSize; k++ {
			res.Double(&res)
			kthBit := (ithWord >> (wordSize - 1 - k)) & 1
			if kthBit == 1 {
				res.Add(&res, p)
			}
		}
	}
	p.Set(&res)
	return p
}

// ClearCofactor sets p to [h]p, where h is the curve cofactor, mapping p
// into the prime-order subgroup, and returns p. It uses a fixed short
// addition chain (see mulByCofactor).
func (p *PointAffine) ClearCofactor() *PointAffine {
	var proj PointProj
	proj.FromAffine(p)
	proj.mulByCofactor(&proj)
	proj.ToAffine(p)
	return p
}

// Code generated by internal/generator BEGIN mulByOrder. DO NOT EDIT.

// mulByOrder computes p = [ℓ]a and returns p, where
//
//	ℓ = 6554484396890773809930967563523245729705921265872317281365359162392183254199
//
// is the 252-bit prime subgroup order. The scalar multiplication is unrolled
// as a short addition chain generated with github.com/mmcloughlin/addchain v0.4.0.
//
// The unified Add/Double formulas are valid for arbitrary curve points
// (including the identity and points outside the ℓ-torsion subgroup), so the
// result is the identity if and only if a is in the prime-order subgroup.
//
// Addition chain:
//
//	_10       = 2*1
//	_100      = 2*_10
//	_110      = _10 + _100
//	_111      = 1 + _110
//	_1010     = _100 + _110
//	_1011     = 1 + _1010
//	_1101     = _10 + _1011
//	_10010    = _111 + _1011
//	_10100    = _10 + _10010
//	_101000   = 2*_10100
//	_110011   = _1011 + _101000
//	_110101   = _10 + _110011
//	_111010   = _111 + _110011
//	_111011   = 1 + _111010
//	_1000001  = _110 + _111011
//	_1001011  = _1010 + _1000001
//	_1001101  = _10 + _1001011
//	_1001111  = _10 + _1001101
//	_1010011  = _100 + _1001111
//	_1100111  = _10100 + _1010011
//	_1101111  = _110101 + _111010
//	_10000001 = _10010 + _1101111
//	_10000011 = _10 + _10000001
//	_10010111 = _10100 + _10000011
//	_10011001 = _10 + _10010111
//	_10011101 = _100 + _10011001
//	_11010111 = _111010 + _10011101
//	_11011011 = _100 + _11010111
//	_11100101 = _1010 + _11011011
//	_11100111 = _10 + _11100101
//	i58       = ((_11100111 << 8 + _11011011) << 9 + _10011101) << 9
//	i78       = ((_10011001 + i58) << 9 + _10011001) << 8 + _11010111
//	i105      = ((i78 << 6 + _110101) << 10 + _10000011) << 9
//	i124      = ((_1100111 + i105) << 8 + _111011) << 8 + 1
//	i165      = ((i124 << 14 + _1001101) << 10 + _111011) << 15
//	i186      = ((_1010011 + i165) << 6 + _1101) << 12 + _1000001
//	i215      = ((i186 << 9 + _1001111) << 8 + _110011) << 10
//	i234      = ((_10000001 + i215) << 11 + _1000001) << 5 + _1101
//	i268      = ((i234 << 12 + _10010111) << 12 + _11100101) << 8
//	i288      = ((_11100111 + i268 + _110) << 8 + _1101111) << 9
//	return      ((_11100101 + i288) << 7 + _1001011) << 4 + _111
//
// Operations: 247 point doublings, 55 point additions.
func (p *PointProj) mulByOrder(a *PointProj) *PointProj {
	// x is a working copy of the input so that p may alias a.
	var x, z PointProj
	var t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21 PointProj
	x.Set(a)

	// Step 1: t4 = [0x2]a
	t4.Double(&x)
	// Step 2: t20 = [0x4]a
	t20.Double(&t4)
	// Step 3: t3 = [0x6]a
	t3.Add(&t4, &t20)
	// Step 4: z = [0x7]a
	z.Add(&x, &t3)
	// Step 5: t1 = [0xa]a
	t1.Add(&t20, &t3)
	// Step 6: t0 = [0xb]a
	t0.Add(&x, &t1)
	// Step 7: t6 = [0xd]a
	t6.Add(&t4, &t0)
	// Step 8: t8 = [0x12]a
	t8.Add(&z, &t0)
	// Step 9: t5 = [0x14]a
	t5.Add(&t4, &t8)
	// Step 10: t2 = [0x28]a
	t2.Double(&t5)
	// Step 11: t9 = [0x33]a
	t9.Add(&t0, &t2)
	// Step 12: t16 = [0x35]a
	t16.Add(&t4, &t9)
	// Step 13: t17 = [0x3a]a
	t17.Add(&z, &t9)
	// Step 14: t12 = [0x3b]a
	t12.Add(&x, &t17)
	// Step 15: t7 = [0x41]a
	t7.Add(&t3, &t12)
	// Step 16: t0 = [0x4b]a
	t0.Add(&t1, &t7)
	// Step 17: t13 = [0x4d]a
	t13.Add(&t4, &t0)
	// Step 18: t10 = [0x4f]a
	t10.Add(&t4, &t13)
	// Step 19: t11 = [0x53]a
	t11.Add(&t20, &t10)
	// Step 20: t14 = [0x67]a
	t14.Add(&t5, &t11)
	// Step 21: t2 = [0x6f]a
	t2.Add(&t16, &t17)
	// Step 22: t8 = [0x81]a
	t8.Add(&t8, &t2)
	// Step 23: t15 = [0x83]a
	t15.Add(&t4, &t8)
	// Step 24: t5 = [0x97]a
	t5.Add(&t5, &t15)
	// Step 25: t18 = [0x99]a
	t18.Add(&t4, &t5)
	// Step 26: t19 = [0x9d]a
	t19.Add(&t20, &t18)
	// Step 27: t17 = [0xd7]a
	t17.Add(&t17, &t19)
	// Step 28: t20 = [0xdb]a
	t20.Add(&t20, &t17)
	// Step 29: t1 = [0xe5]a
	t1.Add(&t1, &t20)
	// Step 30: t4 = [0xe7]a
	t4.Add(&t4, &t1)
	// Step 31: t21 = [0xe700]a
	t21.Double(&t4)
	for range 7 {
		t21.Double(&t21)
	}
	// Step 32: t20 = [0xe7db]a
	t20.Add(&t20, &t21)
	// Step 33: t20 = [0x1cfb600]a
	for range 9 {
		t20.Double(&t20)
	}
	// Step 34: t19 = [0x1cfb69d]a
	t19.Add(&t19, &t20)
	// Step 35: t19 = [0x39f6d3a00]a
	for range 9 {
		t19.Double(&t19)
	}
	// Step 36: t19 = [0x39f6d3a99]a
	t19.Add(&t18, &t19)
	// Step 37: t19 = [0x73eda753200]a
	for range 9 {
		t19.Double(&t19)
	}
	// Step 38: t18 = [0x73eda753299]a
	t18.Add(&t18, &t19)
	// Step 39: t18 = [0x73eda75329900]a
	for range 8 {
		t18.Double(&t18)
	}
	// Step 40: t17 = [0x73eda753299d7]a
	t17.Add(&t17, &t18)
	// Step 41: t17 = [0x1cfb69d4ca675c0]a
	for range 6 {
		t17.Double(&t17)
	}
	// Step 42: t16 = [0x1cfb69d4ca675f5]a
	t16.Add(&t16, &t17)
	// Step 43: t16 = [0x73eda753299d7d400]a
	for range 10 {
		t16.Double(&t16)
	}
	// Step 44: t15 = [0x73eda753299d7d483]a
	t15.Add(&t15, &t16)
	// Step 45: t15 = [0xe7db4ea6533afa90600]a
	for range 9 {
		t15.Double(&t15)
	}
	// Step 46: t14 = [0xe7db4ea6533afa90667]a
	t14.Add(&t14, &t15)
	// Step 47: t14 = [0xe7db4ea6533afa9066700]a
	for range 8 {
		t14.Double(&t14)
	}
	// Step 48: t14 = [0xe7db4ea6533afa906673b]a
	t14.Add(&t12, &t14)
	// Step 49: t14 = [0xe7db4ea6533afa906673b00]a
	for range 8 {
		t14.Double(&t14)
	}
	// Step 50: t14 = [0xe7db4ea6533afa906673b01]a
	t14.Add(&x, &t14)
	// Step 51: t14 = [0x39f6d3a994cebea4199cec04000]a
	for range 14 {
		t14.Double(&t14)
	}
	// Step 52: t13 = [0x39f6d3a994cebea4199cec0404d]a
	t13.Add(&t13, &t14)
	// Step 53: t13 = [0xe7db4ea6533afa906673b01013400]a
	for range 10 {
		t13.Double(&t13)
	}
	// Step 54: t12 = [0xe7db4ea6533afa906673b0101343b]a
	t12.Add(&t12, &t13)
	// Step 55: t12 = [0x73eda753299d7d483339d80809a1d8000]a
	for range 15 {
		t12.Double(&t12)
	}
	// Step 56: t11 = [0x73eda753299d7d483339d80809a1d8053]a
	t11.Add(&t11, &t12)
	// Step 57: t11 = [0x1cfb69d4ca675f520cce7602026876014c0]a
	for range 6 {
		t11.Double(&t11)
	}
	// Step 58: t11 = [0x1cfb69d4ca675f520cce7602026876014cd]a
	t11.Add(&t6, &t11)
	// Step 59: t11 = [0x1cfb69d4ca675f520cce7602026876014cd000]a
	for range 12 {
		t11.Double(&t11)
	}
	// Step 60: t11 = [0x1cfb69d4ca675f520cce7602026876014cd041]a
	t11.Add(&t7, &t11)
	// Step 61: t11 = [0x39f6d3a994cebea4199cec0404d0ec0299a08200]a
	for range 9 {
		t11.Double(&t11)
	}
	// Step 62: t10 = [0x39f6d3a994cebea4199cec0404d0ec0299a0824f]a
	t10.Add(&t10, &t11)
	// Step 63: t10 = [0x39f6d3a994cebea4199cec0404d0ec0299a0824f00]a
	for range 8 {
		t10.Double(&t10)
	}
	// Step 64: t9 = [0x39f6d3a994cebea4199cec0404d0ec0299a0824f33]a
	t9.Add(&t9, &t10)
	// Step 65: t9 = [0xe7db4ea6533afa906673b0101343b00a6682093ccc00]a
	for range 10 {
		t9.Double(&t9)
	}
	// Step 66: t8 = [0xe7db4ea6533afa906673b0101343b00a6682093ccc81]a
	t8.Add(&t8, &t9)
	// Step 67: t8 = [0x73eda753299d7d483339d80809a1d8053341049e6640800]a
	for range 11 {
		t8.Double(&t8)
	}
	// Step 68: t7 = [0x73eda753299d7d483339d80809a1d8053341049e6640841]a
	t7.Add(&t7, &t8)
	// Step 69: t7 = [0xe7db4ea6533afa906673b0101343b00a6682093ccc810820]a
	for range 5 {
		t7.Double(&t7)
	}
	// Step 70: t6 = [0xe7db4ea6533afa906673b0101343b00a6682093ccc81082d]a
	t6.Add(&t6, &t7)
	// Step 71: t6 = [0xe7db4ea6533afa906673b0101343b00a6682093ccc81082d000]a
	for range 12 {
		t6.Double(&t6)
	}
	// Step 72: t5 = [0xe7db4ea6533afa906673b0101343b00a6682093ccc81082d097]a
	t5.Add(&t5, &t6)
	// Step 73: t5 = [0xe7db4ea6533afa906673b0101343b00a6682093ccc81082d097000]a
	for range 12 {
		t5.Double(&t5)
	}
	// Step 74: t5 = [0xe7db4ea6533afa906673b0101343b00a6682093ccc81082d0970e5]a
	t5.Add(&t1, &t5)
	// Step 75: t5 = [0xe7db4ea6533afa906673b0101343b00a6682093ccc81082d0970e500]a
	for range 8 {
		t5.Double(&t5)
	}
	// Step 76: t4 = [0xe7db4ea6533afa906673b0101343b00a6682093ccc81082d0970e5e7]a
	t4.Add(&t4, &t5)
	// Step 77: t3 = [0xe7db4ea6533afa906673b0101343b00a6682093ccc81082d0970e5ed]a
	t3.Add(&t3, &t4)
	// Step 78: t3 = [0xe7db4ea6533afa906673b0101343b00a6682093ccc81082d0970e5ed00]a
	for range 8 {
		t3.Double(&t3)
	}
	// Step 79: t2 = [0xe7db4ea6533afa906673b0101343b00a6682093ccc81082d0970e5ed6f]a
	t2.Add(&t2, &t3)
	// Step 80: t2 = [0x1cfb69d4ca675f520cce7602026876014cd0412799902105a12e1cbdade00]a
	for range 9 {
		t2.Double(&t2)
	}
	// Step 81: t1 = [0x1cfb69d4ca675f520cce7602026876014cd0412799902105a12e1cbdadee5]a
	t1.Add(&t1, &t2)
	// Step 82: t1 = [0xe7db4ea6533afa906673b0101343b00a6682093ccc81082d0970e5ed6f7280]a
	for range 7 {
		t1.Double(&t1)
	}
	// Step 83: t0 = [0xe7db4ea6533afa906673b0101343b00a6682093ccc81082d0970e5ed6f72cb]a
	t0.Add(&t0, &t1)
	// Step 84: t0 = [0xe7db4ea6533afa906673b0101343b00a6682093ccc81082d0970e5ed6f72cb0]a
	for range 4 {
		t0.Double(&t0)
	}
	// Step 85: z = [0xe7db4ea6533afa906673b0101343b00a6682093ccc81082d0970e5ed6f72cb7]a
	z.Add(&z, &t0)

	return p.Set(&z)
}

// Code generated by internal/generator END mulByOrder. DO NOT EDIT.

// Code generated by internal/generator BEGIN mulByCofactor. DO NOT EDIT.

// mulByCofactor computes p = [h]a and returns p, where
//
//	h = 8
//
// is the cofactor. The result lies in the prime-order subgroup for any
// curve point a. The chain was generated with github.com/mmcloughlin/addchain v0.4.0.
//
// Addition chain:
//
//	return  1 << 3
//
// Operations: 3 point doublings, 0 point additions.
func (p *PointProj) mulByCofactor(a *PointProj) *PointProj {
	// x is a working copy of the input so that p may alias a.
	var x, z PointProj
	x.Set(a)

	// Step 1: z = [0x8]a
	z.Double(&x)
	for range 2 {
		z.Double(&z)
	}

	return p.Set(&z)
}

// Code generated by internal/generator END mulByCofactor. DO NOT EDIT.
