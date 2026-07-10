package fourq

import (
	"math/big"
	"math/bits"

	fp2 "github.com/yelhousni/divide-and-pair/fourq/fp2"
)

// PointAffine is a point on FourQ in affine coordinates over Fp2.
type PointAffine struct {
	X, Y fp2.E2
}

// Set sets p to p1 and returns p.
func (p *PointAffine) Set(p1 *PointAffine) *PointAffine {
	p.X.Set(&p1.X)
	p.Y.Set(&p1.Y)
	return p
}

// Equal returns true if p and p1 are equal.
func (p *PointAffine) Equal(p1 *PointAffine) bool {
	return p.X.Equal(&p1.X) && p.Y.Equal(&p1.Y)
}

// SetInfinity sets p to the identity point, encoded as (0, 1).
func (p *PointAffine) SetInfinity() *PointAffine {
	p.X.SetZero()
	p.Y.SetOne()
	return p
}

// IsZero returns true if p is the identity point, encoded as (0, 1).
func (p *PointAffine) IsZero() bool {
	return p.X.IsZero() && p.Y.IsOne()
}

// IsOnCurve checks: a*x² + y² = 1 + d*x²*y² with a = -1.
func (p *PointAffine) IsOnCurve() bool {
	initOnce.Do(initCurveParams)

	var x2, y2, lhs, rhs, tmp fp2.E2

	x2.Square(&p.X)
	y2.Square(&p.Y)

	// lhs = a*x² + y² = -x² + y²
	lhs.Set(&x2)
	mulByA(&lhs)
	lhs.Add(&lhs, &y2)

	// rhs = 1 + d*x²*y²
	tmp.Mul(&x2, &y2)
	tmp.Mul(&tmp, &curveParams.D)
	rhs.SetOne()
	rhs.Add(&rhs, &tmp)

	return lhs.Equal(&rhs)
}

// Neg negates a point: -(x, y) = (-x, y).
func (p *PointAffine) Neg(p1 *PointAffine) *PointAffine {
	p.X.Neg(&p1.X)
	p.Y.Set(&p1.Y)
	return p
}

// Add adds two affine points using unified addition for twisted Edwards.
func (p *PointAffine) Add(p1, p2 *PointAffine) *PointAffine {
	initOnce.Do(initCurveParams)

	var x1y2, y1x2, y1y2, x1x2, dxy fp2.E2

	x1y2.Mul(&p1.X, &p2.Y)
	y1x2.Mul(&p1.Y, &p2.X)
	y1y2.Mul(&p1.Y, &p2.Y)
	x1x2.Mul(&p1.X, &p2.X)
	dxy.Mul(&curveParams.D, &x1x2)
	dxy.Mul(&dxy, &y1y2)

	var numX, numY, denX, denY fp2.E2
	numX.Add(&x1y2, &y1x2)

	// numY = y1y2 - a*x1x2 = y1y2 + x1x2 (since a = -1)
	numY.Add(&y1y2, &x1x2)

	denX.SetOne()
	denX.Add(&denX, &dxy)
	denY.SetOne()
	denY.Sub(&denY, &dxy)

	denX.Inverse(&denX)
	denY.Inverse(&denY)

	p.X.Mul(&numX, &denX)
	p.Y.Mul(&numY, &denY)
	return p
}

// Double doubles a point.
func (p *PointAffine) Double(p1 *PointAffine) *PointAffine {
	return p.Add(p1, p1)
}

// ScalarMultiplication computes [scalar]*p1 using projective coordinates internally.
func (p *PointAffine) ScalarMultiplication(p1 *PointAffine, scalar *big.Int) *PointAffine {
	var proj, res PointProj
	proj.FromAffine(p1)
	res.ScalarMultiplication(&proj, scalar)
	res.ToAffine(p)
	return p
}

// --- Projective coordinates over Fp2 ---

// PointProj point in projective coordinates (X:Y:Z) with x=X/Z, y=Y/Z over Fp2.
type PointProj struct {
	X, Y, Z fp2.E2
}

// Set sets p to p1.
func (p *PointProj) Set(p1 *PointProj) *PointProj {
	p.X.Set(&p1.X)
	p.Y.Set(&p1.Y)
	p.Z.Set(&p1.Z)
	return p
}

// FromAffine converts affine to projective.
func (p *PointProj) FromAffine(a *PointAffine) *PointProj {
	p.X.Set(&a.X)
	p.Y.Set(&a.Y)
	p.Z.SetOne()
	return p
}

// ToAffine converts projective to affine (single Fp2 inversion).
func (p *PointProj) ToAffine(a *PointAffine) *PointAffine {
	var zInv fp2.E2
	zInv.Inverse(&p.Z)
	a.X.Mul(&p.X, &zInv)
	a.Y.Mul(&p.Y, &zInv)
	return a
}

// IsZero returns true if p is the identity (0:1:1).
func (p *PointProj) IsZero() bool {
	return p.X.IsZero() && p.Y.Equal(&p.Z)
}

func (p *PointProj) setInfinity() *PointProj {
	p.X.SetZero()
	p.Y.SetOne()
	p.Z.SetOne()
	return p
}

// Neg negates a projective point: -(X:Y:Z) = (-X:Y:Z).
func (p *PointProj) Neg(p1 *PointProj) *PointProj {
	p.X.Neg(&p1.X)
	p.Y.Set(&p1.Y)
	p.Z.Set(&p1.Z)
	return p
}

// Add adds two projective points using the unified addition formula.
// https://hyperelliptic.org/EFD/g1p/auto-twisted-projective.html#addition-add-2008-bbjlp
func (p *PointProj) Add(p1, p2 *PointProj) *PointProj {
	initOnce.Do(initCurveParams)

	var A, B, C, D, E, F, G, H, I fp2.E2
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
	p.X.Mul(&H, &I).
		Sub(&p.X, &C).
		Sub(&p.X, &D).
		Mul(&p.X, &A).
		Mul(&p.X, &F)
	var aC fp2.E2
	aC.Set(&C)
	mulByA(&aC)
	p.Y.Sub(&D, &aC).
		Mul(&p.Y, &A).
		Mul(&p.Y, &G)
	p.Z.Mul(&F, &G)

	return p
}

// Double doubles a projective point.
// https://hyperelliptic.org/EFD/g1p/auto-twisted-projective.html#doubling-dbl-2008-bbjlp
func (p *PointProj) Double(p1 *PointProj) *PointProj {
	var B, C, D, E, F, H, J fp2.E2

	B.Add(&p1.X, &p1.Y).Square(&B)
	C.Square(&p1.X)
	D.Square(&p1.Y)
	E.Set(&C)
	mulByA(&E)
	F.Add(&E, &D)
	H.Square(&p1.Z)
	J.Sub(&F, &H).Sub(&J, &H)
	p.X.Sub(&B, &C).
		Sub(&p.X, &D).
		Mul(&p.X, &J)
	p.Y.Sub(&E, &D).Mul(&p.Y, &F)
	p.Z.Mul(&F, &J)

	return p
}

// ScalarMultiplication computes [scalar]*p1 in projective coordinates.
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

// Code generated by internal/generator BEGIN mulByOrder. DO NOT EDIT.

// mulByOrder computes p = [ℓ]a and returns p, where
//
//	ℓ = 73846995687063900142583536357581573884798075859800097461294096333596429543
//
// is the 246-bit prime subgroup order. The scalar multiplication is unrolled
// as a short addition chain generated with github.com/mmcloughlin/addchain v0.4.0.
//
// The unified Add/Double formulas are valid for arbitrary curve points
// (including the identity and points outside the ℓ-torsion subgroup), so the
// result is the identity if and only if a is in the prime-order subgroup.
//
// Addition chain:
//
//	_10       = 2*1
//	_11       = 1 + _10
//	_100      = 1 + _11
//	_111      = _11 + _100
//	_1010     = _11 + _111
//	_1100     = _10 + _1010
//	_10110    = _1010 + _1100
//	_11010    = _100 + _10110
//	_110100   = 2*_11010
//	_1001010  = _10110 + _110100
//	_1010110  = _1100 + _1001010
//	_1011010  = _100 + _1010110
//	_1011011  = 1 + _1011010
//	_1110101  = _11010 + _1011011
//	_1111111  = _1010 + _1110101
//	_11001001 = _1001010 + _1111111
//	_11101010 = 2*_1110101
//	i19       = 2*_11101010 + _110100
//	i20       = _1011010 + i19
//	i21       = _1110101 + i20
//	i22       = _1010110 + i21
//	i23       = _11001001 + i22
//	i24       = i20 + i23
//	i25       = _11101010 + i24
//	i26       = _10110 + i25
//	i27       = i21 + i26
//	i28       = _1111111 + i27
//	i29       = i22 + i28
//	i30       = _1011011 + i29
//	i31       = i27 + i30
//	i32       = i26 + i31
//	i33       = i25 + i32
//	i34       = i31 + i33
//	i35       = i24 + i34
//	i36       = i29 + i35
//	i37       = i34 + i36
//	i38       = i19 + i37
//	i39       = i30 + i38
//	i40       = i23 + i39
//	i41       = i28 + i40
//	i106      = ((i40 << 21 + i40) << 21 + i40) << 21
//	i150      = ((i40 + i106) << 21 + i40) << 20 + i36
//	i210      = ((i150 << 14 + i32) << 23 + i33) << 21
//	i247      = ((i41 + i210 + i35) << 16 + i37) << 18
//	i282      = ((i38 + i247) << 16 + i41) << 16 + i39
//	return      i282 << 2 + _11
//
// Operations: 234 point doublings, 51 point additions.
func (p *PointProj) mulByOrder(a *PointProj) *PointProj {
	// x is a working copy of the input so that p may alias a.
	var x, z PointProj
	var t0, t1, t2, t3, t4, t5, t6, t7, t8, t9 PointProj
	x.Set(a)

	// Step 1: t2 = [0x2]a
	t2.Double(&x)
	// Step 2: z = [0x3]a
	z.Add(&x, &t2)
	// Step 3: t0 = [0x4]a
	t0.Add(&x, &z)
	// Step 4: t1 = [0x7]a
	t1.Add(&z, &t0)
	// Step 5: t1 = [0xa]a
	t1.Add(&z, &t1)
	// Step 6: t4 = [0xc]a
	t4.Add(&t2, &t1)
	// Step 7: t6 = [0x16]a
	t6.Add(&t1, &t4)
	// Step 8: t3 = [0x1a]a
	t3.Add(&t0, &t6)
	// Step 9: t2 = [0x34]a
	t2.Double(&t3)
	// Step 10: t5 = [0x4a]a
	t5.Add(&t6, &t2)
	// Step 11: t7 = [0x56]a
	t7.Add(&t4, &t5)
	// Step 12: t4 = [0x5a]a
	t4.Add(&t0, &t7)
	// Step 13: t0 = [0x5b]a
	t0.Add(&x, &t4)
	// Step 14: t3 = [0x75]a
	t3.Add(&t3, &t0)
	// Step 15: t1 = [0x7f]a
	t1.Add(&t1, &t3)
	// Step 16: t8 = [0xc9]a
	t8.Add(&t5, &t1)
	// Step 17: t5 = [0xea]a
	t5.Double(&t3)
	// Step 18: t9 = [0x1d4]a
	t9.Double(&t5)
	// Step 19: t2 = [0x208]a
	t2.Add(&t2, &t9)
	// Step 20: t4 = [0x262]a
	t4.Add(&t4, &t2)
	// Step 21: t3 = [0x2d7]a
	t3.Add(&t3, &t4)
	// Step 22: t7 = [0x32d]a
	t7.Add(&t7, &t3)
	// Step 23: t8 = [0x3f6]a
	t8.Add(&t8, &t7)
	// Step 24: t4 = [0x658]a
	t4.Add(&t4, &t8)
	// Step 25: t5 = [0x742]a
	t5.Add(&t5, &t4)
	// Step 26: t6 = [0x758]a
	t6.Add(&t6, &t5)
	// Step 27: t3 = [0xa2f]a
	t3.Add(&t3, &t6)
	// Step 28: t1 = [0xaae]a
	t1.Add(&t1, &t3)
	// Step 29: t7 = [0xddb]a
	t7.Add(&t7, &t1)
	// Step 30: t0 = [0xe36]a
	t0.Add(&t0, &t7)
	// Step 31: t3 = [0x1865]a
	t3.Add(&t3, &t0)
	// Step 32: t6 = [0x1fbd]a
	t6.Add(&t6, &t3)
	// Step 33: t5 = [0x26ff]a
	t5.Add(&t5, &t6)
	// Step 34: t3 = [0x3f64]a
	t3.Add(&t3, &t5)
	// Step 35: t4 = [0x45bc]a
	t4.Add(&t4, &t3)
	// Step 36: t7 = [0x5397]a
	t7.Add(&t7, &t4)
	// Step 37: t3 = [0x92fb]a
	t3.Add(&t3, &t7)
	// Step 38: t2 = [0x9503]a
	t2.Add(&t2, &t3)
	// Step 39: t0 = [0xa339]a
	t0.Add(&t0, &t2)
	// Step 40: t8 = [0xa72f]a
	t8.Add(&t8, &t0)
	// Step 41: t1 = [0xb1dd]a
	t1.Add(&t1, &t8)
	// Step 42: t9 = [0x14e5e00000]a
	t9.Double(&t8)
	for range 20 {
		t9.Double(&t9)
	}
	// Step 43: t9 = [0x14e5e0a72f]a
	t9.Add(&t8, &t9)
	// Step 44: t9 = [0x29cbc14e5e00000]a
	for range 21 {
		t9.Double(&t9)
	}
	// Step 45: t9 = [0x29cbc14e5e0a72f]a
	t9.Add(&t8, &t9)
	// Step 46: t9 = [0x5397829cbc14e5e00000]a
	for range 21 {
		t9.Double(&t9)
	}
	// Step 47: t9 = [0x5397829cbc14e5e0a72f]a
	t9.Add(&t8, &t9)
	// Step 48: t9 = [0xa72f05397829cbc14e5e00000]a
	for range 21 {
		t9.Double(&t9)
	}
	// Step 49: t8 = [0xa72f05397829cbc14e5e0a72f]a
	t8.Add(&t8, &t9)
	// Step 50: t8 = [0xa72f05397829cbc14e5e0a72f00000]a
	for range 20 {
		t8.Double(&t8)
	}
	// Step 51: t7 = [0xa72f05397829cbc14e5e0a72f05397]a
	t7.Add(&t7, &t8)
	// Step 52: t7 = [0x29cbc14e5e0a72f05397829cbc14e5c000]a
	for range 14 {
		t7.Double(&t7)
	}
	// Step 53: t6 = [0x29cbc14e5e0a72f05397829cbc14e5dfbd]a
	t6.Add(&t6, &t7)
	// Step 54: t6 = [0x14e5e0a72f05397829cbc14e5e0a72efde800000]a
	for range 23 {
		t6.Double(&t6)
	}
	// Step 55: t5 = [0x14e5e0a72f05397829cbc14e5e0a72efde8026ff]a
	t5.Add(&t5, &t6)
	// Step 56: t5 = [0x29cbc14e5e0a72f05397829cbc14e5dfbd004dfe00000]a
	for range 21 {
		t5.Double(&t5)
	}
	// Step 57: t5 = [0x29cbc14e5e0a72f05397829cbc14e5dfbd004dfe0b1dd]a
	t5.Add(&t1, &t5)
	// Step 58: t4 = [0x29cbc14e5e0a72f05397829cbc14e5dfbd004dfe0f799]a
	t4.Add(&t4, &t5)
	// Step 59: t4 = [0x29cbc14e5e0a72f05397829cbc14e5dfbd004dfe0f7990000]a
	for range 16 {
		t4.Double(&t4)
	}
	// Step 60: t3 = [0x29cbc14e5e0a72f05397829cbc14e5dfbd004dfe0f79992fb]a
	t3.Add(&t3, &t4)
	// Step 61: t3 = [0xa72f05397829cbc14e5e0a72f053977ef40137f83de664bec0000]a
	for range 18 {
		t3.Double(&t3)
	}
	// Step 62: t2 = [0xa72f05397829cbc14e5e0a72f053977ef40137f83de664bec9503]a
	t2.Add(&t2, &t3)
	// Step 63: t2 = [0xa72f05397829cbc14e5e0a72f053977ef40137f83de664bec95030000]a
	for range 16 {
		t2.Double(&t2)
	}
	// Step 64: t1 = [0xa72f05397829cbc14e5e0a72f053977ef40137f83de664bec9503b1dd]a
	t1.Add(&t1, &t2)
	// Step 65: t1 = [0xa72f05397829cbc14e5e0a72f053977ef40137f83de664bec9503b1dd0000]a
	for range 16 {
		t1.Double(&t1)
	}
	// Step 66: t0 = [0xa72f05397829cbc14e5e0a72f053977ef40137f83de664bec9503b1dda339]a
	t0.Add(&t0, &t1)
	// Step 67: t0 = [0x29cbc14e5e0a72f05397829cbc14e5dfbd004dfe0f79992fb2540ec7768ce4]a
	for range 2 {
		t0.Double(&t0)
	}
	// Step 68: z = [0x29cbc14e5e0a72f05397829cbc14e5dfbd004dfe0f79992fb2540ec7768ce7]a
	z.Add(&z, &t0)

	return p.Set(&z)
}

// Code generated by internal/generator END mulByOrder. DO NOT EDIT.
