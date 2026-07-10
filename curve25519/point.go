package curve25519

import (
	"math/big"
	"math/bits"

	fp "github.com/yelhousni/divide-and-pair/curve25519/fp"
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

// ScalarMultiplication computes [scalar]*p1 using projective coordinates internally.
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

// setInfinity sets p to the identity (0:1:1)
func (p *PointProj) setInfinity() *PointProj {
	p.X.SetZero()
	p.Y.SetOne()
	p.Z.SetOne()
	return p
}

// Neg negates a projective point: -(X:Y:Z) = (-X:Y:Z)
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
	p.X.Mul(&H, &I).
		Sub(&p.X, &C).
		Sub(&p.X, &D).
		Mul(&p.X, &A).
		Mul(&p.X, &F)
	// Y3 = A*G*(D - a*C)
	var aC fp.Element
	aC.Set(&C)
	mulByA(&aC) // aC = a*C
	p.Y.Sub(&D, &aC).
		Mul(&p.Y, &A).
		Mul(&p.Y, &G)
	p.Z.Mul(&F, &G)

	return p
}

// Double doubles a projective point.
// https://hyperelliptic.org/EFD/g1p/auto-twisted-projective.html#doubling-dbl-2008-bbjlp
func (p *PointProj) Double(p1 *PointProj) *PointProj {
	var B, C, D, E, F, H, J fp.Element

	B.Add(&p1.X, &p1.Y).Square(&B)
	C.Square(&p1.X)
	D.Square(&p1.Y)
	E.Set(&C)
	mulByA(&E) // E = a*C
	F.Add(&E, &D)
	H.Square(&p1.Z)
	J.Sub(&F, &H).Sub(&J, &H) // J = F - 2H
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
//	ℓ = 7237005577332262213973186563042994240857116359379907606001950938285454250989
//	  = 2^252 + 27742317777372353535851937790883648493
//
// is the 253-bit prime subgroup order. The scalar multiplication is unrolled
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
//	_110      = _10 + _100
//	_1000     = _10 + _110
//	_1011     = _11 + _1000
//	_10000    = 2*_1000
//	_100000   = 2*_10000
//	_100110   = _110 + _100000
//	_1000000  = 2*_100000
//	_1010000  = _10000 + _1000000
//	_1010011  = _11 + _1010000
//	_1100011  = _10000 + _1010011
//	_1100111  = _100 + _1100011
//	_1101011  = _100 + _1100111
//	_10010011 = _1000000 + _1010011
//	_10010111 = _100 + _10010011
//	_10111101 = _100110 + _10010111
//	_11010011 = _1000000 + _10010011
//	_11100111 = _1010000 + _10010111
//	_11101101 = _110 + _11100111
//	_11110101 = _1000 + _11101101
//	i160      = ((_1011 + _11110101) << 126 + _1010011) << 9 + _11110101
//	i179      = ((_10 + i160) << 7 + _1100111) << 9 + _11110101
//	i209      = ((i179 << 11 + _10111101) << 8 + _11100111) << 9
//	i232      = ((_1101011 + i209) << 6 + _1011) << 14 + _10010011
//	i263      = ((i232 << 10 + _1100011) << 9 + _10010111) << 10
//	return      ((_11110101 + i263) << 8 + _11010011) << 8 + _11101101
//
// Operations: 248 point doublings, 34 point additions.
func (p *PointProj) mulByOrder(a *PointProj) *PointProj {
	// x is a working copy of the input so that p may alias a.
	var x, z PointProj
	var t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12 PointProj
	x.Set(a)

	// Step 1: t10 = [0x2]a
	t10.Double(&x)
	// Step 2: t4 = [0x3]a
	t4.Add(&x, &t10)
	// Step 3: t2 = [0x4]a
	t2.Add(&x, &t4)
	// Step 4: z = [0x6]a
	z.Add(&t10, &t2)
	// Step 5: t1 = [0x8]a
	t1.Add(&t10, &z)
	// Step 6: t5 = [0xb]a
	t5.Add(&t4, &t1)
	// Step 7: t3 = [0x10]a
	t3.Double(&t1)
	// Step 8: t0 = [0x20]a
	t0.Double(&t3)
	// Step 9: t8 = [0x26]a
	t8.Add(&z, &t0)
	// Step 10: t0 = [0x40]a
	t0.Double(&t0)
	// Step 11: t7 = [0x50]a
	t7.Add(&t3, &t0)
	// Step 12: t11 = [0x53]a
	t11.Add(&t4, &t7)
	// Step 13: t3 = [0x63]a
	t3.Add(&t3, &t11)
	// Step 14: t9 = [0x67]a
	t9.Add(&t2, &t3)
	// Step 15: t6 = [0x6b]a
	t6.Add(&t2, &t9)
	// Step 16: t4 = [0x93]a
	t4.Add(&t0, &t11)
	// Step 17: t2 = [0x97]a
	t2.Add(&t2, &t4)
	// Step 18: t8 = [0xbd]a
	t8.Add(&t8, &t2)
	// Step 19: t0 = [0xd3]a
	t0.Add(&t0, &t4)
	// Step 20: t7 = [0xe7]a
	t7.Add(&t7, &t2)
	// Step 21: z = [0xed]a
	z.Add(&z, &t7)
	// Step 22: t1 = [0xf5]a
	t1.Add(&t1, &z)
	// Step 23: t12 = [0x100]a
	t12.Add(&t5, &t1)
	// Step 24: t12 = [0x4000000000000000000000000000000000]a
	for range 126 {
		t12.Double(&t12)
	}
	// Step 25: t11 = [0x4000000000000000000000000000000053]a
	t11.Add(&t11, &t12)
	// Step 26: t11 = [0x80000000000000000000000000000000a600]a
	for range 9 {
		t11.Double(&t11)
	}
	// Step 27: t11 = [0x80000000000000000000000000000000a6f5]a
	t11.Add(&t1, &t11)
	// Step 28: t10 = [0x80000000000000000000000000000000a6f7]a
	t10.Add(&t10, &t11)
	// Step 29: t10 = [0x40000000000000000000000000000000537b80]a
	for range 7 {
		t10.Double(&t10)
	}
	// Step 30: t9 = [0x40000000000000000000000000000000537be7]a
	t9.Add(&t9, &t10)
	// Step 31: t9 = [0x80000000000000000000000000000000a6f7ce00]a
	for range 9 {
		t9.Double(&t9)
	}
	// Step 32: t9 = [0x80000000000000000000000000000000a6f7cef5]a
	t9.Add(&t1, &t9)
	// Step 33: t9 = [0x40000000000000000000000000000000537be77a800]a
	for range 11 {
		t9.Double(&t9)
	}
	// Step 34: t8 = [0x40000000000000000000000000000000537be77a8bd]a
	t8.Add(&t8, &t9)
	// Step 35: t8 = [0x40000000000000000000000000000000537be77a8bd00]a
	for range 8 {
		t8.Double(&t8)
	}
	// Step 36: t7 = [0x40000000000000000000000000000000537be77a8bde7]a
	t7.Add(&t7, &t8)
	// Step 37: t7 = [0x80000000000000000000000000000000a6f7cef517bce00]a
	for range 9 {
		t7.Double(&t7)
	}
	// Step 38: t6 = [0x80000000000000000000000000000000a6f7cef517bce6b]a
	t6.Add(&t6, &t7)
	// Step 39: t6 = [0x2000000000000000000000000000000029bdf3bd45ef39ac0]a
	for range 6 {
		t6.Double(&t6)
	}
	// Step 40: t5 = [0x2000000000000000000000000000000029bdf3bd45ef39acb]a
	t5.Add(&t5, &t6)
	// Step 41: t5 = [0x80000000000000000000000000000000a6f7cef517bce6b2c000]a
	for range 14 {
		t5.Double(&t5)
	}
	// Step 42: t4 = [0x80000000000000000000000000000000a6f7cef517bce6b2c093]a
	t4.Add(&t4, &t5)
	// Step 43: t4 = [0x2000000000000000000000000000000029bdf3bd45ef39acb024c00]a
	for range 10 {
		t4.Double(&t4)
	}
	// Step 44: t3 = [0x2000000000000000000000000000000029bdf3bd45ef39acb024c63]a
	t3.Add(&t3, &t4)
	// Step 45: t3 = [0x40000000000000000000000000000000537be77a8bde735960498c600]a
	for range 9 {
		t3.Double(&t3)
	}
	// Step 46: t2 = [0x40000000000000000000000000000000537be77a8bde735960498c697]a
	t2.Add(&t2, &t3)
	// Step 47: t2 = [0x1000000000000000000000000000000014def9dea2f79cd65812631a5c00]a
	for range 10 {
		t2.Double(&t2)
	}
	// Step 48: t1 = [0x1000000000000000000000000000000014def9dea2f79cd65812631a5cf5]a
	t1.Add(&t1, &t2)
	// Step 49: t1 = [0x1000000000000000000000000000000014def9dea2f79cd65812631a5cf500]a
	for range 8 {
		t1.Double(&t1)
	}
	// Step 50: t0 = [0x1000000000000000000000000000000014def9dea2f79cd65812631a5cf5d3]a
	t0.Add(&t0, &t1)
	// Step 51: t0 = [0x1000000000000000000000000000000014def9dea2f79cd65812631a5cf5d300]a
	for range 8 {
		t0.Double(&t0)
	}
	// Step 52: z = [0x1000000000000000000000000000000014def9dea2f79cd65812631a5cf5d3ed]a
	z.Add(&z, &t0)

	return p.Set(&z)
}

// Code generated by internal/generator END mulByOrder. DO NOT EDIT.
