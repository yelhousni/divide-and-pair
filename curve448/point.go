package curve448

import (
	"math/big"
	"math/bits"

	fp "github.com/yelhousni/divide-and-pair/curve448/fp"
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
	numY.Sub(&y1y2, &ax1x2) // y1y2 - a*x1x2 = y1y2 - x1x2 (since a=1)

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
//	ℓ = 181709681073901722637330951972001133588410340171829515070372549795146003961539585716195755291692375963310293709091662304773755859649779
//	  = 2^446 - 13818066809895115352007386748515426880336692474882178609894547503885
//
// is the 446-bit prime subgroup order. The scalar multiplication is unrolled
// as a short addition chain generated with github.com/mmcloughlin/addchain v0.4.0.
//
// The unified Add/Double formulas are valid for arbitrary curve points
// (including the identity and points outside the ℓ-torsion subgroup), so the
// result is the identity if and only if a is in the prime-order subgroup.
//
// Addition chain:
//
//	_10      = 2*1
//	_11      = 1 + _10
//	_100     = 1 + _11
//	_101     = 1 + _100
//	_1001    = _100 + _101
//	_1011    = _10 + _1001
//	_1101    = _10 + _1011
//	_1111    = _10 + _1101
//	_10001   = _10 + _1111
//	_10011   = _10 + _10001
//	_10101   = _10 + _10011
//	_10111   = _10 + _10101
//	_11001   = _10 + _10111
//	_11011   = _10 + _11001
//	_11101   = _10 + _11011
//	_11111   = _10 + _11101
//	_111110  = 2*_11111
//	_1111100 = 2*_111110
//	i24      = _1111100 << 5 + _1111100
//	i41      = (i24 << 4 + _111110) << 11 + i24
//	i73      = (i41 << 4 + _111110) << 26 + i41
//	i129     = i73 << 55 + i73
//	x222     = i129 << 110 + i129 + _11
//	i262     = ((x222 << 6 + _11111) << 7 + _11001) << 6
//	i279     = ((_10001 + i262) << 8 + _11111) << 6 + _10011
//	i298     = ((i279 << 5 + _10001) << 8 + _10011) << 4
//	i312     = ((_1011 + i298) << 6 + _11011) << 5 + _1001
//	i331     = ((i312 << 6 + _1101) << 6 + _11101) << 5
//	i343     = ((_10101 + i331) << 5 + _10001) << 4 + _1011
//	i365     = ((i343 << 5 + _1001) << 7 + 1) << 8
//	i375     = 2*((_1011 + i365) << 6 + _11001) + 1
//	i396     = ((i375 << 9 + _10011) << 4 + _1001) << 6
//	i411     = ((_10001 + i396) << 5 + _10111) << 7 + _1011
//	i431     = ((i411 << 7 + _1111) << 6 + _10101) << 5
//	i444     = ((_1001 + i431) << 8 + _11011) << 2 + _11
//	i464     = ((i444 << 5 + _11) << 7 + _101) << 6
//	i478     = ((_1001 + i464) << 6 + _10101) << 5 + _1101
//	i498     = ((i478 << 3 + _11) << 9 + _10001) << 6
//	return     (_1111 + i498) << 4 + _11
//
// Operations: 442 point doublings, 62 point additions.
func (p *PointProj) mulByOrder(a *PointProj) *PointProj {
	// x is a working copy of the input so that p may alias a.
	var x, z PointProj
	var t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15 PointProj
	x.Set(a)

	// Step 1: t12 = [0x2]a
	t12.Double(&x)
	// Step 2: z = [0x3]a
	z.Add(&x, &t12)
	// Step 3: t0 = [0x4]a
	t0.Add(&x, &z)
	// Step 4: t5 = [0x5]a
	t5.Add(&x, &t0)
	// Step 5: t4 = [0x9]a
	t4.Add(&t0, &t5)
	// Step 6: t7 = [0xb]a
	t7.Add(&t12, &t4)
	// Step 7: t2 = [0xd]a
	t2.Add(&t12, &t7)
	// Step 8: t0 = [0xf]a
	t0.Add(&t12, &t2)
	// Step 9: t1 = [0x11]a
	t1.Add(&t12, &t0)
	// Step 10: t9 = [0x13]a
	t9.Add(&t12, &t1)
	// Step 11: t3 = [0x15]a
	t3.Add(&t12, &t9)
	// Step 12: t8 = [0x17]a
	t8.Add(&t12, &t3)
	// Step 13: t10 = [0x19]a
	t10.Add(&t12, &t8)
	// Step 14: t6 = [0x1b]a
	t6.Add(&t12, &t10)
	// Step 15: t11 = [0x1d]a
	t11.Add(&t12, &t6)
	// Step 16: t12 = [0x1f]a
	t12.Add(&t12, &t11)
	// Step 17: t14 = [0x3e]a
	t14.Double(&t12)
	// Step 18: t13 = [0x7c]a
	t13.Double(&t14)
	// Step 19: t15 = [0xf80]a
	t15.Double(&t13)
	for range 4 {
		t15.Double(&t15)
	}
	// Step 20: t13 = [0xffc]a
	t13.Add(&t13, &t15)
	// Step 21: t15 = [0xffc0]a
	t15.Double(&t13)
	for range 3 {
		t15.Double(&t15)
	}
	// Step 22: t15 = [0xfffe]a
	t15.Add(&t14, &t15)
	// Step 23: t15 = [0x7fff000]a
	for range 11 {
		t15.Double(&t15)
	}
	// Step 24: t13 = [0x7fffffc]a
	t13.Add(&t13, &t15)
	// Step 25: t15 = [0x7fffffc0]a
	t15.Double(&t13)
	for range 3 {
		t15.Double(&t15)
	}
	// Step 26: t14 = [0x7ffffffe]a
	t14.Add(&t14, &t15)
	// Step 27: t14 = [0x1fffffff8000000]a
	for range 26 {
		t14.Double(&t14)
	}
	// Step 28: t13 = [0x1fffffffffffffc]a
	t13.Add(&t13, &t14)
	// Step 29: t14 = [0xfffffffffffffe00000000000000]a
	t14.Double(&t13)
	for range 54 {
		t14.Double(&t14)
	}
	// Step 30: t13 = [0xfffffffffffffffffffffffffffc]a
	t13.Add(&t13, &t14)
	// Step 31: t14 = [0x3fffffffffffffffffffffffffff0000000000000000000000000000]a
	t14.Double(&t13)
	for range 109 {
		t14.Double(&t14)
	}
	// Step 32: t13 = [0x3ffffffffffffffffffffffffffffffffffffffffffffffffffffffc]a
	t13.Add(&t13, &t14)
	// Step 33: t13 = [0x3fffffffffffffffffffffffffffffffffffffffffffffffffffffff]a
	t13.Add(&z, &t13)
	// Step 34: t13 = [0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffc0]a
	for range 6 {
		t13.Double(&t13)
	}
	// Step 35: t13 = [0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffdf]a
	t13.Add(&t12, &t13)
	// Step 36: t13 = [0x7ffffffffffffffffffffffffffffffffffffffffffffffffffffffef80]a
	for range 7 {
		t13.Double(&t13)
	}
	// Step 37: t13 = [0x7ffffffffffffffffffffffffffffffffffffffffffffffffffffffef99]a
	t13.Add(&t10, &t13)
	// Step 38: t13 = [0x1fffffffffffffffffffffffffffffffffffffffffffffffffffffffbe640]a
	for range 6 {
		t13.Double(&t13)
	}
	// Step 39: t13 = [0x1fffffffffffffffffffffffffffffffffffffffffffffffffffffffbe651]a
	t13.Add(&t1, &t13)
	// Step 40: t13 = [0x1fffffffffffffffffffffffffffffffffffffffffffffffffffffffbe65100]a
	for range 8 {
		t13.Double(&t13)
	}
	// Step 41: t12 = [0x1fffffffffffffffffffffffffffffffffffffffffffffffffffffffbe6511f]a
	t12.Add(&t12, &t13)
	// Step 42: t12 = [0x7ffffffffffffffffffffffffffffffffffffffffffffffffffffffef99447c0]a
	for range 6 {
		t12.Double(&t12)
	}
	// Step 43: t12 = [0x7ffffffffffffffffffffffffffffffffffffffffffffffffffffffef99447d3]a
	t12.Add(&t9, &t12)
	// Step 44: t12 = [0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffdf3288fa60]a
	for range 5 {
		t12.Double(&t12)
	}
	// Step 45: t12 = [0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffdf3288fa71]a
	t12.Add(&t1, &t12)
	// Step 46: t12 = [0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffdf3288fa7100]a
	for range 8 {
		t12.Double(&t12)
	}
	// Step 47: t12 = [0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffdf3288fa7113]a
	t12.Add(&t9, &t12)
	// Step 48: t12 = [0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffdf3288fa71130]a
	for range 4 {
		t12.Double(&t12)
	}
	// Step 49: t12 = [0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffdf3288fa7113b]a
	t12.Add(&t7, &t12)
	// Step 50: t12 = [0x3fffffffffffffffffffffffffffffffffffffffffffffffffffffff7cca23e9c44ec0]a
	for range 6 {
		t12.Double(&t12)
	}
	// Step 51: t12 = [0x3fffffffffffffffffffffffffffffffffffffffffffffffffffffff7cca23e9c44edb]a
	t12.Add(&t6, &t12)
	// Step 52: t12 = [0x7ffffffffffffffffffffffffffffffffffffffffffffffffffffffef99447d3889db60]a
	for range 5 {
		t12.Double(&t12)
	}
	// Step 53: t12 = [0x7ffffffffffffffffffffffffffffffffffffffffffffffffffffffef99447d3889db69]a
	t12.Add(&t4, &t12)
	// Step 54: t12 = [0x1fffffffffffffffffffffffffffffffffffffffffffffffffffffffbe6511f4e2276da40]a
	for range 6 {
		t12.Double(&t12)
	}
	// Step 55: t12 = [0x1fffffffffffffffffffffffffffffffffffffffffffffffffffffffbe6511f4e2276da4d]a
	t12.Add(&t2, &t12)
	// Step 56: t12 = [0x7ffffffffffffffffffffffffffffffffffffffffffffffffffffffef99447d3889db69340]a
	for range 6 {
		t12.Double(&t12)
	}
	// Step 57: t11 = [0x7ffffffffffffffffffffffffffffffffffffffffffffffffffffffef99447d3889db6935d]a
	t11.Add(&t11, &t12)
	// Step 58: t11 = [0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffdf3288fa7113b6d26ba0]a
	for range 5 {
		t11.Double(&t11)
	}
	// Step 59: t11 = [0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffdf3288fa7113b6d26bb5]a
	t11.Add(&t3, &t11)
	// Step 60: t11 = [0x1fffffffffffffffffffffffffffffffffffffffffffffffffffffffbe6511f4e2276da4d76a0]a
	for range 5 {
		t11.Double(&t11)
	}
	// Step 61: t11 = [0x1fffffffffffffffffffffffffffffffffffffffffffffffffffffffbe6511f4e2276da4d76b1]a
	t11.Add(&t1, &t11)
	// Step 62: t11 = [0x1fffffffffffffffffffffffffffffffffffffffffffffffffffffffbe6511f4e2276da4d76b10]a
	for range 4 {
		t11.Double(&t11)
	}
	// Step 63: t11 = [0x1fffffffffffffffffffffffffffffffffffffffffffffffffffffffbe6511f4e2276da4d76b1b]a
	t11.Add(&t7, &t11)
	// Step 64: t11 = [0x3fffffffffffffffffffffffffffffffffffffffffffffffffffffff7cca23e9c44edb49aed6360]a
	for range 5 {
		t11.Double(&t11)
	}
	// Step 65: t11 = [0x3fffffffffffffffffffffffffffffffffffffffffffffffffffffff7cca23e9c44edb49aed6369]a
	t11.Add(&t4, &t11)
	// Step 66: t11 = [0x1fffffffffffffffffffffffffffffffffffffffffffffffffffffffbe6511f4e2276da4d76b1b480]a
	for range 7 {
		t11.Double(&t11)
	}
	// Step 67: t11 = [0x1fffffffffffffffffffffffffffffffffffffffffffffffffffffffbe6511f4e2276da4d76b1b481]a
	t11.Add(&x, &t11)
	// Step 68: t11 = [0x1fffffffffffffffffffffffffffffffffffffffffffffffffffffffbe6511f4e2276da4d76b1b48100]a
	for range 8 {
		t11.Double(&t11)
	}
	// Step 69: t11 = [0x1fffffffffffffffffffffffffffffffffffffffffffffffffffffffbe6511f4e2276da4d76b1b4810b]a
	t11.Add(&t7, &t11)
	// Step 70: t11 = [0x7ffffffffffffffffffffffffffffffffffffffffffffffffffffffef99447d3889db6935dac6d2042c0]a
	for range 6 {
		t11.Double(&t11)
	}
	// Step 71: t10 = [0x7ffffffffffffffffffffffffffffffffffffffffffffffffffffffef99447d3889db6935dac6d2042d9]a
	t10.Add(&t10, &t11)
	// Step 72: t10 = [0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffdf3288fa7113b6d26bb58da4085b2]a
	t10.Double(&t10)
	// Step 73: t10 = [0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffdf3288fa7113b6d26bb58da4085b3]a
	t10.Add(&x, &t10)
	// Step 74: t10 = [0x1fffffffffffffffffffffffffffffffffffffffffffffffffffffffbe6511f4e2276da4d76b1b4810b6600]a
	for range 9 {
		t10.Double(&t10)
	}
	// Step 75: t9 = [0x1fffffffffffffffffffffffffffffffffffffffffffffffffffffffbe6511f4e2276da4d76b1b4810b6613]a
	t9.Add(&t9, &t10)
	// Step 76: t9 = [0x1fffffffffffffffffffffffffffffffffffffffffffffffffffffffbe6511f4e2276da4d76b1b4810b66130]a
	for range 4 {
		t9.Double(&t9)
	}
	// Step 77: t9 = [0x1fffffffffffffffffffffffffffffffffffffffffffffffffffffffbe6511f4e2276da4d76b1b4810b66139]a
	t9.Add(&t4, &t9)
	// Step 78: t9 = [0x7ffffffffffffffffffffffffffffffffffffffffffffffffffffffef99447d3889db6935dac6d2042d984e40]a
	for range 6 {
		t9.Double(&t9)
	}
	// Step 79: t9 = [0x7ffffffffffffffffffffffffffffffffffffffffffffffffffffffef99447d3889db6935dac6d2042d984e51]a
	t9.Add(&t1, &t9)
	// Step 80: t9 = [0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffdf3288fa7113b6d26bb58da4085b309ca20]a
	for range 5 {
		t9.Double(&t9)
	}
	// Step 81: t8 = [0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffdf3288fa7113b6d26bb58da4085b309ca37]a
	t8.Add(&t8, &t9)
	// Step 82: t8 = [0x7ffffffffffffffffffffffffffffffffffffffffffffffffffffffef99447d3889db6935dac6d2042d984e51b80]a
	for range 7 {
		t8.Double(&t8)
	}
	// Step 83: t7 = [0x7ffffffffffffffffffffffffffffffffffffffffffffffffffffffef99447d3889db6935dac6d2042d984e51b8b]a
	t7.Add(&t7, &t8)
	// Step 84: t7 = [0x3fffffffffffffffffffffffffffffffffffffffffffffffffffffff7cca23e9c44edb49aed63690216cc2728dc580]a
	for range 7 {
		t7.Double(&t7)
	}
	// Step 85: t7 = [0x3fffffffffffffffffffffffffffffffffffffffffffffffffffffff7cca23e9c44edb49aed63690216cc2728dc58f]a
	t7.Add(&t0, &t7)
	// Step 86: t7 = [0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffdf3288fa7113b6d26bb58da4085b309ca37163c0]a
	for range 6 {
		t7.Double(&t7)
	}
	// Step 87: t7 = [0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffdf3288fa7113b6d26bb58da4085b309ca37163d5]a
	t7.Add(&t3, &t7)
	// Step 88: t7 = [0x1fffffffffffffffffffffffffffffffffffffffffffffffffffffffbe6511f4e2276da4d76b1b4810b6613946e2c7aa0]a
	for range 5 {
		t7.Double(&t7)
	}
	// Step 89: t7 = [0x1fffffffffffffffffffffffffffffffffffffffffffffffffffffffbe6511f4e2276da4d76b1b4810b6613946e2c7aa9]a
	t7.Add(&t4, &t7)
	// Step 90: t7 = [0x1fffffffffffffffffffffffffffffffffffffffffffffffffffffffbe6511f4e2276da4d76b1b4810b6613946e2c7aa900]a
	for range 8 {
		t7.Double(&t7)
	}
	// Step 91: t6 = [0x1fffffffffffffffffffffffffffffffffffffffffffffffffffffffbe6511f4e2276da4d76b1b4810b6613946e2c7aa91b]a
	t6.Add(&t6, &t7)
	// Step 92: t6 = [0x7ffffffffffffffffffffffffffffffffffffffffffffffffffffffef99447d3889db6935dac6d2042d984e51b8b1eaa46c]a
	for range 2 {
		t6.Double(&t6)
	}
	// Step 93: t6 = [0x7ffffffffffffffffffffffffffffffffffffffffffffffffffffffef99447d3889db6935dac6d2042d984e51b8b1eaa46f]a
	t6.Add(&z, &t6)
	// Step 94: t6 = [0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffdf3288fa7113b6d26bb58da4085b309ca37163d548de0]a
	for range 5 {
		t6.Double(&t6)
	}
	// Step 95: t6 = [0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffdf3288fa7113b6d26bb58da4085b309ca37163d548de3]a
	t6.Add(&z, &t6)
	// Step 96: t6 = [0x7ffffffffffffffffffffffffffffffffffffffffffffffffffffffef99447d3889db6935dac6d2042d984e51b8b1eaa46f180]a
	for range 7 {
		t6.Double(&t6)
	}
	// Step 97: t5 = [0x7ffffffffffffffffffffffffffffffffffffffffffffffffffffffef99447d3889db6935dac6d2042d984e51b8b1eaa46f185]a
	t5.Add(&t5, &t6)
	// Step 98: t5 = [0x1fffffffffffffffffffffffffffffffffffffffffffffffffffffffbe6511f4e2276da4d76b1b4810b6613946e2c7aa91bc6140]a
	for range 6 {
		t5.Double(&t5)
	}
	// Step 99: t4 = [0x1fffffffffffffffffffffffffffffffffffffffffffffffffffffffbe6511f4e2276da4d76b1b4810b6613946e2c7aa91bc6149]a
	t4.Add(&t4, &t5)
	// Step 100: t4 = [0x7ffffffffffffffffffffffffffffffffffffffffffffffffffffffef99447d3889db6935dac6d2042d984e51b8b1eaa46f185240]a
	for range 6 {
		t4.Double(&t4)
	}
	// Step 101: t3 = [0x7ffffffffffffffffffffffffffffffffffffffffffffffffffffffef99447d3889db6935dac6d2042d984e51b8b1eaa46f185255]a
	t3.Add(&t3, &t4)
	// Step 102: t3 = [0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffdf3288fa7113b6d26bb58da4085b309ca37163d548de30a4aa0]a
	for range 5 {
		t3.Double(&t3)
	}
	// Step 103: t2 = [0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffdf3288fa7113b6d26bb58da4085b309ca37163d548de30a4aad]a
	t2.Add(&t2, &t3)
	// Step 104: t2 = [0x7ffffffffffffffffffffffffffffffffffffffffffffffffffffffef99447d3889db6935dac6d2042d984e51b8b1eaa46f18525568]a
	for range 3 {
		t2.Double(&t2)
	}
	// Step 105: t2 = [0x7ffffffffffffffffffffffffffffffffffffffffffffffffffffffef99447d3889db6935dac6d2042d984e51b8b1eaa46f1852556b]a
	t2.Add(&z, &t2)
	// Step 106: t2 = [0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffdf3288fa7113b6d26bb58da4085b309ca37163d548de30a4aad600]a
	for range 9 {
		t2.Double(&t2)
	}
	// Step 107: t1 = [0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffdf3288fa7113b6d26bb58da4085b309ca37163d548de30a4aad611]a
	t1.Add(&t1, &t2)
	// Step 108: t1 = [0x3fffffffffffffffffffffffffffffffffffffffffffffffffffffff7cca23e9c44edb49aed63690216cc2728dc58f552378c292ab58440]a
	for range 6 {
		t1.Double(&t1)
	}
	// Step 109: t0 = [0x3fffffffffffffffffffffffffffffffffffffffffffffffffffffff7cca23e9c44edb49aed63690216cc2728dc58f552378c292ab5844f]a
	t0.Add(&t0, &t1)
	// Step 110: t0 = [0x3fffffffffffffffffffffffffffffffffffffffffffffffffffffff7cca23e9c44edb49aed63690216cc2728dc58f552378c292ab5844f0]a
	for range 4 {
		t0.Double(&t0)
	}
	// Step 111: z = [0x3fffffffffffffffffffffffffffffffffffffffffffffffffffffff7cca23e9c44edb49aed63690216cc2728dc58f552378c292ab5844f3]a
	z.Add(&z, &t0)

	return p.Set(&z)
}

// Code generated by internal/generator END mulByOrder. DO NOT EDIT.

// Code generated by internal/generator BEGIN mulByCofactor. DO NOT EDIT.

// mulByCofactor computes p = [h]a and returns p, where
//
//	h = 4
//
// is the cofactor. The result lies in the prime-order subgroup for any
// curve point a. The chain was generated with github.com/mmcloughlin/addchain v0.4.0.
//
// Addition chain:
//
//	return  1 << 2
//
// Operations: 2 point doublings, 0 point additions.
func (p *PointProj) mulByCofactor(a *PointProj) *PointProj {
	// x is a working copy of the input so that p may alias a.
	var x, z PointProj
	x.Set(a)

	// Step 1: z = [0x4]a
	z.Double(&x)
	z.Double(&z)

	return p.Set(&z)
}

// Code generated by internal/generator END mulByCofactor. DO NOT EDIT.
