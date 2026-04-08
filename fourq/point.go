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
