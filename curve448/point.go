package curve448

import (
	"math/big"

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
	var res PointAffine
	res.X.SetZero()
	res.Y.SetOne() // identity

	s := new(big.Int).Set(scalar)
	if s.Sign() < 0 {
		var neg PointAffine
		neg.Neg(p1)
		p1 = &neg
		s.Neg(s)
	}

	var tmp PointAffine
	tmp.Set(p1)

	for s.Sign() > 0 {
		if s.Bit(0) == 1 {
			res.Add(&res, &tmp)
		}
		tmp.Double(&tmp)
		s.Rsh(s, 1)
	}

	p.Set(&res)
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
