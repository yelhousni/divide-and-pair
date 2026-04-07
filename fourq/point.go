package fourq

import (
	"math/big"

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
	denX.Add(&denX, &dxy) // 1 + d*x1x2*y1y2
	denY.SetOne()
	denY.Sub(&denY, &dxy) // 1 - d*x1x2*y1y2

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

// ScalarMultiplication computes scalar * p1.
func (p *PointAffine) ScalarMultiplication(p1 *PointAffine, scalar *big.Int) *PointAffine {
	var res PointAffine
	res.SetInfinity()

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
