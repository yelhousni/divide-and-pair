package curve25519

import (
	"math/big"
	"sync"

	fp "github.com/yelhousni/divide-and-pair/curve25519/fp"
)

// CurveParams curve parameters: Ax^2 + y^2 = 1 + Dx^2*y^2
type CurveParams struct {
	A, D     fp.Element
	Cofactor fp.Element
	Order    big.Int
	Base     PointAffine
}

// GetEdwardsCurve returns the Curve25519 twisted Edwards curve on Fp
func GetEdwardsCurve() CurveParams {
	initOnce.Do(initCurveParams)
	var res CurveParams
	res.A.Set(&curveParams.A)
	res.D.Set(&curveParams.D)
	res.Cofactor.Set(&curveParams.Cofactor)
	res.Order.Set(&curveParams.Order)
	res.Base.Set(&curveParams.Base)
	return res
}

var (
	initOnce    sync.Once
	curveParams CurveParams
)

func initCurveParams() {
	// Curve25519 twisted Edwards: -x^2 + y^2 = 1 + d*x^2*y^2
	// a = -1
	// d = -121665/121666 mod q
	// cofactor = 8
	// ell = 2^252 + 27742317777372353535851937790883648493
	curveParams.A.SetOne()
	curveParams.A.Neg(&curveParams.A) // A = -1

	// d = -121665 * 121666^{-1} mod q
	var num, den fp.Element
	num.SetUint64(121665)
	num.Neg(&num)
	den.SetUint64(121666)
	den.Inverse(&den)
	curveParams.D.Mul(&num, &den)

	curveParams.Cofactor.SetUint64(8)
	curveParams.Order.SetString("7237005577332262213973186563042994240857116359379907606001950938285454250989", 10)

	// Base point from RFC 8032 (Ed25519), converted to Edwards (x, y)
	// y = 4/5 mod q
	var four, five fp.Element
	four.SetUint64(4)
	five.SetUint64(5)
	five.Inverse(&five)
	curveParams.Base.Y.Mul(&four, &five)

	// x is recovered from y: x^2 = (1 - y^2) / (a - d*y^2) = (y^2 - 1) / (d*y^2 - a)
	// since a = -1: x^2 = (y^2 - 1) / (d*y^2 + 1)
	var y2, num2, den2 fp.Element
	y2.Square(&curveParams.Base.Y)
	num2.Sub(&y2, new(fp.Element).SetOne()) // y^2 - 1
	den2.Mul(&curveParams.D, &y2)
	den2.Add(&den2, new(fp.Element).SetOne()) // d*y^2 + 1
	den2.Inverse(&den2)
	var x2 fp.Element
	x2.Mul(&num2, &den2)
	curveParams.Base.X.Sqrt(&x2)

	// Choose the positive x (even, per RFC 8032 convention)
	// The standard Ed25519 base point has x with bit 0 = 0
	var xBytes [32]byte
	b := curveParams.Base.X.Bytes()
	copy(xBytes[:], b[:])
	if xBytes[31]&1 == 1 {
		curveParams.Base.X.Neg(&curveParams.Base.X)
	}
}

// mulByA multiplies fp.Element by curveParams.A (= -1)
func mulByA(x *fp.Element) {
	x.Neg(x)
}
