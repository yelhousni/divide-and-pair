package jubjub

import (
	"math/big"
	"sync"

	fp "github.com/yelhousni/divide-and-pair/jubjub/fp"
)

// CurveParams curve parameters: Ax^2 + y^2 = 1 + Dx^2*y^2
type CurveParams struct {
	A, D     fp.Element
	Cofactor fp.Element
	Order    big.Int
	Base     PointAffine
}

// GetEdwardsCurve returns the JubJub twisted Edwards curve on BLS12-381/Fr
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
	// JubJub twisted Edwards: -x^2 + y^2 = 1 + d*x^2*y^2
	// a = -1
	// d = 19257038036680949359750312669786877991949435402254120286184196891950884077233
	// cofactor = 8
	// ell = 6554484396890773809930967563523245729705921265872317281365359162392183254199
	curveParams.A.SetOne()
	curveParams.A.Neg(&curveParams.A) // A = -1

	curveParams.D.SetString("19257038036680949359750312669786877991949435402254120286184196891950884077233")

	curveParams.Cofactor.SetUint64(8)
	curveParams.Order.SetString("6554484396890773809930967563523245729705921265872317281365359162392183254199", 10)

	// Base point from gnark-crypto (generator of prime-order subgroup)
	curveParams.Base.X.SetString("23426137002068529236790192115758361610982344002369094106619281483467893291614")
	curveParams.Base.Y.SetString("39325435222430376843701388596190331198052476467368316772266670064146548432123")
}

// mulByA multiplies fp.Element by curveParams.A (= -1)
func mulByA(x *fp.Element) {
	x.Neg(x)
}
