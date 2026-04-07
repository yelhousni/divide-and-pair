package curve448

import (
	"math/big"
	"sync"

	fp "github.com/yelhousni/divide-and-pair/curve448/fp"
)

// curveData holds curve parameters: Ax^2 + y^2 = 1 + Dx^2*y^2.
type curveData struct {
	A, D     fp.Element
	Cofactor fp.Element
	Order    big.Int
	Base     PointAffine
}

func curveParameters() curveData {
	initOnce.Do(initCurveParams)
	var res curveData
	res.A.Set(&curveParams.A)
	res.D.Set(&curveParams.D)
	res.Cofactor.Set(&curveParams.Cofactor)
	res.Order.Set(&curveParams.Order)
	res.Base.Set(&curveParams.Base)
	return res
}

// Generators returns the generator of the prime-order subgroup.
func Generators() (base PointAffine) {
	initOnce.Do(initCurveParams)
	base.Set(&curveParams.Base)
	return
}

// CurveCoefficients returns the a and d coefficients of
// A*x² + y² = 1 + D*x²*y².
func CurveCoefficients() (a, d fp.Element) {
	initOnce.Do(initCurveParams)
	a.Set(&curveParams.A)
	d.Set(&curveParams.D)
	return
}

// Order returns the prime subgroup order.
func Order() *big.Int {
	initOnce.Do(initCurveParams)
	return new(big.Int).Set(&curveParams.Order)
}

var (
	initOnce    sync.Once
	curveParams curveData
)

func initCurveParams() {
	// Curve448 Edwards: x^2 + y^2 = 1 + d*x^2*y^2
	// a = 1
	// d = -39081
	// cofactor = 4
	// ell = 2^446 - 13818066809895115352007386748515426880336692474882178609894547503885
	curveParams.A.SetOne() // A = 1

	// d = -39081
	curveParams.D.SetUint64(39081)
	curveParams.D.Neg(&curveParams.D)

	curveParams.Cofactor.SetUint64(4)
	curveParams.Order.SetString("181709681073901722637330951972001133588410340171829515070372549795146003961539585716195755291692375963310293709091662304773755859649779", 10)

	// Base point from RFC 8032 / crrl (Pornin's implementation)
	// This base point has prime order ell (it is already in the prime-order subgroup)
	curveParams.Base.X.SetString("224580040295924300187604334099896036246789641632564134246125461686950415467406032909029192869357953282578032075146446173674602635247710")
	curveParams.Base.Y.SetString("298819210078481492676017930443930673437544040154080242095928241372331506189835876003536878655418784733982303233503462500531545062832660")
}

// mulByA multiplies fp.Element by curveParams.A (= 1), i.e., identity
func mulByA(x *fp.Element) {
	// a = 1, nothing to do
	_ = x
}
