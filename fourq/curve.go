package fourq

import (
	"math/big"
	"sync"

	fp2 "github.com/yelhousni/divide-and-pair/fourq/fp2"
)

// curveData holds FourQ twisted Edwards curve parameters over Fp2.
// Curve equation: -x² + y² = 1 + d·x²·y² (a = -1)
type curveData struct {
	A, D  fp2.E2
	Order big.Int // prime subgroup order N
	Base  PointAffine
}

func curveParameters() curveData {
	initOnce.Do(initCurveParams)
	var res curveData
	res.A.Set(&curveParams.A)
	res.D.Set(&curveParams.D)
	res.Order.Set(&curveParams.Order)
	res.Base.Set(&curveParams.Base)
	return res
}

// Generators returns the generator of the prime-order subgroup.
func Generators() (g1Aff PointAffine) {
	initOnce.Do(initCurveParams)
	g1Aff.Set(&curveParams.Base)
	return
}

// CurveCoefficients returns the a and d coefficients of
// -x² + y² = 1 + d·x²·y².
func CurveCoefficients() (a, d fp2.E2) {
	initOnce.Do(initCurveParams)
	a.Set(&curveParams.A)
	d.Set(&curveParams.D)
	return
}

// Order returns the prime subgroup order N.
func Order() *big.Int {
	initOnce.Do(initCurveParams)
	return new(big.Int).Set(&curveParams.Order)
}

var (
	initOnce    sync.Once
	curveParams curveData
)

func initCurveParams() {
	// a = -1 (in Fp2: (-1, 0))
	curveParams.A.A0.SetOne()
	curveParams.A.A0.Neg(&curveParams.A.A0)
	curveParams.A.A1.SetZero()

	// d = d0 + d1*i from CIRCL/FourQ paper
	curveParams.D.A0.SetString("4205857648805777768770")
	curveParams.D.A1.SetString("125317048443780598345676279555970305165")

	// Prime subgroup order N ≈ 2^246
	curveParams.Order.SetString("73846995687063900142583536357581573884798075859800097461294096333596429543", 10)

	// Generator point from CIRCL
	curveParams.Base.X.A0.SetString("34832242333165934151976439273177494442")
	curveParams.Base.X.A1.SetString("40039530084877881816286215037915002870")
	curveParams.Base.Y.A0.SetString("18941146186793715734774048165794132615")
	curveParams.Base.Y.A1.SetString("146361984425930646555497992424795179868")
}

// mulByA multiplies an Fp2 element by a = -1.
func mulByA(x *fp2.E2) {
	x.Neg(x)
}
