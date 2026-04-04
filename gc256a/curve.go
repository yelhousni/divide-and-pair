package gc256a

import (
	"math/big"
	"sync"

	fp "github.com/yelhousni/divide-and-pair/gc256a/fp"
)

// CurveParams curve parameters: Ax^2 + y^2 = 1 + Dx^2*y^2
type CurveParams struct {
	A, D     fp.Element
	Cofactor fp.Element
	Order    big.Int
	Base     PointAffine
}

// GetEdwardsCurve returns the GC256A Edwards curve on Fp
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
	// GC256A (id-tc26-gost-3410-2012-256-paramSetA) Edwards curve:
	//   x^2 + y^2 = 1 + d*x^2*y^2  (a = 1, e = 1 in RFC notation)
	// p = 2^256 - 617
	// d = 0x0605F6B7C183FA81578BC39CFAD518132B9DF62897009AF7E522C32D6DC7BFFB
	// m = 0x01000000000000000000000000000000003F63377F21ED98D70456BD55B0D8319C (group order = 4*q)
	// q = 0x400000000000000000000000000000000FD8CDDFC87B6635C115AF556C360C67 (prime subgroup order)
	// Generator (u, v) from RFC 7836 Appendix A.2:
	//   u = 0x0D
	//   v = 0x60CA1E32AA475B348488C38FAB07649CE7EF8DBE87F22E81F92B2592DBA300E7
	curveParams.A.SetOne()

	curveParams.D.SetString("2724414110474605931834268501164757645998726878473076809432604223414351675387")

	curveParams.Cofactor.SetUint64(4)
	curveParams.Order.SetString("28948022309329048855892746252171976963338560298092253442512153408785530358887", 10)

	// Base point from RFC 7836: (u, v) = (0x0D, 0x60CA1E32...)
	// This point has prime order q.
	curveParams.Base.X.SetUint64(13) // u = 0x0D = 13
	curveParams.Base.Y.SetString("43779144989398987843428779166090436406934195821915183574454224403186176950503")
}

// mulByA multiplies fp.Element by curveParams.A (= 1), i.e., identity
func mulByA(x *fp.Element) {
	_ = x
}
