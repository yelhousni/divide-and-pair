package curve448

import (
	"crypto/rand"
	"math/big"
	"testing"

	"github.com/leanovate/gopter"
	"github.com/leanovate/gopter/prop"
	fp "github.com/yelhousni/divide-and-pair/curve448/fp"
)

func GenBigInt() gopter.Gen {
	return func(genParams *gopter.GenParameters) *gopter.GenResult {
		var s big.Int
		var b [fp.Bytes]byte
		_, err := rand.Read(b[:])
		if err != nil {
			panic(err)
		}
		s.SetBytes(b[:])
		genResult := gopter.NewGenResult(s, gopter.NoShrinker)
		return genResult
	}
}

func TestReceiverIsOperand(t *testing.T) {
	t.Parallel()
	parameters := gopter.DefaultTestParameters()
	parameters.MinSuccessfulTests = 10

	properties := gopter.NewProperties(parameters)

	properties.Property("Add affine: receiver as operand", prop.ForAll(
		func() bool {
			params := curveParameters()
			var p1, p2, p3 PointAffine
			p1.Set(&params.Base)
			p2.Set(&params.Base)
			p3.Set(&params.Base)

			p3.Add(&p1, &p2)
			p1.Add(&p1, &p2)
			if !p3.Equal(&p1) {
				return false
			}
			p1.Set(&params.Base)
			p2.Add(&p1, &p2)
			return p2.Equal(&p3)
		},
	))

	properties.Property("Double affine: receiver as operand", prop.ForAll(
		func() bool {
			params := curveParameters()
			var p1, p2 PointAffine
			p1.Set(&params.Base)
			p2.Set(&params.Base)
			p2.Double(&p1)
			p1.Double(&p1)
			return p2.Equal(&p1)
		},
	))

	properties.Property("Neg affine: receiver as operand", prop.ForAll(
		func() bool {
			params := curveParameters()
			var p1, p2 PointAffine
			p1.Set(&params.Base)
			p2.Neg(&p1)
			p1.Neg(&p1)
			return p2.Equal(&p1)
		},
	))

	properties.Property("ScalarMul affine: receiver as operand", prop.ForAll(
		func() bool {
			params := curveParameters()
			var p1, p2 PointAffine
			p1.Set(&params.Base)
			p2.Set(&params.Base)
			s := big.NewInt(10)
			p2.ScalarMultiplication(&p1, s)
			p1.ScalarMultiplication(&p1, s)
			return p2.Equal(&p1)
		},
	))

	properties.TestingRun(t, gopter.ConsoleReporter(false))
}

func TestField(t *testing.T) {
	t.Parallel()
	parameters := gopter.DefaultTestParameters()
	parameters.MinSuccessfulTests = 10

	properties := gopter.NewProperties(parameters)
	genS := GenBigInt()

	properties.Property("MulByA(x) should match Mul(x, curve.A)", prop.ForAll(
		func(s big.Int) bool {
			params := curveParameters()
			var z1, z2 fp.Element
			z1.SetBigInt(&s)
			z2.Mul(&z1, &params.A)
			mulByA(&z1)
			return z1.Equal(&z2)
		},
		genS,
	))

	properties.TestingRun(t, gopter.ConsoleReporter(false))
}

func TestOps(t *testing.T) {
	parameters := gopter.DefaultTestParameters()
	parameters.MinSuccessfulTests = 10

	properties := gopter.NewProperties(parameters)
	genS1 := GenBigInt()
	genS2 := GenBigInt()

	properties.Property("P+0=P", prop.ForAll(
		func(s1 big.Int) bool {
			params := curveParameters()
			var p1, p2, zero PointAffine
			p1.ScalarMultiplication(&params.Base, &s1)
			zero.X.SetZero()
			zero.Y.SetOne()
			p2.Add(&p1, &zero)
			return p2.IsOnCurve() && p2.Equal(&p1)
		},
		genS1,
	))

	properties.Property("P+(-P)=O", prop.ForAll(
		func(s1 big.Int) bool {
			params := curveParameters()
			var p1, p2 PointAffine
			p1.ScalarMultiplication(&params.Base, &s1)
			p2.Neg(&p1)
			p1.Add(&p1, &p2)
			return p1.IsOnCurve() && p1.IsZero()
		},
		genS1,
	))

	properties.Property("P+P=2*P", prop.ForAll(
		func(s big.Int) bool {
			params := curveParameters()
			var p1, p2 PointAffine
			p1.ScalarMultiplication(&params.Base, &s)
			p2.Set(&p1)
			p1.Add(&p1, &p2)
			p2.Double(&p2)
			return p1.IsOnCurve() && p1.Equal(&p2)
		},
		genS1,
	))

	properties.Property("[a]P+[b]P = [a+b]P", prop.ForAll(
		func(s1, s2 big.Int) bool {
			params := curveParameters()
			var p1, p2, p3 PointAffine
			p1.ScalarMultiplication(&params.Base, &s1)
			p2.ScalarMultiplication(&params.Base, &s2)
			p3.Set(&params.Base)
			p2.Add(&p1, &p2)
			s1.Add(&s1, &s2)
			p3.ScalarMultiplication(&params.Base, &s1)
			return p2.IsOnCurve() && p3.Equal(&p2)
		},
		genS1, genS2,
	))

	properties.Property("[5]P=[2][2]P+P", prop.ForAll(
		func(s1 big.Int) bool {
			params := curveParameters()
			var p1, p2 PointAffine
			p1.ScalarMultiplication(&params.Base, &s1)
			p2.Double(&p1)
			p2.Double(&p2)
			p2.Add(&p2, &p1)
			p1.ScalarMultiplication(&p1, big.NewInt(5))
			return p2.IsOnCurve() && p2.Equal(&p1)
		},
		genS1,
	))

	properties.Property("[0]P = O", prop.ForAll(
		func(s1 big.Int) bool {
			params := curveParameters()
			var p1 PointAffine
			p1.ScalarMultiplication(&params.Base, &s1)
			p1.ScalarMultiplication(&p1, big.NewInt(0))
			return p1.IsZero()
		},
		genS1,
	))

	properties.Property("[1]P = P", prop.ForAll(
		func(s1 big.Int) bool {
			params := curveParameters()
			var p1, p2 PointAffine
			p1.ScalarMultiplication(&params.Base, &s1)
			p2.ScalarMultiplication(&p1, big.NewInt(1))
			return p1.Equal(&p2)
		},
		genS1,
	))

	properties.Property("[-1]P = Neg(P)", prop.ForAll(
		func(s1 big.Int) bool {
			params := curveParameters()
			var p1, p2 PointAffine
			p1.ScalarMultiplication(&params.Base, &s1)
			p2.ScalarMultiplication(&p1, big.NewInt(-1))
			p1.Neg(&p1)
			return p1.Equal(&p2)
		},
		genS1,
	))

	properties.Property("[order]P = O", prop.ForAll(
		func(s1 big.Int) bool {
			params := curveParameters()
			var p1 PointAffine
			p1.ScalarMultiplication(&params.Base, &s1)
			p1.ScalarMultiplication(&p1, &params.Order)
			return p1.IsZero()
		},
		genS1,
	))

	properties.Property("Projective round-trip", prop.ForAll(
		func(s1 big.Int) bool {
			params := curveParameters()
			var p1, p2 PointAffine
			var proj PointProj
			p1.ScalarMultiplication(&params.Base, &s1)
			proj.FromAffine(&p1)
			proj.ToAffine(&p2)
			return p1.Equal(&p2)
		},
		genS1,
	))

	properties.TestingRun(t, gopter.ConsoleReporter(false))
}

func TestMulByOrder(t *testing.T) {
	t.Parallel()

	// The identity maps to the identity.
	var id, res PointProj
	id.setInfinity()
	res.mulByOrder(&id)
	if !res.IsZero() {
		t.Fatal("mulByOrder(O) should be the identity")
	}

	parameters := gopter.DefaultTestParameters()
	parameters.MinSuccessfulTests = 50
	properties := gopter.NewProperties(parameters)
	genS := GenBigInt()

	properties.Property("mulByOrder agrees with ScalarMultiplication by the order", prop.ForAll(
		func(s big.Int) bool {
			params := curveParameters()

			// Subgroup point p = [s]Base and non-subgroup point q = p + (0,-1).
			var offset, p, q PointAffine
			offset.X.SetZero()
			offset.Y.SetOne()
			offset.Y.Neg(&offset.Y)
			p.ScalarMultiplication(&params.Base, &s)
			q.Add(&p, &offset)

			for _, pt := range []*PointAffine{&p, &q} {
				var proj, want, got PointProj
				proj.FromAffine(pt)
				want.ScalarMultiplication(&proj, &params.Order)
				got.mulByOrder(&proj)
				if !projEqualMulByOrder(&want, &got) {
					return false
				}

				// The receiver may alias the argument.
				proj.mulByOrder(&proj)
				if !projEqualMulByOrder(&want, &proj) {
					return false
				}
			}
			return true
		},
		genS,
	))

	properties.TestingRun(t, gopter.ConsoleReporter(false))
}

// projEqualMulByOrder compares two projective points as curve points, i.e. up
// to the projective factor.
func projEqualMulByOrder(a, b *PointProj) bool {
	if a.IsZero() || b.IsZero() {
		return a.IsZero() == b.IsZero()
	}
	var aAff, bAff PointAffine
	a.ToAffine(&aAff)
	b.ToAffine(&bAff)
	return aAff.Equal(&bAff)
}

func TestClearCofactor(t *testing.T) {
	t.Parallel()

	// The identity maps to the identity.
	var id PointAffine
	id.X.SetZero()
	id.Y.SetOne()
	id.ClearCofactor()
	if !id.IsZero() {
		t.Fatal("ClearCofactor(O) should be the identity")
	}

	// The order-2 point (0,-1) is killed by cofactor clearing.
	var two PointAffine
	two.X.SetZero()
	two.Y.SetOne()
	two.Y.Neg(&two.Y)
	two.ClearCofactor()
	if !two.IsZero() {
		t.Fatal("ClearCofactor((0,-1)) should be the identity")
	}

	parameters := gopter.DefaultTestParameters()
	parameters.MinSuccessfulTests = 50
	properties := gopter.NewProperties(parameters)
	genS := GenBigInt()

	properties.Property("ClearCofactor equals [h]P and lands in the subgroup", prop.ForAll(
		func(s big.Int) bool {
			params := curveParameters()

			// Subgroup point p = [s]Base and non-subgroup point q = p + (0,-1).
			var offset, p, q PointAffine
			offset.X.SetZero()
			offset.Y.SetOne()
			offset.Y.Neg(&offset.Y)
			p.ScalarMultiplication(&params.Base, &s)
			q.Add(&p, &offset)

			for _, pt := range []*PointAffine{&p, &q} {
				var got, want PointAffine
				got.Set(pt).ClearCofactor()
				want.ScalarMultiplication(pt, big.NewInt(4))
				if !got.Equal(&want) {
					return false
				}
				if !got.IsOnCurve() || !got.isInSubGroupNaive() {
					return false
				}
			}
			return true
		},
		genS,
	))

	properties.TestingRun(t, gopter.ConsoleReporter(false))
}

func BenchmarkClearCofactor(b *testing.B) {
	params := curveParameters()
	k, _ := rand.Int(rand.Reader, &params.Order)
	var p PointAffine
	p.ScalarMultiplication(&params.Base, k)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		var q PointAffine
		q.Set(&p)
		q.ClearCofactor()
	}
}
