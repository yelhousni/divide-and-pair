package fp2

import "math/big"

// E2 is an element of Fp2 = Fp[i]/(i² + 1) where p = 2^127 - 1.
// An element is represented as A0 + A1*i.
type E2 struct {
	A0, A1 Element
}

func (z *E2) SetZero() *E2 {
	z.A0.SetZero()
	z.A1.SetZero()
	return z
}

func (z *E2) SetOne() *E2 {
	z.A0.SetOne()
	z.A1.SetZero()
	return z
}

func (z *E2) Set(x *E2) *E2 {
	z.A0.Set(&x.A0)
	z.A1.Set(&x.A1)
	return z
}

func (z *E2) Equal(x *E2) bool {
	return z.A0.Equal(&x.A0) && z.A1.Equal(&x.A1)
}

func (z *E2) IsZero() bool {
	return z.A0.IsZero() && z.A1.IsZero()
}

func (z *E2) IsOne() bool {
	return z.A0.IsOne() && z.A1.IsZero()
}

// Neg sets z = -x and returns z.
func (z *E2) Neg(x *E2) *E2 {
	z.A0.Neg(&x.A0)
	z.A1.Neg(&x.A1)
	return z
}

// Add sets z = x + y and returns z.
func (z *E2) Add(x, y *E2) *E2 {
	z.A0.Add(&x.A0, &y.A0)
	z.A1.Add(&x.A1, &y.A1)
	return z
}

// Sub sets z = x - y and returns z.
func (z *E2) Sub(x, y *E2) *E2 {
	z.A0.Sub(&x.A0, &y.A0)
	z.A1.Sub(&x.A1, &y.A1)
	return z
}

// Double sets z = 2*x and returns z.
func (z *E2) Double(x *E2) *E2 {
	z.A0.Double(&x.A0)
	z.A1.Double(&x.A1)
	return z
}

// Mul sets z = x * y (Karatsuba) and returns z.
// (a+bi)(c+di) = (ac - bd) + (ad + bc)i
// Karatsuba: (ac-bd) + ((a+b)(c+d) - ac - bd)i
func (z *E2) Mul(x, y *E2) *E2 {
	var ac, bd, apb, cpd Element
	ac.Mul(&x.A0, &y.A0)
	bd.Mul(&x.A1, &y.A1)
	apb.Add(&x.A0, &x.A1)
	cpd.Add(&y.A0, &y.A1)

	z.A1.Mul(&apb, &cpd)
	z.A1.Sub(&z.A1, &ac)
	z.A1.Sub(&z.A1, &bd)
	z.A0.Sub(&ac, &bd)
	return z
}

// Square sets z = x² and returns z.
// (a+bi)² = (a²-b²) + 2ab*i = (a+b)(a-b) + 2ab*i
func (z *E2) Square(x *E2) *E2 {
	var apb, amb, twoab Element
	apb.Add(&x.A0, &x.A1)
	amb.Sub(&x.A0, &x.A1)
	twoab.Mul(&x.A0, &x.A1)
	twoab.Double(&twoab)

	z.A0.Mul(&apb, &amb)
	z.A1.Set(&twoab)
	return z
}

// MulByElement sets z = x * a (where a is in Fp) and returns z.
func (z *E2) MulByElement(x *E2, a *Element) *E2 {
	z.A0.Mul(&x.A0, a)
	z.A1.Mul(&x.A1, a)
	return z
}

// Conjugate sets z = x* (complex conjugate) and returns z.
func (z *E2) Conjugate(x *E2) *E2 {
	z.A0.Set(&x.A0)
	z.A1.Neg(&x.A1)
	return z
}

// Norm returns the norm a² + b² in Fp.
func (z *E2) Norm() Element {
	var a2, b2 Element
	a2.Square(&z.A0)
	b2.Square(&z.A1)
	a2.Add(&a2, &b2)
	return a2
}

// Inverse sets z = 1/x and returns z.
// 1/(a+bi) = (a-bi)/(a²+b²)
func (z *E2) Inverse(x *E2) *E2 {
	norm := x.Norm()
	norm.Inverse(&norm)
	z.A0.Mul(&x.A0, &norm)
	var negB Element
	negB.Neg(&x.A1)
	z.A1.Mul(&negB, &norm)
	return z
}

// Exp sets z = x^k (binary method) and returns z.
func (z *E2) Exp(x *E2, k *big.Int) *E2 {
	if k.Sign() == 0 {
		return z.SetOne()
	}
	e := new(big.Int).Set(k)
	if e.Sign() < 0 {
		e.Neg(e)
		x = new(E2).Inverse(x)
	}

	z.Set(x)
	for i := e.BitLen() - 2; i >= 0; i-- {
		z.Square(z)
		if e.Bit(i) == 1 {
			z.Mul(z, x)
		}
	}
	return z
}

// Legendre returns the Legendre symbol of the norm of z.
// Returns 1 if z is a non-zero square in Fp2, -1 if not, 0 if z is zero.
func (z *E2) Legendre() int {
	norm := z.Norm()
	return norm.Legendre()
}

// SetRandom sets z to a random Fp2 element.
func (z *E2) SetRandom() (*E2, error) {
	if _, err := z.A0.SetRandom(); err != nil {
		return nil, err
	}
	if _, err := z.A1.SetRandom(); err != nil {
		return nil, err
	}
	return z, nil
}

// CyclotomicSquare sets z = x² assuming x is on the norm-1 torus (x·x̄ = 1).
// For x = a+bi with a²+b² = 1: x² = (2a²-1) + ((a+b)²-1)i.
// Cost: 2S instead of 2M.
func (z *E2) CyclotomicSquare(x *E2) *E2 {
	var a2, apb2 Element
	a2.Square(&x.A0)        // a²
	apb2.Add(&x.A0, &x.A1)  // a+b
	apb2.Square(&apb2)      // (a+b)²

	z.A0.Double(&a2)         // 2a²
	z.A0.Sub(&z.A0, &one)   // 2a² - 1

	z.A1.Sub(&apb2, &one)   // (a+b)² - 1 = 2ab
	return z
}

var one Element

func init() {
	one.SetOne()
}

// String returns the string representation.
func (z *E2) String() string {
	return z.A0.String() + " + " + z.A1.String() + "*i"
}
