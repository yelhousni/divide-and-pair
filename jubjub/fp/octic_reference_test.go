package fp

import "math/big"

// Arithmetic in Z[t]/(t^4+1), with t a primitive eighth root of unity.
// Coefficients are in ascending powers of t.
type octicInteger [4]big.Int

// This primary prime has norm p. Under reduction modulo it, t maps to
// 7^((p-1)/8) in Fp. It is a primary associate of gcd(p, t - 7^((p-1)/8)).
var octicPrime = func() (v octicInteger) {
	v[0].SetUint64(6451692399033376397)
	v[1].SetUint64(11914338290826914742)
	v[2].SetUint64(1094064715653685168)
	v[3].SetUint64(11306383988181076382)
	v[3].Neg(&v[3])
	return
}()

// z must be distinct from the inputs.
func (z *octicInteger) mul(a, b *octicInteger, term *big.Int) {
	for i := range z {
		z[i].SetInt64(0)
	}
	for i := range a {
		for j := range b {
			term.Mul(&a[i], &b[j])
			if i+j < 4 {
				z[i+j].Add(&z[i+j], term)
			} else {
				z[i+j-4].Sub(&z[i+j-4], term)
			}
		}
	}
}

// a*sigma_5(a) = A+B*i. Multiplication by sigma_5(a)*(A-B*i)
// implements division by a with rational denominator A²+B².
func (a *octicInteger) gaussianNorm(A, B, h *big.Int) {
	A.Mul(&a[1], &a[3])
	A.Lsh(A, 1)
	A.Add(A, h.Mul(&a[0], &a[0]))
	A.Sub(A, h.Mul(&a[2], &a[2]))
	B.Mul(&a[0], &a[2])
	B.Lsh(B, 1)
	B.Add(B, h.Mul(&a[3], &a[3]))
	B.Sub(B, h.Mul(&a[1], &a[1]))
}

type octicDivisor struct {
	adj        octicInteger
	norm, half big.Int
}

func (d *octicDivisor) set(a *octicInteger, A, B, h *big.Int) {
	d.norm.Mul(A, A)
	d.norm.Add(&d.norm, h.Mul(B, B))
	d.half.Rsh(&d.norm, 1)
	for j := 0; j < 2; j++ {
		d.adj[j].Mul(&a[j], A)
		d.adj[j].Add(&d.adj[j], h.Mul(&a[j+2], B))
		d.adj[j+2].Mul(&a[j+2], A)
		d.adj[j+2].Sub(&d.adj[j+2], h.Mul(&a[j], B))
	}
	d.adj[1].Neg(&d.adj[1])
	d.adj[3].Neg(&d.adj[3])
}

var octicPrimeDivisor = func() (d octicDivisor) {
	var A, B, h big.Int
	octicPrime.gaussianNorm(&A, &B, &h)
	d.set(&octicPrime, &A, &B, &h)
	return
}()

// Round the four power-basis coefficients independently. For an error e
// with coefficients in [-1/2,1/2], N(e) <= (sum e_j²)² <= 1. Equality
// in the latter bound requires four half-integers; at those corners N(e)=1/2.
// Thus the remainder has strictly smaller absolute norm than the divisor.
func (d *octicDivisor) remainder(r, a, b, q *octicInteger, h *big.Int) {
	q.mul(a, &d.adj, h)
	for i := range q {
		q[i].Add(&q[i], &d.half)
		q[i].DivMod(&q[i], &d.norm, h)
	}
	r.mul(q, b, h)
	for i := range r {
		r[i].Sub(&a[i], &r[i])
	}
}

func (a *octicInteger) divRamified() {
	a[0].Add(&a[0], &a[1])
	a[0].Sub(&a[0], &a[2])
	a[0].Add(&a[0], &a[3])
	a[0].Rsh(&a[0], 1)
	for i := 1; i < 4; i++ {
		a[i].Sub(&a[i], &a[i-1])
	}
}

// Write the input as t^k * (1+sqrt(2))^l * primary and leave its primary
// associate in a (Joye, Proposition 1). The input has odd absolute norm.
func (a *octicInteger) makePrimary(s *octicInteger) (k, l uint64) {
	for a[0].Bit(0) == 0 || a[2].Bit(0) != 0 {
		a[0], a[1], a[2], a[3] = a[1], a[2], a[3], a[0]
		a[3].Neg(&a[3])
		k++
	}
	l = uint64(a[1].Bit(0))
	if l != 0 {
		// Multiply by (1+sqrt(2))^-1 = -1 + t - t³.
		s[0].Sub(&a[1], &a[3])
		s[1].Add(&a[0], &a[2])
		s[2].Add(&a[1], &a[3])
		s[3].Sub(&a[2], &a[0])
		for i := range a {
			a[i].Sub(&s[i], &a[i])
		}
	}
	// The real coefficient is odd and the other three are even.
	sum := uint(1)
	for i := range a {
		sum += 2 * a[i].Bit(1)
	}
	if sum&2 != 0 {
		for i := range a {
			a[i].Neg(&a[i])
		}
		k += 4
	}
	return
}

func (a *octicInteger) low() (v [4]uint64) {
	for i := range v {
		if words := a[i].Bits(); len(words) != 0 {
			v[i] = uint64(words[0])
		}
		if a[i].Sign() < 0 {
			v[i] = -v[i]
		}
	}
	return
}

func (x *Element) octicSymbolReference() int {
	if x.IsZero() {
		return -1
	}
	var slots [3]octicInteger
	a, b, r := &slots[0], &slots[1], &slots[2]
	var q octicInteger
	var A, B, h big.Int
	var divisor octicDivisor
	// The prime denominator and nonzero field input are coprime. Every
	// swap keeps them coprime, so division by a nonunit cannot yield zero.
	x.BigInt(&r[0])
	for i := range b {
		b[i].Set(&octicPrime[i])
	}
	octicPrimeDivisor.remainder(a, r, b, &q, &h)
	den := octicLawsOf(b.low())
	var j uint64
	for {
		if a[0].Bit(0)^a[1].Bit(0)^a[2].Bit(0)^a[3].Bit(0) == 0 {
			a.divRamified()
			j += den.I
			continue
		}
		k, l := a.makePrimary(&q)
		v := a.low()
		j += k*den.K + l*den.L
		a.gaussianNorm(&A, &B, &h)
		// For primary a, A is 1 mod 4, so A²+B²=1 iff (A,B)=(1,0).
		if B.Sign() == 0 && A.Cmp(big.NewInt(1)) == 0 {
			// Algorithm 1 recovers the primary unit's two residual exponents
			// modulo 8 from its constant and cubic coefficients.
			return int((j + (v[0]-1)*den.K + v[3]*den.L) & 7)
		}
		num := octicLawsOf(v)
		j += num.reciprocity(den)
		divisor.set(a, &A, &B, &h)
		divisor.remainder(r, b, a, &q, &h)
		a, b, r = r, a, b
		den = num
	}
}
