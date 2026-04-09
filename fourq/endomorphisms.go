package fourq

import (
	"math/big"
	"sync"

	fp2 "github.com/yelhousni/divide-and-pair/fourq/fp2"
)

// This file implements the two efficiently computable endomorphisms on FourQ
// and an endomorphism-based subgroup membership test, following Costello-Longa
// "FourQ: four-dimensional decompositions on a Q-curve over the Mersenne prime".
//
// The endomorphisms are:
//   ψ (psi): Q-curve endomorphism of degree 2p
//   φ (phi): CM endomorphism of degree 5p (discriminant D = -40)
//
// Both are computed as ψ = τ̂ ∘ innerψ ∘ τ and φ = τ̂ ∘ innerφ ∘ τ, where
// τ: E → Ê and τ̂: Ê → E are dual 4-isogenies (Proposition 1), and the
// inner maps operate on the isogenous curve Ê.

var (
	endoInitOnce sync.Once

	// √d̂ for τ and τ̂ isogenies, where d̂ = -1/(1+d)
	sqrtDhat fp2.E2

	// Inner ψ constants (Section 3.1)
	cPsi1, cPsi2, cPsi3, cPsi4 fp2.E2

	// Inner φ constants (Section 3.2)
	cPhi1, cPhi2, cPhi3, cPhi4, cPhi5  fp2.E2
	cPhi6, cPhi7, cPhi8, cPhi9, cPhi10 fp2.E2

	// Eigenvalues on E(Fp2)[N]
	lambdaPsi    big.Int
	lambdaPhi    big.Int
	lambdaPhiPsi big.Int

	// Short lattice vector for SMT (from LLL, ~61 bits)
	endoC1, endoC2, endoC3, endoC4 big.Int
)

func initEndoConstants() {
	initOnce.Do(initCurveParams)

	// All constants verified in Sage (sage/fourq_endo_constants.sage)
	sqrtDhat.A0.SetString("120525532476903946900736407295642634213")
	sqrtDhat.A1.SetString("110680464442257309687")

	cPsi1.A0.SetString("55340232221128654846")
	cPsi1.A1.SetString("82748376373132255413681934895176861263")
	cPsi2.A0.SetString("170141183460469231621006839273626796022")
	cPsi2.A1.SetString("68193373840384168448244632122363004318")
	cPsi3.A0.SetString("1826227663297245609844")
	cPsi3.A1.SetString("0")
	cPsi4.A0.SetString("1051464412201444442036")
	cPsi4.A1.SetString("47353979519174373375585696079102794256")

	cPhi1.A0.SetString("170141183460469231621006839273626796040")
	cPhi1.A1.SetString("120525532476903946900736407295642634213")
	cPhi2.A0.SetString("92233720368547758087")
	cPhi2.A1.SetString("131306912742858181648727312260439119609")
	cPhi3.A0.SetString("276701161105643274261")
	cPhi3.A1.SetString("160666015865631300014011952927357137809")
	cPhi4.A0.SetString("36893488147419103235")
	cPhi4.A1.SetString("107027644557995218531204623577807990436")
	cPhi5.A0.SetString("55340232221128654851")
	cPhi5.A1.SetString("24279268184862963117522688682631129173")
	cPhi6.A0.SetString("184467440737095516175")
	cPhi6.A1.SetString("92472642025247131565767320804994133491")
	cPhi7.A0.SetString("1660206966633859645560")
	cPhi7.A1.SetString("74020502950125156999236689470520806275")
	cPhi8.A0.SetString("2213609288845146194095")
	cPhi8.A1.SetString("41136875617835130850915686116875720836")
	cPhi9.A0.SetString("3135946492530623774960")
	cPhi9.A1.SetString("41635071732389019719735756359456329456")
	cPhi10.A0.SetString("39844967199212631493615")
	cPhi10.A1.SetString("21045324596686230484035983431638590725")
}

func init() {
	lambdaPsi.SetString("43760231755807040276284855770911078252536368422635318376310714077319867016", 10)
	lambdaPhi.SetString("12098939722099758392970036154455447385486035337534694534042314319425271908", 10)
	lambdaPhiPsi.SetString("42306631464858389077121485826103052927889336227980061597005228525569104570", 10)

	endoC1.SetInt64(-650487742939046294)
	endoC2.SetInt64(1397215820276968864)
	endoC3.SetInt64(-523086274270593807)
	endoC4.SetInt64(598824378691085905)
}

// multiScalarMulIsZero checks if [s1]P1 + [s2]P2 + [s3]P3 + [s4]P4 == O
// using Shamir's trick (joint double-and-add) in projective coordinates.
// This processes all 4 scalars in a single pass of max(bitlen(si)) doublings,
// precomputing all 2^4-1=15 possible sums of the base points.
func multiScalarMulIsZero(P1, P2, P3, P4 *PointAffine, s1, s2, s3, s4 *big.Int) bool {
	// Convert to projective and make all scalars positive
	var pts [4]PointProj
	var scalars [4]big.Int
	inputs := [4]*PointAffine{P1, P2, P3, P4}
	inputScalars := [4]*big.Int{s1, s2, s3, s4}
	for i := range inputs {
		pts[i].FromAffine(inputs[i])
		scalars[i].Set(inputScalars[i])
		if scalars[i].Sign() < 0 {
			pts[i].Neg(&pts[i])
			scalars[i].Neg(&scalars[i])
		}
	}

	// Precompute all 15 non-trivial sums in projective
	var table [15]PointProj
	table[0].Set(&pts[0])             // 0001
	table[1].Set(&pts[1])             // 0010
	table[2].Add(&pts[0], &pts[1])    // 0011
	table[3].Set(&pts[2])             // 0100
	table[4].Add(&pts[0], &pts[2])    // 0101
	table[5].Add(&pts[1], &pts[2])    // 0110
	table[6].Add(&table[2], &pts[2])  // 0111
	table[7].Set(&pts[3])             // 1000
	table[8].Add(&pts[0], &pts[3])    // 1001
	table[9].Add(&pts[1], &pts[3])    // 1010
	table[10].Add(&table[2], &pts[3]) // 1011
	table[11].Add(&pts[2], &pts[3])   // 1100
	table[12].Add(&table[4], &pts[3]) // 1101
	table[13].Add(&table[5], &pts[3]) // 1110
	table[14].Add(&table[6], &pts[3]) // 1111

	// Find max bit length
	maxBits := 0
	for i := range scalars {
		if b := scalars[i].BitLen(); b > maxBits {
			maxBits = b
		}
	}

	// Joint double-and-add (MSB to LSB) in projective
	var result PointProj
	result.setInfinity()

	for i := maxBits - 1; i >= 0; i-- {
		result.Double(&result)

		mask := 0
		for j := range scalars {
			if scalars[j].Bit(i) == 1 {
				mask |= 1 << j
			}
		}
		if mask != 0 {
			result.Add(&result, &table[mask-1])
		}
	}

	return result.IsZero()
}

// mulByI sets z = z * i where i² = -1.
// (a + b·i) * i = -b + a·i
func mulByI(z *fp2.E2) {
	z.A0, z.A1 = z.A1, z.A0
	z.A0.Neg(&z.A0)
}

// tau computes the 4-isogeny τ: E → Ê (Proposition 1).
//
//	τ(x,y) = ( 2xy / ((x²+y²)·√d̂), (x²−y²+2) / (y²−x²) )
func tau(p *PointAffine, rx, ry *fp2.E2) {
	var x2, y2, xy, x2py2, x2my2 fp2.E2

	x2.Square(&p.X)
	y2.Square(&p.Y)
	xy.Mul(&p.X, &p.Y)
	x2py2.Add(&x2, &y2)
	x2my2.Sub(&x2, &y2)

	// rx = 2xy / ((x²+y²)·√d̂)
	var num, den fp2.E2
	num.Double(&xy)
	den.Mul(&x2py2, &sqrtDhat)
	den.Inverse(&den)
	rx.Mul(&num, &den)

	// ry = (x²−y²+2) / (y²−x²)
	var two fp2.E2
	two.A0.SetUint64(2)
	num.Add(&x2my2, &two)
	den.Neg(&x2my2) // y²−x²
	den.Inverse(&den)
	ry.Mul(&num, &den)
}

// tauDual computes the dual 4-isogeny τ̂: Ê → E (Proposition 1).
//
//	τ̂(x,y) = ( 2xy·√d̂ / (x²−y²+2), (y²−x²) / (y²+x²) )
func tauDual(x, y *fp2.E2, result *PointAffine) {
	var x2, y2, xy, x2py2, x2my2 fp2.E2

	x2.Square(x)
	y2.Square(y)
	xy.Mul(x, y)
	x2py2.Add(&x2, &y2)
	x2my2.Sub(&x2, &y2)

	// result.X = 2xy·√d̂ / (x²−y²+2)
	var num, den, two fp2.E2
	num.Double(&xy)
	num.Mul(&num, &sqrtDhat)
	two.A0.SetUint64(2)
	den.Add(&x2my2, &two)
	den.Inverse(&den)
	result.X.Mul(&num, &den)

	// result.Y = (y²−x²) / (y²+x²)
	num.Neg(&x2my2) // y²−x²
	x2py2.Inverse(&x2py2)
	result.Y.Mul(&num, &x2py2)
}

// innerPsi computes (δψ_Wδ⁻¹): Ê → Ê (Section 3.1).
//
//	X' = 2i·x^p·cPsi1 / (y^p · ((x^p)²·cPsi3 + cPsi4))
//	Y' = (cPsi2 − (x^p)²) / (cPsi2 + (x^p)²)
//
// where x^p = conj(x), y^p = conj(y).
func innerPsi(x, y, rx, ry *fp2.E2) {
	var xp, yp, xp2 fp2.E2
	xp.Conjugate(x)
	yp.Conjugate(y)
	xp2.Square(&xp)

	// rx = 2i·xp·cPsi1 / (yp · (xp²·cPsi3 + cPsi4))
	var num, den fp2.E2
	num.Mul(&xp, &cPsi1)
	num.Double(&num)
	mulByI(&num)
	den.Mul(&xp2, &cPsi3)
	den.Add(&den, &cPsi4)
	den.Mul(&den, &yp)
	den.Inverse(&den)
	rx.Mul(&num, &den)

	// ry = (cPsi2 − xp²) / (cPsi2 + xp²)
	var numY, denY fp2.E2
	numY.Sub(&cPsi2, &xp2)
	denY.Add(&cPsi2, &xp2)
	denY.Inverse(&denY)
	ry.Mul(&numY, &denY)
}

// innerPhi computes (δφ_Wδ⁻¹): Ê → Ê (Section 3.2).
//
//	x_φ = (cPhi1·x·(y²−cPhi2·y+cPhi3)·(y²+cPhi2·y+cPhi3))^p
//	      / ((y²+cPhi4·y+cPhi5)·(y²−cPhi4·y+cPhi5))^p
//	y_φ = (cPhi6·(5y⁴+cPhi7·y²+cPhi8))^p
//	      / (5y·(y⁴+cPhi9·y²+cPhi10))^p
func innerPhi(x, y, rx, ry *fp2.E2) {
	var y2, y4, tmp fp2.E2
	y2.Square(y)
	y4.Square(&y2)

	// x_φ numerator: cPhi1·x·(y²−cPhi2·y+cPhi3)·(y²+cPhi2·y+cPhi3)
	var t1, t2, numX fp2.E2
	// t1 = y² − cPhi2·y + cPhi3
	tmp.Mul(&cPhi2, y)
	t1.Sub(&y2, &tmp)
	t1.Add(&t1, &cPhi3)
	// t2 = y² + cPhi2·y + cPhi3
	t2.Add(&y2, &tmp)
	t2.Add(&t2, &cPhi3)
	numX.Mul(&cPhi1, x)
	numX.Mul(&numX, &t1)
	numX.Mul(&numX, &t2)

	// x_φ denominator: (y²+cPhi4·y+cPhi5)·(y²−cPhi4·y+cPhi5)
	var t3, t4, denX fp2.E2
	tmp.Mul(&cPhi4, y)
	t3.Add(&y2, &tmp)
	t3.Add(&t3, &cPhi5)
	t4.Sub(&y2, &tmp)
	t4.Add(&t4, &cPhi5)
	denX.Mul(&t3, &t4)

	// Apply Frobenius (conjugation) then divide
	numX.Conjugate(&numX)
	denX.Conjugate(&denX)
	denX.Inverse(&denX)
	rx.Mul(&numX, &denX)

	// y_φ numerator: cPhi6·(5·y⁴+cPhi7·y²+cPhi8)
	var numY fp2.E2
	var five fp2.E2
	five.A0.SetUint64(5)
	numY.Mul(&five, &y4)
	tmp.Mul(&cPhi7, &y2)
	numY.Add(&numY, &tmp)
	numY.Add(&numY, &cPhi8)
	numY.Mul(&cPhi6, &numY)

	// y_φ denominator: 5·y·(y⁴+cPhi9·y²+cPhi10)
	var denY fp2.E2
	denY.Set(&y4)
	tmp.Mul(&cPhi9, &y2)
	denY.Add(&denY, &tmp)
	denY.Add(&denY, &cPhi10)
	denY.Mul(&denY, y)
	denY.Mul(&denY, &five)

	// Apply Frobenius then divide
	numY.Conjugate(&numY)
	denY.Conjugate(&denY)
	denY.Inverse(&denY)
	ry.Mul(&numY, &denY)
}

// Psi computes the Q-curve endomorphism ψ: E → E.
// ψ = τ̂ ∘ innerψ ∘ τ. Cost: 19M + 8S + 14.5A in Fp2.
func (p *PointAffine) Psi(q *PointAffine) *PointAffine {
	endoInitOnce.Do(initEndoConstants)

	if q.IsZero() {
		return p.SetInfinity()
	}

	var tx, ty, mx, my fp2.E2
	tau(q, &tx, &ty)
	innerPsi(&tx, &ty, &mx, &my)
	tauDual(&mx, &my, p)
	return p
}

// Phi computes the CM endomorphism φ: E → E.
// φ = τ̂ ∘ innerφ ∘ τ. Cost: 30M + 11S + 20.5A in Fp2.
func (p *PointAffine) Phi(q *PointAffine) *PointAffine {
	endoInitOnce.Do(initEndoConstants)

	if q.IsZero() {
		return p.SetInfinity()
	}

	var tx, ty, mx, my fp2.E2
	tau(q, &tx, &ty)
	innerPhi(&tx, &ty, &mx, &my)
	tauDual(&mx, &my, p)
	return p
}

// isInSubGroupEndo tests subgroup membership using the endomorphism-based
// GLV decomposition. It checks:
//
//	[c1]P + [c2]φ(P) + [c3]ψ(P) + [c4]ψ(φ(P)) == O
//
// where (c1,c2,c3,c4) = (-650487742939046294, 1397215820276968864,
// -523086274270593807, 598824378691085905) is a short vector (~61 bits)
// in the lattice of zero decompositions.
//
// Cost: endomorphisms (68M+27S+49.5A) + 4-dim MSM with ~61-bit scalars.
func (p *PointAffine) isInSubGroupEndo() bool {
	endoInitOnce.Do(initEndoConstants)
	initOnce.Do(initCurveParams)

	if p.IsZero() {
		return true
	}

	// The endomorphisms ψ, φ are composed via the 4-isogeny τ: E → Ê.
	// Points in ker(τ) (i.e., low-order torsion killed by τ) are mapped
	// to the identity, making the test degenerate. We guard against this
	// by checking that [392]P (cofactor clearing) is not identity — if it
	// is, P has order dividing 392, which means P ∉ E[N] (since N is prime
	// and gcd(N, 392) = 1).
	var cofP PointAffine
	cofP.ScalarMultiplication(p, big.NewInt(392))
	if cofP.IsZero() {
		return false
	}

	// Compute endomorphism images using proper formulas
	var phiP, psiP, psiPhiP PointAffine
	phiP.Phi(p)
	psiP.Psi(p)
	psiPhiP.Psi(&phiP) // ψ(φ(P))

	// 4-dimensional multi-scalar multiplication using Shamir's trick:
	// Process all 4 scalars simultaneously in a single pass of ~61 doublings.
	return multiScalarMulIsZero(p, &phiP, &psiP, &psiPhiP, &endoC1, &endoC2, &endoC3, &endoC4)
}
