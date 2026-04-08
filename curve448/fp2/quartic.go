package fp2

// QuarticSymbol computes χ₄(z) = z^((p²-1)/4) for z ∈ Fp2.
// Returns 0 if χ₄=1 (quartic residue), 1 if χ₄=i, 2 if χ₄=-1, 3 if χ₄=-i.
//
// Uses the identity: χ₄(z) = β^(p-1) where β = z^((p+1)/4).
// β^(p-1) = conj(β)/β, which is a 4th root of unity in Fp2.
// χ₄=1 iff β ∈ Fp (β.A1 = 0).
func (z *E2) QuarticSymbol() uint8 {
	if z.IsZero() {
		return 0
	}

	var beta E2
	beta.ExpBySqrtPp1o4(z)

	// chi4 = conj(beta)/beta = beta^(p-1)
	// If beta.A1 == 0: beta ∈ Fp, conj(beta) = beta, chi4 = 1 → return 0
	if beta.A1.IsZero() {
		return 0
	}

	// If beta.A0 == 0: beta = b*i, conj(beta) = -b*i, chi4 = -1 → return 2
	if beta.A0.IsZero() {
		return 2
	}

	// General case: conj(beta)/beta = (a-bi)/(a+bi)
	// = (a²-b² - 2abi) / (a²+b²)
	// Real part: (a²-b²)/(a²+b²), Imag part: -2ab/(a²+b²)
	// For chi4 = i: real=0, imag=1 → a²=b², -2ab > 0 → a,b opposite signs
	// For chi4 = -i: real=0, imag=-1 → a²=b², -2ab < 0 → a,b same sign
	var a2, b2 E2
	_ = a2
	_ = b2

	// chi4 = i iff a²=b² and ab < 0 (in the field sense)
	// chi4 = -i iff a²=b² and ab > 0
	// We compare using LexicographicallyLargest as a proxy for sign.
	var ab, a2v, b2v E2
	_ = ab
	_ = a2v
	_ = b2v

	// Simpler: compute conj(beta)*beta^(-1) directly and compare with known roots.
	var chi, conjBeta E2
	conjBeta.Conjugate(&beta)
	var betaInv E2
	betaInv.Inverse(&beta)
	chi.Mul(&conjBeta, &betaInv)

	// chi should be one of: (1,0), (0,1), (-1,0), (0,-1)
	if chi.A1.IsZero() {
		// chi is real: either 1 or -1
		if chi.A0.IsOne() {
			return 0
		}
		return 2
	}
	if chi.A0.IsZero() {
		// chi is purely imaginary: i or -i
		if chi.A1.IsOne() {
			return 1
		}
		return 3
	}

	// Shouldn't reach here for valid inputs
	return 0
}
