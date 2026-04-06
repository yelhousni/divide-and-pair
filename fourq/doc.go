// Package fourq implements the FourQ twisted Edwards curve
// (-x² + y² = 1 + d·x²·y² with a = -1) over Fp² where p = 2¹²⁷ - 1
// (Mersenne prime) and Fp² = Fp[i]/(i² + 1).
// The curve has cofactor 392 = 2³ × 7² and prime subgroup order ≈ 2²⁴⁶.
package fourq
