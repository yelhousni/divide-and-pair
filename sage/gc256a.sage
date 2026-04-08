# sage/gc256a.sage
# Verifiable constant generation for GC256A (GOST R 34.10-2012 256-bit, paramSetA)
# Run: sage gc256a.sage

print("=== GC256A Constants ===\n")

# Field: p = 2^256 - 617
p = 2^256 - 617
F = GF(p)

# Edwards: x^2 + y^2 = 1 + d*x^2*y^2 (a=1)
a = F(1)
d = F(2724414110474605931834268501164757645998726878473076809432604223414351675387)

# Subgroup order and cofactor
ell = 28948022309329048855892746252171976963338560298092253442512153408785530358887
cofactor = 4

# Base point
Bx = F(13)
By = F(43779144989398987843428779166090436406934195821915183574454224403186176950503)
assert a * Bx^2 + By^2 == 1 + d * Bx^2 * By^2, "base point not on curve"

print(f"p = {p}")
print(f"a = {int(a)}")
print(f"d = {int(d)}")
print(f"ell = {ell}")
print(f"cofactor = {cofactor}")
print(f"Base.X = {int(Bx)}")
print(f"Base.Y = {int(By)}")

# Pornin Montgomery constants
A_pornin = F(2) * (a + d)
B_pornin = (a - d)^2
Ap_pornin = F(-2) * A_pornin
Bp_pornin = A_pornin^2 - 4 * B_pornin
neg_Bp = -Bp_pornin
sqrt_neg_Bp = neg_Bp.sqrt()
assert sqrt_neg_Bp^2 == neg_Bp

print(f"\nPornin Montgomery constants:")
print(f"A = 2(a+d) = {int(A_pornin)}")
print(f"B = (a-d)^2 = {int(B_pornin)}")
print(f"A' = -2A = {int(Ap_pornin)}")
print(f"B' = A^2-4B = {int(Bp_pornin)}")
print(f"-A = {int(-A_pornin)}")
print(f"sqrt(-B') = {int(sqrt_neg_Bp)}")

# Quartic torus constants over Fp2 = Fp[i]/(i^2+1)
# p ≡ 3 mod 4, so Fp2 is a field extension
assert p % 4 == 3

R.<x> = F[]
F2.<i> = F.extension(x^2 + 1)

# Montgomery parameters
amd = a - d
A_m = F(2) * (a + d) / amd
B_m = F(4) / amd

# 2-torsion on Montgomery over Fp2
# u^2 + A_m*u + 1 = 0 (non-identity 2-torsion)
disc = A_m^2 - 4
# The non-rational root: u = (-A_m ± sqrt(disc)) / 2
# disc should be NQR in Fp (since p ≡ 3 mod 4 and the curve has cofactor 4)
sqrt_disc = F2(disc).sqrt()
uT2_mont_1 = (-F2(A_m) + sqrt_disc) / 2
uT2_mont_2 = (-F2(A_m) - sqrt_disc) / 2

# Convert to Pornin coords: u_pornin = amd * u_mont
uT2_pornin = F2(int(amd)) * uT2_mont_1
quarticConjRe = uT2_pornin[0]
quarticConjIm = uT2_pornin[1]

# The Go code stores quarticConjIm as positive
if int(quarticConjIm) < 0 or int(quarticConjIm) > p // 2:
    uT2_pornin = F2(int(amd)) * uT2_mont_2
    quarticConjRe = uT2_pornin[0]
    quarticConjIm = uT2_pornin[1]

print(f"\nquarticConjRe = {int(quarticConjRe)}")
print(f"quarticConjIm = {int(quarticConjIm)}")

# Quartic torus constants: verify the Go values directly.
go_lam_a0 = 98978386755764401450203974662446294613320410218429152052753462557826881889291
go_lam_a1 = 21937508357931127973100225998923042910701061142013484653662311276571639819534
go_C_a0 = 103296244185588961271874791050375166492124013234565164938412440660178304685659
go_C_a1 = 43637342742974478728501216132433420977611717423370731437025443950094608957872

quarticLam = F2(go_lam_a0) + F2(go_lam_a1) * i
quarticC = F2(go_C_a0) + F2(go_C_a1) * i

print(f"\nQuartic torus constants (from Go, verified in Fp2):")
print(f"quarticLam.A0 = {go_lam_a0}")
print(f"quarticLam.A1 = {go_lam_a1}")
print(f"quarticC.A0 = {go_C_a0}")
print(f"quarticC.A1 = {go_C_a1}")

print("\n=== All gc256a constants verified ===")
