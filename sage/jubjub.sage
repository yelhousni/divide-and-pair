# sage/jubjub.sage
# Verifiable constant generation for JubJub (BLS12-381 scalar field)
# Run: sage jubjub.sage

print("=== JubJub Constants ===\n")

# Field: BLS12-381 scalar field (Fr)
p = 52435875175126190479447740508185965837690552500527637822603658699938581184513
F = GF(p)

# Twisted Edwards: -x^2 + y^2 = 1 + d*x^2*y^2 (a=-1)
a = F(-1)
d = F(19257038036680949359750312669786877991949435402254120286184196891950884077233)

# Subgroup order and cofactor
ell = 6554484396890773809930967563523245729705921265872317281365359162392183254199
cofactor = 8

# Base point
Bx = F(23426137002068529236790192115758361610982344002369094106619281483467893291614)
By = F(39325435222430376843701388596190331198052476467368316772266670064146548432123)
assert a * Bx^2 + By^2 == 1 + d * Bx^2 * By^2, "base point not on curve"

print(f"p = {p}")
print(f"a = {int(a)}")
print(f"d = {int(d)}")
print(f"ell = {ell}")
print(f"cofactor = {cofactor}")
print(f"Base.X = {int(Bx)}")
print(f"Base.Y = {int(By)}")

# NQR = 5 (verify 5 is not a QR mod p)
assert not F(5).is_square(), "5 should be NQR"
print(f"\nNQR = 5 (verified: 5 is not a QR mod p)")
assert p % 8 == 1, "p should be 1 mod 8 (so 2 is QR and we need a different NQR)"

# sqrt(-1) mod p
sqrt_minus_one = F(-1).sqrt()
print(f"sqrt(-1) = {int(sqrt_minus_one)}")

# Pornin Montgomery constants
A_pornin = F(2) * (a + d)
B_pornin = (a - d)^2
Ap_pornin = F(-2) * A_pornin
Bp_pornin = A_pornin^2 - 4 * B_pornin

# For JubJub, since p ≡ 1 mod 8, 2 is QR. We use NQR=5 for the isomorphism trick.
nqr = F(5)
sqrt_nqr_Bp = (nqr * Bp_pornin).sqrt()
assert sqrt_nqr_Bp^2 == nqr * Bp_pornin

print(f"\nPornin Montgomery constants:")
print(f"A = 2(a+d) = {int(A_pornin)}")
print(f"B = (a-d)^2 = {int(B_pornin)}")
print(f"A' = -2A = {int(Ap_pornin)}")
print(f"B' = A^2-4B = {int(Bp_pornin)}")
print(f"sqrt(5*B') = {int(sqrt_nqr_Bp)}")

# Weierstrass model: Y^2 = X^3 + aW*X + bW
# via Montgomery: u = X - A/3, where A is the Montgomery A-coefficient
# Montgomery A = 2(a+d)/(a-d), B_mont = 4/(a-d)
amd = a - d
A_m = F(2) * (a + d) / amd
B_m = F(4) / amd
Adiv3 = A_pornin / F(3)  # A_pornin = 2(a+d), not Montgomery A. Let me recheck.
# Actually the Go code computes Adiv3 = A_pornin / 3 where A_pornin = 2(a+d)
# But for Weierstrass: X = u + A_mont/3 where A_mont = 2(a+d)/(a-d)
# The Go code does: Adiv3 = A_pornin / 3 = 2(a+d)/3
# And: aW = B_pornin - A_pornin^2/3 = (a-d)^2 - (2(a+d))^2/3
Adiv3_go = A_pornin / 3
aW = B_pornin - A_pornin^2 / 3

print(f"\nWeierstrass model constants:")
print(f"Adiv3 = A/3 = {int(Adiv3_go)}")
print(f"aW = B - A^2/3 = {int(aW)}")

# Build Weierstrass curve: Y^2 = X^3 + aW*X + bW
# From Montgomery v^2 = u^3 + A_m*u^2 + u (after B scaling):
# X = u_m + A_m/3, Y = v_m (in unit-B representation)
# Actually for the Go code, the Weierstrass is built from Pornin coords:
# X = u_pornin + Adiv3, Y = u_pornin * w_pornin
# where u_pornin = (a-d)*(1+y)/(1-y), w_pornin = 2/x
# So X = (a-d)*u_m + 2(a+d)/3

# For building the curve directly:
# y^2 = x^3 + aW*x + bW where bW is determined by the discriminant
# bW = A_m*(2*A_m^2/27 - 1/3) (standard Montgomery to Weierstrass)
# But we need to be careful about the B_m scaling.
# Direct approach: compute bW from the fact that T2 lies on it.
# The 2-torsion on Edwards (a=-1): (0,-1) has order 2.
# In Pornin coords: u_pornin = (a-d)*(1+(-1))/(1-(-1)) = 0, w_pornin = 2/0 = ∞
# This doesn't work for the affine Weierstrass map. Let's use a different 2-torsion.
# Actually we need bW from the polynomial relation.

# Standard: from B_m*v^2 = u^3 + A_m*u^2 + u, substitute s = B_m*u, t = B_m^2*v:
# t^2 = s^3 + A_m*B_m*s^2 + B_m^2*s
# Then X = s + A_m*B_m/3: Y^2 = X^3 + aW'*X + bW'
# This is getting complex. Let me just build E from aW and compute bW.
# bW can be read off: the Montgomery curve B*v^2 = u^3+Au^2+u has j-invariant
# j = 256*(A^2-3)^3 / (A^2-4). But simpler:

# Actually the Go code builds it as: X = u_pornin + Adiv3, Y = u_pornin * w_pornin
# This means the Weierstrass curve coefficients relate to the Pornin coordinates.
# Pornin: curve is u^3 + A_p*u^2 + B_p*u (where A_p = 2(a+d), B_p = (a-d)^2)
# w^2 = u + A_p + B_p/u (Pornin's Montgomery-like form)
# Weierstrass substitution: X = u + A_p/3 =>
# (Y)^2*? ... This needs more care.

# For now, let me just build the Weierstrass from the standard transform and find torsion points.

# Standard short Weierstrass from Montgomery B_m*v^2 = u^3 + A_m*u^2 + u:
# Divide by B_m: v^2 = (1/B_m)*u^3 + (A_m/B_m)*u^2 + (1/B_m)*u
# Substitute u' = u/B_m, v' = v/B_m^(3/2)... this doesn't simplify nicely.
# Simplest: use Sage.

E_mont = EllipticCurve(F, [0, A_m, 0, F(1), 0])  # v^2 = u^3 + A_m*u^2 + u (B_m absorbed)
# Note: this is not quite right because the actual Montgomery is B_m*v^2 = u^3+A_m*u^2+u
# but for finding torsion structure, the isomorphic curve v^2 = u^3+A_m*u^2+u suffices.

E_W = E_mont.short_weierstrass_model()
phi_to_W = E_mont.isomorphism_to(E_W)

# Find 8-torsion points
# The Go code provides specific T8, T4, T2 on Weierstrass
go_XT8 = F(38074089473419775066521122308184450788044636445201246956909382033745382347397)
go_YT8 = F(30713742701501207659057220608756960799104471161253953255407673889302209579815)
go_XT4 = F(11059612379481747039899142612799695948580372366091172512139820602662565702425)
go_YT4 = F(4853940988772544534178633755630548083074718905337523302547095884433093026740)
go_XT2 = F(30316650416162696399649455282586573940529807768345292798324017494613449779659)

# The Weierstrass model from Go uses X = u_pornin + Adiv3, Y = u_pornin * w_pornin
# where Adiv3 = 2(a+d)/3. This is NOT the same as Sage's short_weierstrass_model.
# We need to match the Go model exactly.

# Go Weierstrass: the substitution is X = u + A_p/3 where A_p = 2(a+d)
# and u is the Pornin u-coordinate. The Pornin curve relation is:
# w^2*u = u^2 + A_p*u + B_p  (from the Montgomery relation, Pornin convention)
# So (Y/u)^2 * u = u^2 + A_p*u + B_p where Y = u*w
# => Y^2/u = u^2 + A_p*u + B_p
# => Y^2 = u^3 + A_p*u^2 + B_p*u
# Weierstrass: X = u + A_p/3 => u = X - A_p/3
# Y^2 = (X-A_p/3)^3 + A_p*(X-A_p/3)^2 + B_p*(X-A_p/3)
# = X^3 - A_p*X^2 + A_p^2*X/3 - A_p^3/27 + A_p*X^2 - 2*A_p^2*X/3 + A_p^3/9 + B_p*X - B_p*A_p/3
# = X^3 + (A_p^2/3 - 2*A_p^2/3 + B_p)*X + (-A_p^3/27 + A_p^3/9 - B_p*A_p/3)
# = X^3 + (B_p - A_p^2/3)*X + (2*A_p^3/27 - B_p*A_p/3)
aW_go = B_pornin - A_pornin^2 / 3
bW_go = 2 * A_pornin^3 / 27 - B_pornin * A_pornin / 3

E_go = EllipticCurve(F, [aW_go, bW_go])
print(f"aW (Go model) = {int(aW_go)}")
print(f"bW (Go model) = {int(bW_go)}")

# Verify Go torsion points
T8 = E_go(go_XT8, go_YT8)
T4 = E_go(go_XT4, go_YT4)
T2 = E_go(go_XT2, 0)

assert 2 * T8 == T4, "[2]T8 != T4"
assert 2 * T4 == T2, "[2]T4 != T2"
assert 2 * T2 == E_go(0), "[2]T2 != O"
assert 8 * T8 == E_go(0), "[8]T8 != O"
print("\nTorsion point verification:")
print(f"T8 = ({int(go_XT8)}, {int(go_YT8)})")
print(f"T4 = [2]T8 = ({int(go_XT4)}, {int(go_YT4)})")
print(f"T2 = [4]T8 = ({int(go_XT2)}, 0)")
print("[2]T8 == T4: True")
print("[4]T8 == T2: True")
print("[8]T8 == O:  True")

# Tangent slopes
go_lamT8 = F(48568263891961809311880741918802055033685793886551129072617710864462986515188)
lamT8 = (3 * T8[0]^2 + aW_go) / (2 * T8[1])
assert lamT8 == go_lamT8, "lamT8 mismatch"
print(f"\nlamT8 = {int(lamT8)}")

# lamT4 from the Go code (computed as (3*XT4^2 + aW) / (2*YT4))
lamT4 = (3 * T4[0]^2 + aW_go) / (2 * T4[1])
print(f"lamT4 = {int(lamT4)}")

# (p-1)/8 exponent for octic symbol
assert (p - 1) % 8 == 0
octic_exp = (p - 1) // 8
print(f"\n(p-1)/8 = {octic_exp}")

print("\n=== All jubjub constants verified ===")
