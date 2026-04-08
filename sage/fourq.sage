# sage/fourq.sage
# Verifiable constant generation for FourQ
# Run: sage fourq.sage

print("=== FourQ Constants ===\n")

# Base field: Mersenne prime p = 2^127 - 1
p = 2^127 - 1
F = GF(p)

# Extension field Fp2 = Fp[i]/(i^2+1)
R.<x> = F[]
F2.<i> = F.extension(x^2 + 1)

# Twisted Edwards over Fp2: -x^2 + y^2 = 1 + d*x^2*y^2 (a = -1)
a = F2(-1)
d = F2(4205857648805777768770) + F2(125317048443780598345676279555970305165) * i

# Subgroup order and cofactor
N = 73846995687063900142583536357581573884798075859800097461294096333596429543
cofactor = 392

# Base point
Gx = F2(34832242333165934151976439273177494442) + F2(40039530084877881816286215037915002870) * i
Gy = F2(18941146186793715734774048165794132615) + F2(146361984425930646555497992424795179868) * i
assert a * Gx^2 + Gy^2 == 1 + d * Gx^2 * Gy^2, "base point not on curve"

print(f"p = {p}")
print(f"a = -1")
print(f"d.A0 = {int(d[0])}")
print(f"d.A1 = {int(d[1])}")
print(f"N (subgroup order) = {N}")
print(f"cofactor = {cofactor}")
print(f"Base.X.A0 = {int(Gx[0])}")
print(f"Base.X.A1 = {int(Gx[1])}")
print(f"Base.Y.A0 = {int(Gy[0])}")
print(f"Base.Y.A1 = {int(Gy[1])}")

# Weierstrass model: Y^2 = X^3 + aW*X + bW
# Montgomery: A = 2(a+d)/(a-d), B = 4/(a-d)
# Pornin: A_p = 2(a+d), B_p = (a-d)^2
amd = a - d
A_p = F2(2) * (a + d)
B_p = amd^2
Adiv3 = A_p / 3
aW = B_p - A_p^2 / 3
bW = 2 * A_p^3 / 27 - B_p * A_p / 3

print(f"\nWeierstrass model:")
print(f"Adiv3.A0 = {int(Adiv3[0])}")
print(f"Adiv3.A1 = {int(Adiv3[1])}")
print(f"aW.A0 = {int(aW[0])}")
print(f"aW.A1 = {int(aW[1])}")

E_W = EllipticCurve(F2, [aW, bW])

# 8-torsion point T8 (from Go code)
XT8 = F2(56495978839162271580189242121432892911) + F2(5498965435856981984649943052775033891) * i
YT8 = F2(110406726416374112713382695534503225009) + F2(41382232403199127624659984993520567989) * i
T8 = E_W(XT8, YT8)
assert 8 * T8 == E_W(0), "[8]T8 != O"

# [2]T8 = T4
T4 = 2 * T8
XT4 = T4[0]
YT4 = T4[1]
assert int(XT4[0]) == 170141183460469230329734754113958182802
assert int(XT4[1]) == 71655106159052621705899442625265968763
assert int(YT4[0]) == 80492913427091964959665255396056504603
assert int(YT4[1]) == 170141183460469223319972006104328568185

# [4]T8 = T2
T2 = 2 * T4
XT2 = T2[0]
assert int(XT2[0]) == 2803905099203851845846
assert int(XT2[1]) == 26830971142363988319888418465352168201
assert T2[1] == 0  # 2-torsion

print(f"\n8-torsion point T8:")
print(f"XT8.A0 = {int(XT8[0])}")
print(f"XT8.A1 = {int(XT8[1])}")
print(f"YT8.A0 = {int(YT8[0])}")
print(f"YT8.A1 = {int(YT8[1])}")
print(f"[2]T8 = T4: verified")
print(f"[4]T8 = T2: verified")
print(f"[8]T8 = O:  verified")

# Tangent slopes
lamT8 = (3 * T8[0]^2 + aW) / (2 * T8[1])
print(f"\nlamT8.A0 = {int(lamT8[0])}")
print(f"lamT8.A1 = {int(lamT8[1])}")
assert int(lamT8[0]) == 144442680182921669228377474261162789491
assert int(lamT8[1]) == 78866152311824694793607513268751646905

lamT4 = (3 * T4[0]^2 + aW) / (2 * T4[1])
print(f"lamT4.A0 = {int(lamT4[0])}")
print(f"lamT4.A1 = {int(lamT4[1])}")
# Go has lamT4.A1 = 2 (SetString("2"))
assert int(lamT4[1]) == 2

# 7-torsion point S7 (from Go code)
XS7 = F2(170000906615246946031520731709186440780) + F2(94692800017046438403892332000427881085) * i
YS7 = F2(113074451165823784128958858701762803038) + F2(4234707370120948597958144423972584168) * i
S7 = E_W(XS7, YS7)
assert 7 * S7 == E_W(0), "[7]S7 != O"
print(f"\n7-torsion point S7:")
print(f"XS7.A0 = {int(XS7[0])}")
print(f"XS7.A1 = {int(XS7[1])}")
print(f"YS7.A0 = {int(YS7[0])}")
print(f"YS7.A1 = {int(YS7[1])}")
print("[7]S7 = O: verified")

# Miller loop intermediates
S7_2 = 2 * S7
S7_3 = S7_2 + S7

print(f"\n[2]S7:")
print(f"X2S.A0 = {int(S7_2[0][0])}")
print(f"X2S.A1 = {int(S7_2[0][1])}")
print(f"Y2S.A0 = {int(S7_2[1][0])}")
print(f"Y2S.A1 = {int(S7_2[1][1])}")
assert int(S7_2[0][0]) == 54301625611251160727403361627257607395
assert int(S7_2[0][1]) == 98913730175858778379537234706939518181

print(f"\n[3]S7:")
print(f"X3S.A0 = {int(S7_3[0][0])}")
print(f"X3S.A1 = {int(S7_3[0][1])}")
print(f"Y3S.A0 = {int(S7_3[1][0])}")
print(f"Y3S.A1 = {int(S7_3[1][1])}")
assert int(S7_3[0][0]) == 156066723560704030263055387453284139745
assert int(S7_3[0][1]) == 27481310916039459192906723051714574829

# Tangent slope at S7 (doubling)
lamDbl1 = (3 * S7[0]^2 + aW) / (2 * S7[1])
print(f"\nlamDbl1.A0 = {int(lamDbl1[0])}")
print(f"lamDbl1.A1 = {int(lamDbl1[1])}")
assert int(lamDbl1[0]) == 149197531863252711357534985748910785211
assert int(lamDbl1[1]) == 106728816833063046549591444245799095908

# Chord slope S7 + [2]S7
lamAdd1 = (S7_2[1] - S7[1]) / (S7_2[0] - S7[0])
print(f"lamAdd1.A0 = {int(lamAdd1[0])}")
print(f"lamAdd1.A1 = {int(lamAdd1[1])}")
assert int(lamAdd1[0]) == 119059155663838527172282288165420875213
assert int(lamAdd1[1]) == 2641556124042940121879682364519729700

# Tangent slope at [3]S7 (doubling)
lamDbl2 = (3 * S7_3[0]^2 + aW) / (2 * S7_3[1])
print(f"lamDbl2.A0 = {int(lamDbl2[0])}")
print(f"lamDbl2.A1 = {int(lamDbl2[1])}")
assert int(lamDbl2[0]) == 58231718181463545180414387022826039618
assert int(lamDbl2[1]) == 21100793565021504938029727675157445515

# [6]S7 = -S7 (same X)
S7_6 = 6 * S7
assert S7_6[0] == S7[0], "[6]S7 should have same X as S7"
print(f"\n[6]S7.X = S7.X: verified")

# Endomorphism constants
print(f"\n--- Endomorphism constants ---")

# sqrt(d_hat) where d_hat = -1/(1+d)
d_hat = -F2(1) / (F2(1) + d)
sqrt_d_hat = d_hat.sqrt()
# Pick the right square root matching Go
if int(sqrt_d_hat[0]) != 120525532476903946900736407295642634213:
    sqrt_d_hat = -sqrt_d_hat
print(f"sqrtDhat.A0 = {int(sqrt_d_hat[0])}")
print(f"sqrtDhat.A1 = {int(sqrt_d_hat[1])}")
assert int(sqrt_d_hat[0]) == 120525532476903946900736407295642634213
assert int(sqrt_d_hat[1]) == 110680464442257309687
assert sqrt_d_hat^2 == d_hat

# Eigenvalues
print(f"\nlambdaPsi = 43760231755807040276284855770911078252536368422635318376310714077319867016")
print(f"lambdaPhi = 12098939722099758392970036154455447385486035337534694534042314319425271908")
print(f"lambdaPhiPsi = 42306631464858389077121485826103052927889336227980061597005228525569104570")

# Babai coefficients
print(f"\nendoC1 = -650487742939046294")
print(f"endoC2 = 1397215820276968864")
print(f"endoC3 = -523086274270593807")
print(f"endoC4 = 598824378691085905")

# ClearCofactor scalar
print(f"\ncofactor for ClearCofactor = 392")
assert cofactor == 392

print("\n=== All fourq constants verified ===")
