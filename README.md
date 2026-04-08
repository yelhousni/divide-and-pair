# Divide-and-Pair: Faster subgroup membership testing for elliptic curves

[![License](https://img.shields.io/badge/license-MIT-blue)](LICENSE) [![Go Report Card](https://goreportcard.com/badge/github.com/yelhousni/divide-and-pair)](https://goreportcard.com/report/github.com/yelhousni/divide-and-pair) [![Go Reference](https://pkg.go.dev/badge/github.com/yelhousni/divide-and-pair.svg)](https://pkg.go.dev/github.com/yelhousni/divide-and-pair)

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="ladder-dark.svg">
    <img src="ladder-light.svg" width="200">
  </picture>
</p>

Companion code for the article *"Divide-and-Pair: Faster subgroup membership testing for elliptic curves"* by Y. Dai, Y. El Housni, D. Koshelev and K. Reijinders.

## Curves

| Curve | Field | a | Cofactor | d = gcd(h, q−1) | Halvings | Final check |
|-------|-------|---|----------|------------------|----------|-------------|
| **Curve25519** | p = 2²⁵⁵ − 19 | −1 | 8 | 4 | 1 | Quartic symbol |
| **JubJub** | BLS12-381 Fr | −1 | 8 | 8 | 0 | Octic symbol |
| **FourQ** | Fp² (p = 2¹²⁷ − 1) | −1 | 392 | 392 | 0 | Octic + septic symbols|
| **Curve448** | p = 2⁴⁴⁸ − 2²²⁴ − 1 | 1 | 4 | 4 (over Fp²) | 0 | Quartic symbol (Fp²) |
| **GC256A** | p = 2²⁵⁶ − 617 | 1 | 4 | 4 (over Fp²) | 0 | Quartic symbol (Fp²) |

## Methods

- **Naive**: scalar multiplication by the subgroup order ℓ and checking for the identity point.
- **Pornin** ([ePrint 2022/1164](https://eprint.iacr.org/2022/1164)): division-free halvings on an isogenous Montgomery curve followed by a Legendre symbol check. Implementation follows [crrl](https://github.com/pornin/crrl).
- **Quartic** (Curve25519): replaces one halving + Legendre with a single quartic residue symbol, computed via Weilert's Euclidean algorithm over Z[i].
- **QuarticExp** (Curve25519): same as Quartic but uses addition-chain exponentiation for the quartic symbol.
- **OcticExp** (JubJub): 0 halvings + 1 octic residuosity check. Since d = gcd(8, p−1) = 8 for the BLS12-381 scalar field (p ≡ 1 mod 8), no halvings are needed. The degree-8 Miller function is evaluated on the Weierstrass model via a 3-step doubling chain, and the octic symbol χ₈(f) = f^((p−1)/8) is checked division-free.
- **Tate** (FourQ): 0 divisions + octic check + septic check.
  - **Octic**: degree-8 Miller function on the Weierstrass model, then torus/trace compression: g = conj(f₈)/f₈ is projected to its trace t ∈ Fp (projectively, no Fp²-inversion), followed by 124 iterations of t → t²−2 (1 Fp-squaring per step). Check trace == 2.
  - **Septic**: degree-7 Miller function with precomputed intermediates. Norm-accumulated: Norm(f₇) = ∏ Norm(ℓᵢ) / ∏ Norm(vⱼ) computed entirely in Fp from individual line norms — no Fp² inversions needed. Then χ₇(f₇) = Norm(f₇)^((p−1)/7) via a single Fp exponentiation.
- **Quartic over Fp²** (Curve448, GC256A): since p ≡ 3 mod 4, the quartic character does not exist over Fp, but does over Fp² (embedding degree k = 2 for the 4-torsion). A non-rational 4-torsion point T₄ ∈ E(Fp²)\E(Fp) is used to evaluate a degree-4 Tate pairing, reducing subgroup membership to a single quartic residuosity check in Fp². The torus approach is used: g = conj(α)/α is projected to its trace t ∈ Fp, then a PRAC differential addition chain (Montgomery 1992) evaluates the Lucas V-sequence V_e(t) for e = (p+1)/4.
- **Endomorphism** (FourQ): subgroup membership from the FourQ ψ and φ endomorphisms, reducing the order check to a 4-dimensional multi-scalar relation with short scalars.

## Benchmarks

AWS r7a (AMD EPYC 9R14), Go 1.24, `go test -bench=. -benchtime=3s`.
Curve25519 and JubJub use amd64 assembly for field arithmetic.

### Curve25519
| Method | Time | Speedup vs Naive |
|--------|------|-----------------|
| Naive | 1,563 µs | 1× |
| Pornin | 27 µs | 58× |
| QuarticExp | 16 µs | 98× |
| **Quartic (GCD)** | **15 µs** | **104×** |

### JubJub
| Method | Time | Speedup vs Naive |
|--------|------|-----------------|
| Naive | 1,762 µs | 1× |
| Pornin | 41 µs | 43× |
| **OcticExp** | **8.7 µs** | **203×** |

### FourQ
| Method | Time | Speedup vs Naive |
|--------|------|-----------------|
| Naive | 908 µs | 1× |
| Endomorphism | 242 µs | 3.8× |
| **Tate (torus octic + Norm septic)** | **4.7 µs** | **193×** |

### Curve448
| Method | Time | Speedup vs Naive |
|--------|------|-----------------|
| Naive | 5,364 µs | 1× |
| Pornin (1 halving + Legendre) | 111 µs | 48× |
| **Quartic (torus/PRAC)** | **48 µs** | **112×** |

### GC256A
| Method | Time | Speedup vs Naive |
|--------|------|-----------------|
| Naive | 1,176 µs | 1× |
| Pornin (1 halving + Legendre) | 25 µs | 47× |
| **Quartic (torus/PRAC)** | **16 µs** | **72×** |

## References

- Pornin, [*Point-Halving and Subgroup Membership in Twisted Edwards Curves*](https://eprint.iacr.org/2022/1164), 2022.
- Koshelev, [*Subgroup membership testing on elliptic curves via the Tate pairing*](https://eprint.iacr.org/2022/037), 2022.
- Costello and Longa, [*FourQ: four-dimensional decompositions on a Q-curve over the Mersenne prime*](https://eprint.iacr.org/2015/565), ASIACRYPT 2015.
- Weilert, [*Fast Computation of the Biquadratic Residue Symbol*](https://doi.org/10.1007/s00145-002-0131-7), J. Cryptology, 2003.
- Montgomery, [*Evaluating recurrences of form Xm+n = f(Xm, Xn, Xm-n) via Lucas chains*](https://cr.yp.to/bib/1992/montgomery-lucas.pdf). (PRAC chains)
- [RFC 8032](https://www.rfc-editor.org/rfc/rfc8032) — Curve25519 and Ed448.
- [RFC 7836](https://www.rfc-editor.org/rfc/rfc7836) — GC256A (GOST R 34.10-2012).
- McLoughlin, [addchain](https://github.com/mmcloughlin/addchain). Software to generate short addition chains (Go).
- Bernstein,  https://cr.yp.to/2024/dacbench-20240609.tar.gz. Software to generate PRAC differential addition chains (Python).
