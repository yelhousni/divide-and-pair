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
| **GC256A** | p = 2²⁵⁶ − 617 | 1 | 4 | 2 | 1 | Quadratic symbol |

## Methods

- **Naive**: scalar multiplication by the subgroup order ℓ and checking for the identity point.
- **Pornin** ([ePrint 2022/1164](https://eprint.iacr.org/2022/1164)): division-free halvings on an isogenous Montgomery curve followed by a Legendre symbol check. Implementation follows [crrl](https://github.com/pornin/crrl).
- **Quartic** (Curve25519): replaces one halving + Legendre with a single quartic residue symbol, computed via Weilert's Euclidean algorithm over Z[i].
- **QuarticExp** (Curve25519): same as Quartic but uses addition-chain exponentiation for the quartic symbol.
- **OcticExp** (JubJub): 0 halvings + 1 octic residuosity check. Since d = gcd(8, p−1) = 8 for the BLS12-381 scalar field (p ≡ 1 mod 8), no halvings are needed. The degree-8 Miller function is evaluated on the Weierstrass model via a 3-step doubling chain, and the octic symbol χ₈(f) = f^((p−1)/8) is checked division-free.
- **Tate** (FourQ): 0 divisions + octic check + septic check.
  - **Octic**: degree-8 Miller function on the Weierstrass model, then torus/trace compression: g = conj(f₈)/f₈ is projected to its trace t ∈ Fp (projectively, no Fp²-inversion), followed by 124 iterations of t → t²−2 (1 Fp-squaring per step). Check trace == 2.
  - **Septic**: degree-7 Miller function with precomputed intermediates. Norm-accumulated: Norm(f₇) = ∏ Norm(ℓᵢ) / ∏ Norm(vⱼ) computed entirely in Fp from individual line norms — no Fp² inversions needed. Then χ₇(f₇) = Norm(f₇)^((p−1)/7) via a single Fp exponentiation.
- **Quartic over Fp²** (Curve448): since p ≡ 3 mod 4, the quartic character does not exist over Fp, but does over Fp² (embedding degree k = 2 for the 4-torsion). A non-rational 4-torsion point T₄ ∈ E(Fp²)\E(Fp) is used to evaluate a degree-4 Tate pairing, reducing subgroup membership to a single quartic residuosity check in Fp². The torus approach is used: g = conj(α)/α is projected to its trace t ∈ Fp, then a PRAC differential addition chain (Montgomery 1992) evaluates the Lucas V-sequence V_e(t) for e = (p+1)/4. PRAC uses 632 field ops vs 890 for a binary Montgomery ladder (29% saving).
- **Endomorphism** (FourQ): subgroup membership from the FourQ ψ and φ endomorphisms, reducing the order check to a 4-dimensional multi-scalar relation with short scalars.

## Benchmarks

AWS r7a (AMD EPYC 9R14), Go 1.24, `go test -bench=. -benchtime=3s`:

### Curve25519
| Method | Time | Speedup vs Naive |
|--------|------|-----------------|
| Naive | 1,568 µs | 1× |
| Pornin | 36 µs | 44× |
| QuarticExp | 22 µs | 72× |
| **Quartic (GCD)** | **25 µs** | **62×** |

### JubJub
| Method | Time | Speedup vs Naive |
|--------|------|-----------------|
| Naive | 1,820 µs | 1× |
| Pornin | 45 µs | 40× |
| **OcticExp** | **9.7 µs** | **188×** |

### FourQ
| Method | Time | Speedup vs Naive |
|--------|------|-----------------|
| Naive | 908 µs | 1× |
| Endomorphism | 243 µs | 3.7× |
| **Tate (torus octic + Norm septic)** | **4.8 µs** | **190×** |

### Curve448
| Method | Time | Speedup vs Naive |
|--------|------|-----------------|
| Naive | 5,368 µs | 1× |
| Pornin (1 halving + Legendre) | 111 µs | 48× |
| **Quartic (torus/PRAC)** | **48 µs** | **112×** |

### GC256A
| Method | Time | Speedup vs Naive |
|--------|------|-----------------|
| Naive | 1,170 µs | 1× |
| **Pornin** | **21 µs** | **57×** |

## Usage

```go
import "github.com/yelhousni/divide-and-pair/curve25519"

g := curve25519.Generators()
var p curve25519.PointAffine
p.ScalarMultiplication(&g, k)

p.IsInSubGroup()           // 1 halving + quartic symbol (GCD)
```

```go
import "github.com/yelhousni/divide-and-pair/jubjub"

g := jubjub.Generators()
var p jubjub.PointAffine
p.ScalarMultiplication(&g, k)

p.IsInSubGroup()           // 0 halvings + octic symbol (exp)
```

```go
import "github.com/yelhousni/divide-and-pair/fourq"

g := fourq.Generators()
var p fourq.PointAffine
p.ScalarMultiplication(&g, k)

p.IsInSubGroup()           // torus octic + Norm septic
```

```go
import "github.com/yelhousni/divide-and-pair/curve448"

g := curve448.Generators()
var p curve448.PointAffine
p.ScalarMultiplication(&g, k)

p.IsInSubGroup()           // quartic symbol over Fp² (torus/PRAC)
```

## References

- Pornin, [*Point-Halving and Subgroup Membership in Twisted Edwards Curves*](https://eprint.iacr.org/2022/1164), 2022.
- Koshelev, [*Subgroup membership testing on elliptic curves via the Tate pairing*](https://eprint.iacr.org/2022/037), 2022.
- Costello and Longa, [*FourQ: four-dimensional decompositions on a Q-curve over the Mersenne prime*](https://eprint.iacr.org/2015/565), ASIACRYPT 2015.
- Weilert, [*Fast Computation of the Biquadratic Residue Symbol*](https://doi.org/10.1007/s00145-002-0131-7), J. Cryptology, 2003.
- Montgomery, [*Speeding the Pollard and Elliptic Curve Methods of Factorization*](https://doi.org/10.1090/S0025-5718-1987-0866113-7), Math. Comp., 1987. (PRAC chains)
- [RFC 8032](https://www.rfc-editor.org/rfc/rfc8032) — Curve25519 and Ed448.
- [RFC 7836](https://www.rfc-editor.org/rfc/rfc7836) — GC256A (GOST R 34.10-2012).
