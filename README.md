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
  - **Octic**: degree-8 Miller function on the Weierstrass model. Three approaches implemented:
    - *Exp1* (β approach): β = f₈^((p+1)/8), check β.A1 == 0. Avoids Fp²-inversion by turning the Frobenius into a free zero-check.
    - *Exp2* (g approach): g = conj(f₈)/f₈ on the norm-1 torus, then cyclotomic squarings (2S vs 2M per step).
    - *Exp3* (torus/trace): compress g to its trace t ∈ Fp, then 124 iterations of t → t²−2 (1S per step, projective).
  - **Septic**: degree-7 Miller function with precomputed intermediates. Optimized via Norm reduction: χ₇(f₇) = Norm(f₇)^((p−1)/7), computed as a single Fp exponentiation (vs Fp² in the original). Three Fp²-inversions batched into one via Montgomery's trick.
- **Quartic over Fp²** (Curve448): since p ≡ 3 mod 4, the quartic character does not exist over Fp, but does over Fp² (embedding degree k = 2 for the 4-torsion). A non-rational 4-torsion point T₄ ∈ E(Fp²)\E(Fp) is used to evaluate a degree-4 Tate pairing, reducing subgroup membership to a single quartic residuosity check in Fp². Four approaches implemented:
  - *Exp1* (β approach): β = α^((p+1)/4) in Fp², check β.A1 == 0.
  - *Exp2* (g approach): g = conj(α)/α, cyclotomic Fp² exponentiation.
  - *Exp3* (torus/Lucas): Montgomery ladder on the trace t ∈ Fp with t → t²−2 squarings.
  - *Exp4* (torus/PRAC): PRAC differential addition chain (Montgomery 1992) for the Lucas V-sequence. Uses 632 field ops vs 890 for the binary ladder (29% saving).
- **Endomorphism** (FourQ): subgroup membership from the FourQ ψ and φ endomorphisms, reducing the order check to a 4-dimensional multi-scalar relation with short scalars.

## Benchmarks

Apple M5, Go 1.25, `go test -bench=.`:

### Curve25519
| Method | Time | Speedup vs Naive |
|--------|------|-----------------|
| Naive | 893 µs | 1× |
| Pornin | 17 µs | 53× |
| QuarticExp | 16 µs | 55× |
| **Quartic (GCD)** | **15 µs** | **61×** |

### JubJub
| Method | Time | Speedup vs Naive |
|--------|------|-----------------|
| Naive | 968 µs | 1× |
| Pornin | 27 µs | 35× |
| **OcticExp** | **5.2 µs** | **186×** |

### FourQ
| Method | Time | Speedup vs Naive |
|--------|------|-----------------|
| Naive | 461 µs | 1× |
| Endomorphism | 133 µs | 3.5× |
| Tate (g + generic sq) | 3.2 µs | 143× |
| **Tate Exp1 (β approach)** | **2.9 µs** | **159×** |

### Curve448
| Method | Time | Speedup vs Naive |
|--------|------|-----------------|
| Naive | 3,445 µs | 1× |
| Pornin (1 halving + Legendre) | 52 µs | 66× |
| Quartic Exp1 (β, Fp² add chain) | 40 µs | 86× |
| Quartic Exp3 (torus, binary ladder) | 38 µs | 91× |
| **Quartic Exp4 (torus, PRAC)** | **30 µs** | **115×** |

### GC256A
| Method | Time | Speedup vs Naive |
|--------|------|-----------------|
| Naive | 780 µs | 1× |
| **Pornin** | **18 µs** | **44×** |

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

p.IsInSubGroup()           // octic (β approach) + septic (Norm + Fp exp)
```

```go
import "github.com/yelhousni/divide-and-pair/curve448"

g := curve448.Generators()
var p curve448.PointAffine
p.ScalarMultiplication(&g, k)

p.IsInSubGroup()           // quartic symbol over Fp² (PRAC torus)
```

## References

- Pornin, [*Point-Halving and Subgroup Membership in Twisted Edwards Curves*](https://eprint.iacr.org/2022/1164), 2022.
- Koshelev, [*Subgroup membership testing on elliptic curves via the Tate pairing*](https://eprint.iacr.org/2022/037), 2022.
- Costello and Longa, [*FourQ: four-dimensional decompositions on a Q-curve over the Mersenne prime*](https://eprint.iacr.org/2015/565), ASIACRYPT 2015.
- Weilert, [*Fast Computation of the Biquadratic Residue Symbol*](https://doi.org/10.1007/s00145-002-0131-7), J. Cryptology, 2003.
- Montgomery, [*Speeding the Pollard and Elliptic Curve Methods of Factorization*](https://doi.org/10.1090/S0025-5718-1987-0866113-7), Math. Comp., 1987. (PRAC chains)
- [RFC 8032](https://www.rfc-editor.org/rfc/rfc8032) — Curve25519 and Ed448.
- [RFC 7836](https://www.rfc-editor.org/rfc/rfc7836) — GC256A (GOST R 34.10-2012).
