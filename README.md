# Divide-and-Pair: Subgroup Membership Testing

Go implementation of fast subgroup membership testing on twisted Edwards curves using the **Divide-and-Pair** framework, which combines point halvings with Tate pairing-based residuosity checks.

## Curves

| Curve | Field | a | Cofactor | d = gcd(h, p−1) | Halvings | Final check |
|-------|-------|---|----------|------------------|----------|-------------|
| **Curve25519** | p = 2²⁵⁵ − 19 | −1 | 8 | 4 | 1 | Quartic symbol |
| **JubJub** | BLS12-381 Fr | −1 | 8 | 8 | 0 | Octic symbol |
| **Curve448** | p = 2⁴⁴⁸ − 2²²⁴ − 1 | 1 | 4 | 2 | 1 | Legendre symbol |
| **GC256A** | p = 2²⁵⁶ − 617 | 1 | 4 | 2 | 1 | Legendre symbol |

## Methods

- **Naive**: scalar multiplication by the subgroup order ℓ and checking for the identity point.
- **Pornin** ([ePrint 2022/1164](https://eprint.iacr.org/2022/1164)): division-free halvings on an isogenous Montgomery curve followed by a Legendre symbol check. Implementation follows [crrl](https://github.com/pornin/crrl).
- **Quartic** (Curve25519): replaces one halving + Legendre with a single quartic residue symbol, computed via Weilert's Euclidean algorithm over Z[i].
- **QuarticExp** (Curve25519): same as Quartic but uses addition-chain exponentiation for the quartic symbol.
- **OcticExp** (JubJub): 0 halvings + 1 octic residuosity check. Since d = gcd(8, p−1) = 8 for the BLS12-381 scalar field (p ≡ 1 mod 8), no halvings are needed. The degree-8 Miller function is evaluated on the Weierstrass model via a 3-step doubling chain, and the octic symbol χ₈(f) = f^((p−1)/8) is checked division-free.

## Benchmarks

Apple M5, Go 1.25, `go test -bench=.`:

### Curve25519
| Method | Time | Speedup vs Naive |
|--------|------|-----------------|
| Naive | 1,065 µs | 1× |
| Pornin | 27 µs | 39× |
| QuarticExp | 19 µs | 56× |
| **Quartic (GCD)** | **17 µs** | **63×** |

### JubJub
| Method | Time | Speedup vs Naive |
|--------|------|-----------------|
| Naive | 1,072 µs | 1× |
| Pornin | 28 µs | 38× |
| **OcticExp** | **5.7 µs** | **188×** |

### Curve448
| Method | Time | Speedup vs Naive |
|--------|------|-----------------|
| Naive | 3,985 µs | 1× |
| **Pornin** | **81 µs** | **49×** |

### GC256A
| Method | Time | Speedup vs Naive |
|--------|------|-----------------|
| Naive | 867 µs | 1× |
| **Pornin** | **16 µs** | **55×** |

## Usage

```go
import "github.com/yelhousni/divide-and-pair/curve25519"

params := curve25519.GetEdwardsCurve()
var p curve25519.PointAffine
p.ScalarMultiplication(&params.Base, k)

p.IsInSubGroupNaive()      // scalar mult by ℓ
p.IsInSubGroupPornin()     // halvings + Legendre
p.IsInSubGroupQuartic()    // 1 halving + quartic symbol (GCD)
p.IsInSubGroupQuarticExp() // 1 halving + quartic symbol (exp)
```

```go
import "github.com/yelhousni/divide-and-pair/jubjub"

params := jubjub.GetEdwardsCurve()
var p jubjub.PointAffine
p.ScalarMultiplication(&params.Base, k)

p.IsInSubGroupOcticExp()   // 0 halvings + octic symbol (exp)
```

## References

- Pornin test, [*Point-Halving and Subgroup Membership in Twisted Edwards Curves*](https://eprint.iacr.org/2022/1164), 2022.
- Koshelev test, [*Subgroup membership testing on elliptic curves via the Tate pairing*](https://eprint.iacr.org/2022/037), 2022.
- Weilert, [*Fast Computation of the Biquadratic Residue Symbol*](https://doi.org/10.1007/s00145-002-0131-7), J. Cryptology, 2003.
- [RFC 8032](https://www.rfc-editor.org/rfc/rfc8032) — Curve25519 and Ed448.
- [RFC 7836](https://www.rfc-editor.org/rfc/rfc7836) — GC256A (GOST R 34.10-2012).
