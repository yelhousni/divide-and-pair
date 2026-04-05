# Divide-and-Pair: Subgroup Membership Testing

Go implementation of fast subgroup membership testing on twisted Edwards curves using the **Divide-and-Pair** framework, which combines point halvings with Tate pairing-based residuosity checks.

## Curves

| Curve | Field | a | Cofactor | d = gcd(h, p−1) | Halvings | Final check |
|-------|-------|---|----------|------------------|----------|-------------|
| **Ed25519** | p = 2²⁵⁵ − 19 | −1 | 8 | 4 | 1 | Quartic symbol |
| **JubJub** | BLS12-381 Fr | −1 | 8 | 4 | 1 | Quartic symbol |
| **Curve448** | p = 2⁴⁴⁸ − 2²²⁴ − 1 | 1 | 4 | 2 | 1 | Legendre symbol |
| **GC256A** | p = 2²⁵⁶ − 617 | 1 | 4 | 2 | 1 | Legendre symbol |

## Methods

- **Naive**: scalar multiplication by the subgroup order ℓ and checking for the identity point.
- **Pornin** ([ePrint 2022/1164](https://eprint.iacr.org/2022/1164)): division-free halvings on an isogenous Montgomery curve followed by a Legendre symbol check. Implementation follows [crrl](https://github.com/pornin/crrl).
- **Quartic** (Ed25519): replaces one halving + Legendre with a single quartic residue symbol, computed via Weilert's Euclidean algorithm over Z[i].
- **QuarticExp** (Ed25519, JubJub): same as Quartic but uses addition-chain exponentiation for the quartic symbol.

## Benchmarks

Apple M5, Go 1.25, `go test -bench=.`:

### Ed25519
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
| **QuarticExp** | **18 µs** | **60×** |

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

## References

- Pornin, [*Point-Halving and Subgroup Membership in Twisted Edwards Curves*](https://eprint.iacr.org/2022/1164), 2022.
- Weilert, [*Fast Computation of the Biquadratic Residue Symbol*](https://doi.org/10.1007/s00145-002-0131-7), J. Cryptology, 2003.
- [RFC 8032](https://www.rfc-editor.org/rfc/rfc8032) — Ed25519 and Ed448.
- [RFC 7836](https://www.rfc-editor.org/rfc/rfc7836) — GC256A (GOST R 34.10-2012).
