# Benchmarks: naive subgroup membership test and cofactor clearing

Timings for `isInSubGroupNaive` (scalar multiplication by the subgroup order ℓ
via a fixed addition chain, followed by an identity check) and
`ClearCofactor` (multiplication by the cofactor h via a fixed chain of
doublings/additions, affine in/out).

## Setup

- **Machine:** AWS EC2, AMD EPYC 9R14 (Zen 4), 32 vCPUs, Ubuntu 24.04.3 LTS
  (kernel 6.17.0-1010-aws)
- **Toolchain:** Go 1.26.4 linux/amd64, built with `-tags purego`
  (the field arithmetic in this repository is pure Go; no assembly backend)
- **Command:**

  ```bash
  go test -tags purego -run '^$' \
    -bench 'BenchmarkIsInSubGroupNaive$|BenchmarkClearCofactor$' \
    -count 10 ./curve25519 ./jubjub ./fourq ./curve448 ./gc256a
  ```

- **Reported values:** benchstat means over 10 runs (variation ≤ ±4% on
  ClearCofactor, ≤ ±1% on the naive test).
- **Input points:** all benchmarks (`benchSubgroup` and
  `BenchmarkClearCofactor`) cycle over 256 random input points. Timing a
  single fixed point instead lets the branch predictor learn the
  data-dependent branch pattern of the field arithmetic and understates the
  cost on random inputs (by ~13% on curve25519 and up to ~85% on fourq).
- **Date:** 2026-07-13.

## Results

| Curve      | Field    | ℓ bits | h   | isInSubGroupNaive | ClearCofactor |
|------------|----------|--------|-----|-------------------|---------------|
| curve25519 | Fp, p255 | 253    | 8   | 66.2 µs           | 2.97 µs       |
| jubjub     | Fp, p255 | 252    | 8   | 62.8 µs           | 2.84 µs       |
| fourq      | Fp², p127| 246    | 392 | 102.0 µs          | 4.63 µs       |
| curve448   | Fp, p448 | 446    | 4   | 263.2 µs          | 29.5 µs       |
| gc256a     | Fp, p256 | 255    | 4   | 65.9 µs           | 2.17 µs       |

## Notes

- `isInSubGroupNaive` uses the generated `mulByOrder` addition chain
  (curve25519: 248 doublings + 34 additions; jubjub: 247 + 55;
  fourq: 234 + 51; curve448: 442 + 62; gc256a: 250 + 34), all in projective
  coordinates.
- `ClearCofactor` uses the generated `mulByCofactor` chain
  (h = 8: 3 doublings; h = 4: 2 doublings; h = 392: 8 doublings +
  2 additions) plus one projective→affine conversion, so a single field
  inversion dominates its cost. This is why curve448 (one 448-bit inversion)
  is an order of magnitude slower than the ~250-bit curves despite doing
  fewer doublings.
- Both operations run on public inputs in this context; the implementations
  are not constant-time.
