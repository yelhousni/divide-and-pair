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

## Quartic characters over Fp²

Medians of seven 350 ms Go benchmarks on Apple M5, Go 1.26.0 darwin/arm64,
with `purego` and `-cpu=1`. Each field comparison uses the same 256 seeded
extension-field elements. Each subgroup comparison uses the same 256 seeded
subgroup points. The latter produces its own distribution of Miller values.

| Curve | Operation | Torus PRAC | Gaussian GCD |
| --- | --- | ---: | ---: |
| Curve448 | Field operation | 48.248 µs | 16.850 µs |
| Curve448 | Subgroup membership | 49.152 µs | 17.548 µs |
| GC256A | Field operation | 13.412 µs | 7.707 µs |
| GC256A | Subgroup membership | 11.950 µs | 7.944 µs |

The torus column tests whether the character equals one. The GCD computes
the full index in 0–3 and is the default for subgroup membership on both
curves. Both GCD paths allocated zero bytes per call. The torus path allocated
zero bytes on Curve448 and 600 bytes in 12 allocations per call on GC256A.

```sh
go test -p 1 -cpu=1 -tags purego -run '^$' -bench '^(BenchmarkQuarticExtension|BenchmarkIsInSubGroupQuartic(GCD)?)$' -benchmem -benchtime=350ms -count=7 ./curve448 ./gc256a
```

## FourQ octic character over Fp²

Medians of seven 350 ms Go benchmarks on Apple M5, Go 1.26.0 darwin/arm64,
with `purego` and `-cpu=1`. The field benchmarks cycle over the same 256
seeded extension-field elements. The subgroup benchmarks use the same 256
seeded subgroup points and perform both septic checks.

| Operation | Existing torus | Joye |
| --- | ---: | ---: |
| Field operation | 1.558 µs | 2.256 µs |
| Subgroup membership | 4.779 µs | 5.711 µs |

The torus column tests whether the character equals one. Joye computes the
full index in 0–7. Both methods allocated zero bytes per call. The default
subgroup test uses the torus.

```sh
go test -p 1 -cpu=1 -tags purego -run '^$' -bench '^(BenchmarkOcticExtension|BenchmarkIsInSubGroupTate(Joye)?)$' -benchmem -benchtime=350ms -count=7 ./fourq
```
