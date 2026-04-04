module github.com/yelhousni/divide-and-pair

go 1.25.7

tool (
	github.com/klauspost/asmfmt/cmd/asmfmt
	golang.org/x/tools/cmd/goimports
)

require (
	github.com/bits-and-blooms/bitset v1.24.4
	github.com/consensys/gnark-crypto v0.0.0-00010101000000-000000000000
	github.com/leanovate/gopter v0.2.11
	github.com/stretchr/testify v1.11.1
)

require (
	github.com/consensys/bavard v0.2.2-0.20260118153501-cba9f5475432 // indirect
	github.com/davecgh/go-spew v1.1.1 // indirect
	github.com/klauspost/asmfmt v1.3.2 // indirect
	github.com/mmcloughlin/addchain v0.4.0 // indirect
	github.com/pmezard/go-difflib v1.0.0 // indirect
	github.com/rogpeppe/go-internal v1.14.1 // indirect
	golang.org/x/mod v0.34.0 // indirect
	golang.org/x/sync v0.20.0 // indirect
	golang.org/x/sys v0.42.0 // indirect
	golang.org/x/telemetry v0.0.0-20260311193753-579e4da9a98c // indirect
	golang.org/x/tools v0.43.0 // indirect
	gopkg.in/yaml.v3 v3.0.1 // indirect
	rsc.io/tmplfunc v0.0.3 // indirect
)

replace github.com/consensys/gnark-crypto => /Users/youssefelhousni/workspace/consensys/curves/consensys:gnark-crypto/gnark-crypto
