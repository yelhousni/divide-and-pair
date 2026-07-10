// Command generator emits, for each curve package, a mul_by_order.go file
// containing a mulByOrder method on PointProj: a scalar multiplication by the
// curve's prime subgroup order ℓ, unrolled as a short addition chain found
// with github.com/mmcloughlin/addchain.
//
// The generated files are written both to the main curve packages and to
// their anonymized copies under artifact/. Run from anywhere:
//
//	go run -C internal/generator .
//
// or via `go generate ./...` from the repository root.
package main

import (
	"log"
	"math/big"
	"os"
	"path/filepath"
	"runtime"

	"github.com/mmcloughlin/addchain/acc"
	"github.com/mmcloughlin/addchain/acc/ir"
	"github.com/mmcloughlin/addchain/acc/pass"
	"github.com/mmcloughlin/addchain/acc/printer"
	"github.com/mmcloughlin/addchain/alg/ensemble"
	"github.com/mmcloughlin/addchain/alg/exec"
)

type curveConfig struct {
	// Package is the curve package name; the file is written to <root>/<Package>/
	// and <root>/artifact/<Package>/.
	Package string

	// Order is the prime subgroup order ℓ in decimal, copied verbatim from the
	// package's curve.go.
	Order string

	// OrderComment optionally gives a human-readable special form of the order.
	OrderComment string
}

var curves = []curveConfig{
	{
		Package:      "curve25519",
		Order:        "7237005577332262213973186563042994240857116359379907606001950938285454250989",
		OrderComment: "2^252 + 27742317777372353535851937790883648493",
	},
	{
		Package: "jubjub",
		Order:   "6554484396890773809930967563523245729705921265872317281365359162392183254199",
	},
	{
		Package: "fourq",
		Order:   "73846995687063900142583536357581573884798075859800097461294096333596429543",
	},
	{
		Package:      "curve448",
		Order:        "181709681073901722637330951972001133588410340171829515070372549795146003961539585716195755291692375963310293709091662304773755859649779",
		OrderComment: "2^446 - 13818066809895115352007386748515426880336692474882178609894547503885",
	},
	{
		Package: "gc256a",
		Order:   "28948022309329048855892746252171976963338560298092253442512153408785530358887",
	},
}

func main() {
	log.SetFlags(0)
	log.SetPrefix("generator: ")
	root := repoRoot()

	for _, c := range curves {
		n, ok := new(big.Int).SetString(c.Order, 10)
		if !ok {
			log.Fatalf("%s: invalid order %q", c.Package, c.Order)
		}

		log.Printf("%s: searching addition chain for %d-bit order", c.Package, n.BitLen())
		prog, script := searchChain(n)

		block, err := emitMethod(c, n, prog, script)
		if err != nil {
			log.Fatalf("%s: %v", c.Package, err)
		}

		for _, dir := range []string{c.Package, filepath.Join("artifact", c.Package)} {
			path := filepath.Join(root, dir, "point.go")
			src, err := os.ReadFile(path)
			if err != nil {
				log.Fatal(err)
			}
			out, err := splice(src, block)
			if err != nil {
				log.Fatalf("%s: %v", path, err)
			}
			if err := os.WriteFile(path, out, 0o644); err != nil {
				log.Fatal(err)
			}
			// Clean up the previous generation layout, which kept the method
			// in its own file.
			_ = os.Remove(filepath.Join(root, dir, "mul_by_order.go"))
			log.Printf("%s: updated %s (%d doublings, %d additions)",
				c.Package, path, prog.Program.Doubles(), prog.Program.Adds())
		}
	}
}

// searchChain runs the default addchain algorithm ensemble on n, picks the
// shortest resulting program and returns its allocated intermediate
// representation together with the chain script.
func searchChain(n *big.Int) (*ir.Program, string) {
	algorithms := ensemble.Ensemble()
	ex := exec.NewParallel()
	results := ex.Execute(n, algorithms)

	best := 0
	for i, r := range results {
		if r.Err != nil {
			log.Fatal(r.Err)
		}
		if len(results[i].Program) < len(results[best].Program) {
			best = i
		}
	}

	// Round-trip through the acc script representation, then allocate
	// temporaries and evaluate the chain (same pipeline as gnark-crypto's
	// internal/generator/addchain).
	p, err := acc.Decompile(results[best].Program)
	if err != nil {
		log.Fatal(err)
	}
	script, err := acc.Build(p)
	if err != nil {
		log.Fatal(err)
	}
	prog, err := acc.Translate(script)
	if err != nil {
		log.Fatal(err)
	}
	allocator := pass.Allocator{Input: "x", Output: "z", Format: "t%d"}
	if err := pass.Exec(prog, allocator, pass.Func(pass.Eval)); err != nil {
		log.Fatal(err)
	}

	// Sanity checks: the chain must evaluate to n and end in the output
	// variable z.
	if last := prog.Chain[len(prog.Chain)-1]; last.Cmp(n) != 0 {
		log.Fatalf("chain evaluates to %s, want %s", last, n)
	}
	if out := prog.Instructions[len(prog.Instructions)-1].Output.String(); out != "z" {
		log.Fatalf("chain output variable is %q, want \"z\"", out)
	}

	scriptStr, err := printer.String(script)
	if err != nil {
		log.Fatal(err)
	}
	return prog, scriptStr
}

// repoRoot resolves the repository root from this source file's location so
// that output paths are independent of the working directory.
func repoRoot() string {
	_, thisFile, _, ok := runtime.Caller(0)
	if !ok {
		log.Fatal("unable to locate generator source via runtime.Caller")
	}
	return filepath.Join(filepath.Dir(thisFile), "..", "..")
}
