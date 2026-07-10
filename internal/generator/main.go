// Command generator splices, into each curve package's point.go, fixed-scalar
// multiplication methods on PointProj — mulByOrder (scalar = the prime
// subgroup order ℓ) and mulByCofactor (scalar = the cofactor h) — unrolled as
// short addition chains found with github.com/mmcloughlin/addchain. Each
// method lives between BEGIN/END generation markers and is replaced in place
// on regeneration.
//
// The methods are written both to the main curve packages and to their
// anonymized copies under artifact/. Run from anywhere:
//
//	go run -C internal/generator .
//
// or via `go generate ./...` from the repository root.
package main

import (
	"fmt"
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

	// Cofactor is the curve cofactor h in decimal.
	Cofactor string
}

var curves = []curveConfig{
	{
		Package:      "curve25519",
		Cofactor:     "8",
		Order:        "7237005577332262213973186563042994240857116359379907606001950938285454250989",
		OrderComment: "2^252 + 27742317777372353535851937790883648493",
	},
	{
		Package:  "jubjub",
		Cofactor: "8",
		Order:    "6554484396890773809930967563523245729705921265872317281365359162392183254199",
	},
	{
		Package:  "fourq",
		Cofactor: "392",
		Order:    "73846995687063900142583536357581573884798075859800097461294096333596429543",
	},
	{
		Package:      "curve448",
		Cofactor:     "4",
		Order:        "181709681073901722637330951972001133588410340171829515070372549795146003961539585716195755291692375963310293709091662304773755859649779",
		OrderComment: "2^446 - 13818066809895115352007386748515426880336692474882178609894547503885",
	},
	{
		Package:  "gc256a",
		Cofactor: "4",
		Order:    "28948022309329048855892746252171976963338560298092253442512153408785530358887",
	},
}

func main() {
	log.SetFlags(0)
	log.SetPrefix("generator: ")
	root := repoRoot()

	for _, c := range curves {
		order := mustInt(c.Package, "order", c.Order)
		cofactor := mustInt(c.Package, "cofactor", c.Cofactor)

		methods := []struct {
			spec methodSpec
			n    *big.Int
		}{
			{methodSpec{Name: "mulByOrder", Doc: mulByOrderDoc(c, order)}, order},
			{methodSpec{Name: "mulByCofactor", Doc: mulByCofactorDoc(c)}, cofactor},
		}

		blocks := make(map[string]string, len(methods))
		for _, m := range methods {
			log.Printf("%s: searching addition chain for %s (%d-bit scalar)",
				c.Package, m.spec.Name, m.n.BitLen())
			prog, script := searchChain(m.n)
			block, err := emitMethod(m.spec, prog, script)
			if err != nil {
				log.Fatalf("%s: %s: %v", c.Package, m.spec.Name, err)
			}
			blocks[m.spec.Name] = block
			log.Printf("%s: %s: %d doublings, %d additions",
				c.Package, m.spec.Name, prog.Program.Doubles(), prog.Program.Adds())
		}

		for _, dir := range []string{c.Package, filepath.Join("artifact", c.Package)} {
			path := filepath.Join(root, dir, "point.go")
			src, err := os.ReadFile(path)
			if err != nil {
				log.Fatal(err)
			}
			for _, m := range methods {
				src, err = splice(src, m.spec.Name, blocks[m.spec.Name])
				if err != nil {
					log.Fatalf("%s: %v", path, err)
				}
			}
			if err := os.WriteFile(path, src, 0o644); err != nil {
				log.Fatal(err)
			}
			log.Printf("%s: updated %s", c.Package, path)
		}
	}
}

// mustInt parses a decimal integer from the curve config.
func mustInt(pkg, what, s string) *big.Int {
	n, ok := new(big.Int).SetString(s, 10)
	if !ok {
		log.Fatalf("%s: invalid %s %q", pkg, what, s)
	}
	return n
}

// mulByOrderDoc returns the leading doc comment of the mulByOrder method.
func mulByOrderDoc(c curveConfig, n *big.Int) []string {
	doc := []string{
		"// mulByOrder computes p = [ℓ]a and returns p, where",
		"//",
		"//\tℓ = " + c.Order,
	}
	if c.OrderComment != "" {
		doc = append(doc, "//\t  = "+c.OrderComment)
	}
	return append(doc,
		"//",
		fmt.Sprintf("// is the %d-bit prime subgroup order. The scalar multiplication is unrolled", n.BitLen()),
		fmt.Sprintf("// as a short addition chain generated with %s.", addchainCitation()),
		"//",
		"// The unified Add/Double formulas are valid for arbitrary curve points",
		"// (including the identity and points outside the ℓ-torsion subgroup), so the",
		"// result is the identity if and only if a is in the prime-order subgroup.",
	)
}

// mulByCofactorDoc returns the leading doc comment of the mulByCofactor method.
func mulByCofactorDoc(c curveConfig) []string {
	return []string{
		"// mulByCofactor computes p = [h]a and returns p, where",
		"//",
		"//\th = " + c.Cofactor,
		"//",
		"// is the cofactor. The result lies in the prime-order subgroup for any",
		fmt.Sprintf("// curve point a. The chain was generated with %s.", addchainCitation()),
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
