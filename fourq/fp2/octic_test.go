package fp2

import (
	"math/big"
	"math/rand"
	"testing"

	"github.com/yelhousni/divide-and-pair/fourq/fp"
	"github.com/yelhousni/divide-and-pair/internal/octic"
)

func TestOcticSymbol(t *testing.T) {
	p := fp.Modulus()
	exp := new(big.Int).Mul(p, p)
	exp.Sub(exp, big.NewInt(1)).Rsh(exp, 3)
	var zero E2
	if zero.OcticSymbol() != -1 {
		t.Fatal("zero")
	}
	var seen [8]bool
	check := func(x E2) {
		t.Helper()
		var want E2
		want.Exp(&x, exp)
		j := x.OcticSymbol()
		if j < 0 || j >= 8 || octicRoots[j] != want || x.octicSymbolExp() != j {
			t.Fatalf("x=%s: index %d, Euler %s", x.String(), j, want.String())
		}
		seen[j] = true
	}
	for _, r := range octicRoots {
		check(r)
	}
	for k := 0; k < 127; k++ {
		var x E2
		x.A0.SetBigInt(new(big.Int).Lsh(big.NewInt(1), uint(k)))
		x.A1.SetOne()
		check(x)
		x.A1.Neg(&x.A1)
		check(x)
	}
	// Raw 3 produces a quotient outside the signed descent bound.
	small := E2{A0: fp.Element{3}}
	if _, ok := octic.Descent(octicLift(&small), octicPrime); ok {
		t.Fatal("expected the large-quotient fallback")
	}
	check(small)
	rng := rand.New(rand.NewSource(8))
	completed := 0
	for range 4096 {
		var x E2
		x.A0.SetBigInt(new(big.Int).Rand(rng, p))
		x.A1.SetBigInt(new(big.Int).Rand(rng, p))
		check(x)
		if _, ok := octic.Descent(octicLift(&x), octicPrime); ok {
			completed++
		}
	}
	t.Logf("Joye completed %d/4096 random inputs", completed)
	if completed < 4000 {
		t.Fatal("frequent fallback")
	}
	for j, ok := range seen {
		if !ok {
			t.Fatalf("missed index %d", j)
		}
	}
}
