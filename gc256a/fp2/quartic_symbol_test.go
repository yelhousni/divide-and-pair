package fp2

import (
	"math/big"
	"math/rand"
	"testing"

	"github.com/yelhousni/divide-and-pair/gc256a/fp"
)

func TestQuarticSymbol(t *testing.T) {
	p := fp.Modulus()
	exp := new(big.Int).Mul(p, p)
	exp.Sub(exp, big.NewInt(1)).Rsh(exp, 2)
	var roots [4]E2
	roots[0].SetOne()
	roots[1].A1.SetOne()
	roots[2].Neg(&roots[0])
	roots[3].Neg(&roots[1])
	var zero E2
	if zero.QuarticSymbol() != -1 {
		t.Fatal("zero")
	}
	var seen [4]bool
	check := func(x E2) {
		t.Helper()
		var want E2
		want.Exp(&x, exp)
		j := x.QuarticSymbol()
		if j < 0 || j >= 4 || roots[j] != want {
			t.Fatalf("x=%s: index %d, Euler %s", x.String(), j, want.String())
		}
		seen[j] = true
	}
	for _, x := range roots {
		check(x)
	}
	for k := 1; k < fp.Bits; k += 31 {
		var x E2
		x.A0.SetBigInt(new(big.Int).Lsh(big.NewInt(1), uint(k)))
		x.A1.SetOne()
		check(x)
		x.A1.Neg(&x.A1)
		check(x)
	}
	// Raw representatives with several all-zero low limbs exercise batched
	// removal of powers of 1+i at, immediately below and above word boundaries.
	for k := 1; k < fp.Limbs; k++ {
		for _, offset := range []int{-1, 0, 1} {
			bit := 64*k + offset
			var x E2
			x.A0[bit/64] = uint64(1) << uint(bit%64)
			x.A1 = x.A0
			check(x)
			x.A0.SetZero()
			check(x)
		}
	}
	rng := rand.New(rand.NewSource(4))
	for range 2048 {
		var x E2
		x.A0.SetBigInt(new(big.Int).Rand(rng, p))
		x.A1.SetBigInt(new(big.Int).Rand(rng, p))
		check(x)
	}
	for j, ok := range seen {
		if !ok {
			t.Fatalf("missed index %d", j)
		}
	}
}
