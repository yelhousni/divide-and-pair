package fp

import (
	"math/big"
	"math/rand"
	"testing"

	"github.com/yelhousni/divide-and-pair/internal/octic"
)

func TestOcticSymbol(t *testing.T) {
	var generator, root Element
	generator.SetUint64(7)
	root.ExpByOcticExp(generator)
	var roots [8]Element
	roots[0].SetOne()
	for i := 1; i < len(roots); i++ {
		roots[i].Mul(&roots[i-1], &root)
	}
	if roots[4].IsOne() {
		t.Fatal("7 does not define a primitive eighth root")
	}
	// Verify the hard-coded prime, its orientation and the terminal unit law.
	if octicPrimeDivisor.norm.Cmp(Modulus()) != 0 {
		t.Fatal("incorrect prime norm")
	}
	var eval, coeff Element
	for i := 3; i >= 0; i-- {
		coeff.SetBigInt(&octicPrime[i])
		eval.Mul(&eval, &root).Add(&eval, &coeff)
	}
	if !eval.IsZero() {
		t.Fatal("incorrect embedding of the cyclotomic field")
	}
	var zero Element
	if zero.OcticSymbol() != -1 {
		t.Fatal("symbol of zero")
	}
	check := func(x Element) {
		t.Helper()
		var want Element
		want.ExpByOcticExp(x)
		j := x.OcticSymbol()
		if j < 0 || j >= 8 || !want.Equal(&roots[j]) {
			t.Fatalf("x=%s: symbol %d disagrees with Euler's criterion", x.String(), j)
		}
	}
	var x Element
	x.SetOne()
	for j := range roots {
		if got := x.OcticSymbol(); got != j {
			t.Fatalf("7^%d has symbol %d", j, got)
		}
		x.Mul(&x, &generator)
	}
	// Units need a final correction even after primary normalization.
	// Mix powers of epsilon = 1+t-t³ with the torsion units.
	var epsilon, t3, u, unit Element
	t3.Mul(&roots[2], &root)
	epsilon.SetOne().Add(&epsilon, &root).Sub(&epsilon, &t3)
	u.SetOne()
	for i := 0; i < 16; i++ {
		for j := range roots {
			unit.Mul(&u, &roots[j])
			check(unit)
		}
		u.Mul(&u, &epsilon)
	}
	for k := 0; k < 256; k++ {
		x.SetBigInt(new(big.Int).Lsh(big.NewInt(1), uint(k)))
		check(x)
		x.Neg(&x)
		check(x)
	}
	// A small raw representative forces the wide quotient guard.
	small := Element{3}
	if _, ok := octic.Descent(octicLift(&small), octicMod); ok {
		t.Fatal("expected large quotient for raw 3")
	}
	check(small)
	rng := rand.New(rand.NewSource(1))
	for i := 0; i < 1024; i++ {
		x.SetBigInt(new(big.Int).Rand(rng, Modulus()))
		if !x.IsZero() {
			check(x)
		}
	}
}

func BenchmarkOcticCharacter(b *testing.B) {
	var inputs [256]Element
	rng := rand.New(rand.NewSource(1))
	for i := range inputs {
		inputs[i].SetBigInt(new(big.Int).Rand(rng, Modulus()))
	}
	b.Run("Joye", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			inputs[i%len(inputs)].OcticSymbol()
		}
	})
	b.Run("JoyeBigInt", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			inputs[i%len(inputs)].octicSymbolReference()
		}
	})
	b.Run("Exp", func(b *testing.B) {
		var z Element
		for i := 0; i < b.N; i++ {
			z.ExpByOcticExp(inputs[i%len(inputs)])
		}
	})
}

func TestOcticLiftAndDescent(t *testing.T) {
	p := Modulus()
	toBig := func(w octic.Word) *big.Int {
		h := new(big.Int).Lsh(big.NewInt(w.Hi), 64)
		return h.Add(h, new(big.Int).SetUint64(w.Lo))
	}
	for i, c := range octicRecip {
		n := new(big.Int).Lsh(new(big.Int).Abs(&octicPrimeDivisor.adj[i]), 256)
		n.Add(n, &octicPrimeDivisor.half).Quo(n, p)
		got := new(big.Int)
		for j := 2; j >= 0; j-- {
			got.Lsh(got, 64).Add(got, new(big.Int).SetUint64(c[j]))
		}
		if got.Cmp(n) != 0 || (octicPrimeDivisor.adj[i].Sign() > 0) != (i == 1) || toBig(octicMod[i]).Cmp(&octicPrime[i]) != 0 {
			t.Fatalf("incorrect lift constants at coordinate %d", i)
		}
	}
	rng := rand.New(rand.NewSource(983))
	completed := 0
	for i := 0; i < 4096; i++ {
		var x Element
		x.SetBigInt(new(big.Int).Rand(rng, p))
		if x.IsZero() {
			continue
		}
		a := octicLift(&x)
		var eval, coeff Element
		for j := 3; j >= 0; j-- {
			c := toBig(a[j])
			if c.BitLen() > 65 {
				t.Fatalf("lift too large at input %d", i)
			}
			coeff.SetBigInt(c)
			eval.Mul(&eval, &octicRoots[1]).Add(&eval, &coeff)
		}
		var raw Element
		n := new(big.Int)
		for j := 3; j >= 0; j-- {
			n.Lsh(n, 64).Add(n, new(big.Int).SetUint64(x[j]))
		}
		raw.SetBigInt(n)
		if raw != eval {
			t.Fatalf("incorrect lift at input %d", i)
		}
		j, ok := octic.Descent(a, octicMod)
		if ok {
			completed++
		} else {
			j = x.OcticSymbol()
		}
		var expected Element
		expected.ExpByOcticExp(x)
		if expected != octicRoots[j] {
			t.Fatalf("incorrect symbol at input %d", i)
		}
		if i < 64 && j != x.octicSymbolReference() {
			t.Fatalf("disagreement with exact Joye at input %d", i)
		}
	}
	t.Logf("fixed-width descent completed %d/4096 random inputs", completed)
	if completed < 4000 {
		t.Fatal("frequent fallback on random inputs")
	}
}
