package fp

import (
	"math/big"
	"math/rand"
	"testing"
)

func septicTestInput(rng *rand.Rand) Element {
	x := Element{rng.Uint64(), rng.Uint64() >> 1}
	if x == qElement {
		return Element{}
	}
	return x
}

func TestSepticPaperExample(t *testing.T) {
	// Section 8, N(beta)=463, independently checked in the paper by Euler's criterion.
	a := septicInteger{0, 128, 80, 44, 161, -21, 189}
	b := septicInteger{0, -3, -2, 0, 0, 0, 0}
	if got, ok := septicDescent(a, b, b.primary()); !ok || got != 6 {
		t.Fatalf("paper example: got (%d,%t), want 6", got, ok)
	}
}

func TestSepticLiftConstants(t *testing.T) {
	// Since 2^128=2p+2, rounding gives exactly C=2p/pi for this pi.
	// Check the identity in the full cyclotomic ring, independently of
	// the production low-word multiplication.
	var product [7]big.Int
	var magnitude, term big.Int
	for i, c := range septicRecip {
		magnitude.SetUint64(c.hi).Lsh(&magnitude, 64)
		magnitude.Add(&magnitude, term.SetUint64(c.lo))
		if magnitude.BitLen() > 125 {
			t.Fatal("reciprocal too large for rounding identity")
		}
		if c.negative {
			magnitude.Neg(&magnitude)
		}
		for j := 1; j <= 6; j++ {
			k := (i + 1 + j) % 7
			term.SetInt64(septicPrime[j]).Mul(&term, &magnitude)
			product[k].Add(&product[k], &term)
		}
	}
	want := new(big.Int).Lsh(Modulus(), 1)
	for i := 1; i < 7; i++ {
		if term.Sub(&product[0], &product[i]).Cmp(want) != 0 {
			t.Fatal("incorrect prime or reciprocal constant")
		}
	}
}

func TestSepticSymbol(t *testing.T) {
	var root Element
	root.ExpBySepticFp(NewElement(9))
	if root != septicRoots[1] || root.IsOne() {
		t.Fatal("incorrect root")
	}
	var product Element
	product.Mul(&septicRoots[6], &root)
	if !product.IsOne() {
		t.Fatal("root has incorrect order")
	}
	eval := func(a septicInteger) Element {
		var v, c Element
		for i := 6; i >= 0; i-- {
			c.SetInt64(a[i])
			v.Mul(&v, &root).Add(&v, &c)
		}
		return v
	}
	if eval(septicPrime) != (Element{}) {
		t.Fatal("incorrect cyclotomic prime")
	}
	var zero Element
	if zero.SepticSymbol() != -1 {
		t.Fatal("incorrect symbol of zero")
	}
	var counts [7]int
	completed := 0
	check := func(x Element) {
		t.Helper()
		lift := septicLift(x.Bits())
		if lift.width() > 24 || eval(lift) != x {
			t.Fatalf("incorrect lift for %s", x.String())
		}
		j, ok := septicDescent(lift, septicPrime, septicPrimePrimary)
		got := x.SepticSymbol()
		if ok {
			completed++
			if j != got {
				t.Fatalf("public result %d differs from descent %d", got, j)
			}
		}
		var want Element
		want.ExpBySepticFp(x)
		if got < 0 || got >= 7 || septicRoots[got] != want {
			t.Fatalf("x=%s: got (%d,%t), want %s", x.String(), got, ok, want.String())
		}
		counts[got]++
	}
	rng := rand.New(rand.NewSource(7))
	for range 4096 {
		check(septicTestInput(rng))
	}
	t.Logf("descent completed %d/4096 random inputs", completed)
	if completed < 4000 {
		t.Fatal("frequent fallback on random inputs")
	}
	for j, n := range counts {
		if n == 0 {
			t.Errorf("missed index %d", j)
		}
	}
	for i := 1; i <= 128; i++ {
		x := NewElement(uint64(i))
		check(x)
		x.Neg(&x)
		check(x)
	}
	x := NewElement(1)
	for range 127 {
		check(x)
		x.Double(&x)
	}
	for _, x := range septicRoots {
		check(x)
	}
}

func BenchmarkSepticSymbol(b *testing.B) {
	var inputs [256]Element
	rng := rand.New(rand.NewSource(7))
	for i := range inputs {
		inputs[i] = septicTestInput(rng)
	}
	b.Run("Euclidean", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			inputs[i%len(inputs)].SepticSymbol()
		}
	})
	b.Run("Exp", func(b *testing.B) {
		var z Element
		for i := 0; i < b.N; i++ {
			z.ExpBySepticFp(inputs[i%len(inputs)])
		}
	})
}
