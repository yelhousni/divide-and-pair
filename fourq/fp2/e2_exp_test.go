package fp2

import (
	"math/big"
	"testing"

	fp "github.com/yelhousni/divide-and-pair/fourq/fp"
)

func TestExpBySeptic(t *testing.T) {
	p := fp.Modulus()
	exp := new(big.Int).Sub(p, big.NewInt(1))
	exp.Div(exp, big.NewInt(7))

	var x E2
	x.A0.SetUint64(42)
	x.A1.SetUint64(137)

	var got, want E2
	got.ExpBySeptic(&x)
	want.Exp(&x, exp)
	if !got.Equal(&want) {
		t.Fatal("ExpBySeptic does not match Exp")
	}
}

func BenchmarkE2ExpNaive(b *testing.B) {
	p := fp.Modulus()
	exp := new(big.Int).Sub(p, big.NewInt(1))
	exp.Div(exp, big.NewInt(7))

	var x, z E2
	x.A0.SetUint64(42)
	x.A1.SetUint64(137)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		z.Exp(&x, exp)
	}
}

func BenchmarkE2ExpAddChain(b *testing.B) {
	var x, z E2
	x.A0.SetUint64(42)
	x.A1.SetUint64(137)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		z.ExpBySeptic(&x)
	}
}
