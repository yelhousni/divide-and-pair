package fp2

// ExpBySqrtPp1o4 raises z to the (p+1)/4 power in Fp2.
// (p+1)/4 = 0x3fffffffffffffffffffffffffffffffffffffffffffffffffffffffc0000000000000000000000000000000000000000000000000000000
//         = 2^446 - 2^222 = 2^222 * (2^224 - 1)
//
// Uses the same addition chain as fp.ExpBySqrtPp1o4:
// 445 squarings and 9 multiplications (in Fp2).
func (z *E2) ExpBySqrtPp1o4(x *E2) *E2 {
	// addition chain:
	//
	//	_10       = 2*1
	//	_11       = 1 + _10
	//	_1100     = _11 << 2
	//	_1111     = _11 + _1100
	//	_11110000 = _1111 << 4
	//	_11111111 = _1111 + _11110000
	//	x16       = _11111111 << 8 + _11111111
	//	x32       = x16 << 16 + x16
	//	x64       = x32 << 32 + x32
	//	x128      = x64 << 64 + x64
	//	x192      = x128 << 64 + x64
	//	x224      = x192 << 32 + x32
	//	return      x224 << 222
	//
	// Operations: 445 squares 9 multiplies
	var t0, t1 E2

	// Step 1: z = x^0x2
	z.Square(x)

	// Step 2: z = x^0x3
	z.Mul(x, z)

	// Step 4: t0 = x^0xc
	t0.Square(z)
	for s := 1; s < 2; s++ {
		t0.Square(&t0)
	}

	// Step 5: z = x^0xf
	z.Mul(z, &t0)

	// Step 9: t0 = x^0xf0
	t0.Square(z)
	for s := 1; s < 4; s++ {
		t0.Square(&t0)
	}

	// Step 10: z = x^0xff
	z.Mul(z, &t0)

	// Step 18: t0 = x^0xff00
	t0.Square(z)
	for s := 1; s < 8; s++ {
		t0.Square(&t0)
	}

	// Step 19: z = x^0xffff
	z.Mul(z, &t0)

	// Step 35: t0 = x^0xffff0000
	t0.Square(z)
	for s := 1; s < 16; s++ {
		t0.Square(&t0)
	}

	// Step 36: z = x^0xffffffff
	z.Mul(z, &t0)

	// Step 68: t0 = x^0xffffffff00000000
	t0.Square(z)
	for s := 1; s < 32; s++ {
		t0.Square(&t0)
	}

	// Step 69: t0 = x^0xffffffffffffffff
	t0.Mul(z, &t0)

	// Step 133: t1 = x^0xffffffffffffffff0000000000000000
	t1.Square(&t0)
	for s := 1; s < 64; s++ {
		t1.Square(&t1)
	}

	// Step 134: t1 = x^0xffffffffffffffffffffffffffffffff
	t1.Mul(&t0, &t1)

	// Step 198: t1 = x^0xffffffffffffffffffffffffffffffff0000000000000000
	for range 64 {
		t1.Square(&t1)
	}

	// Step 199: t0 = x^0xffffffffffffffffffffffffffffffffffffffffffffffff
	t0.Mul(&t0, &t1)

	// Step 231: t0 = x^0xffffffffffffffffffffffffffffffffffffffffffffffff00000000
	for range 32 {
		t0.Square(&t0)
	}

	// Step 232: z = x^0xffffffffffffffffffffffffffffffffffffffffffffffffffffffff
	z.Mul(z, &t0)

	// Step 454: z = x^0x3fffffffffffffffffffffffffffffffffffffffffffffffffffffffc0000000000000000000000000000000000000000000000000000000
	for range 222 {
		z.Square(z)
	}

	return z
}
