//go:build !purego

package fp

var supportAdx = true

//go:noescape
func MulBy3(x *Element)

//go:noescape
func MulBy5(x *Element)

//go:noescape
func MulBy13(x *Element)

//go:noescape
func mul(res, x, y *Element)

//go:noescape
func fromMont(res *Element)

//go:noescape
func reduce(res *Element)

//go:noescape
func Butterfly(a, b *Element)

func (z *Element) Mul(x, y *Element) *Element {
	mul(z, x, y)
	return z
}

func (z *Element) Square(x *Element) *Element {
	mul(z, x, x)
	return z
}

