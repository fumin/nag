// Package field implements [finite field] arithmetic.
//
// [finite field]: https://en.wikipedia.org/wiki/Finite_field
package field

import (
	"math/big"
	"strconv"

	"github.com/fumin/nag"
)

// An IrreduciblePoly is an irreducible polynomial for the [construction] of a prime field extension.
//
// [construction]: https://en.wikipedia.org/wiki/Finite_field#Explicit_construction
type IrreduciblePoly struct {
	*nag.Polynomial[*prime]
}

// NewIrreduciblePoly returns an irreducible polynomial for the finite field GF(p^n), where p is a prime number and n >= 1.
func NewIrreduciblePoly(p, n int) *IrreduciblePoly {
	k := newPrime(p, 0)
	xn := nag.NewPolynomial(k, nag.Deglex, nag.PolynomialTerm[*prime]{Coefficient: k.NewOne(), Monomial: make([]nag.Symbol, n)})
	xn.SymbolStringer = func(nag.Symbol) string { return "x" }
	order := expi(p, n)
	for i := 1; i < order; i++ {
		poly := poly0(xn).Set(xn)
		ithPoly(poly, i, p)

		factors := factorize(poly)
		if !(len(factors) == 1 && factors[0].N == 1) {
			continue
		}
		// Only check primitivity if for small orders, since this involves iterating through every element in the field.
		if order < 1024 && !primitive(poly) {
			continue
		}
		return &IrreduciblePoly{poly}
	}
	panic("unable to find irreducible polynomial")
}

// Ext returns the i'th element in the finite field GF(p^n).
func (irr *IrreduciblePoly) Ext(i int) *PrimeExt {
	irrp := irr.Polynomial
	x := &PrimeExt{irr: poly0(irrp).Set(irrp), poly: poly0(irrp)}

	p, e := x.characteristic(), x.primePower()
	order := int(exp(nag.NewRat(int64(p), 1), e).Num().Int64())
	ith := i % order
	ithPoly(x.poly, ith, p)

	return x
}

// A PrimeExt is an element in the finite field GF(p^n), where p is a prime number and n >= 1.
type PrimeExt struct {
	irr  *nag.Polynomial[*prime]
	poly *nag.Polynomial[*prime]
}

// NewZero returns the additive identity 0.
func (x *PrimeExt) NewZero() *PrimeExt {
	y := &PrimeExt{}
	y.irr = poly0(x.irr).Set(x.irr)
	y.poly = poly0(y.irr)
	return y
}

// NewOne returns the multiplicative identity 1.
func (x *PrimeExt) NewOne() *PrimeExt {
	y := &PrimeExt{}
	y.irr = poly0(x.irr).Set(x.irr)
	y.poly = poly1(y.irr)
	return y
}

// Equal reports whether x and y are equal.
func (x *PrimeExt) Equal(y *PrimeExt) bool {
	if !x.irr.Equal(y.irr) {
		return false
	}
	if !x.poly.Equal(y.poly) {
		return false
	}
	return true
}

// Add sets z to the sum x+y and returns z.
func (z *PrimeExt) Add(x, y *PrimeExt) *PrimeExt {
	z.irr.Set(x.irr)
	z.poly.Add(x.poly, y.poly)
	_, remainder := divide(z.poly, z.irr)
	z.poly.Set(remainder)
	return z
}

// Sub sets z to the difference x-y and returns z.
func (z *PrimeExt) Sub(x, y *PrimeExt) *PrimeExt {
	z.irr.Set(x.irr)

	// Compute negative of y.poly.
	k := z.irr.Field()
	neg1 := k.NewZero()
	neg1.Sub(neg1, k.NewOne())
	negY := mulScalar(y.poly, neg1)

	z.poly.Add(x.poly, negY)
	_, remainder := divide(z.poly, z.irr)
	z.poly.Set(remainder)
	return z
}

// Mul sets z to the product x*y and returns z.
func (z *PrimeExt) Mul(x, y *PrimeExt) *PrimeExt {
	z.irr.Set(x.irr)

	xp := x.poly
	if xp == z.poly {
		xp = poly0(xp).Set(xp)
	}
	yp := y.poly
	if yp == z.poly {
		yp = poly0(yp).Set(yp)
	}
	z.poly.Mul(xp, yp)

	_, remainder := divide(z.poly, z.irr)
	z.poly.Set(remainder)
	return z
}

// Div sets z to the quotient x/y and returns z.
func (z *PrimeExt) Div(x, y *PrimeExt) *PrimeExt {
	yInv := z.NewZero().Inv(y)
	z.Mul(x, yInv)
	return z
}

// Inv sets z to 1/x and returns z.
func (z *PrimeExt) Inv(x *PrimeExt) *PrimeExt {
	z.irr.Set(x.irr)
	z.poly.Set(inverse(x.poly, z.irr))
	return z
}

// String returns the integer representation of x.
func (x *PrimeExt) String() string {
	p := x.characteristic()
	var ith int
	for c, w := range x.poly.Terms() {
		ith += int(c.i.Int64()) * expi(p, len(w))
	}
	return strconv.Itoa(ith)
}

func (x *PrimeExt) characteristic() int {
	return x.irr.LeadingTerm().Coefficient.characteristic()
}

func (x *PrimeExt) primePower() int {
	return len(x.irr.LeadingTerm().Monomial)
}

type prime struct {
	// Order is the number of elements in the finite field which must be a prime number.
	order *big.Int
	// I is the integer representation of the element.
	i *big.Int
}

func newPrime(p, i int) *prime {
	return &prime{order: big.NewInt(int64(p)), i: big.NewInt(int64(i))}
}

func (x *prime) NewZero() *prime {
	return &prime{order: big.NewInt(0).Set(x.order), i: big.NewInt(0)}
}

func (x *prime) NewOne() *prime {
	return &prime{order: big.NewInt(0).Set(x.order), i: big.NewInt(1)}
}

func (x *prime) Equal(y *prime) bool {
	return (x.order.Cmp(y.order) == 0) && (x.i.Cmp(y.i) == 0)
}

func (z *prime) Add(x, y *prime) *prime {
	z.order.Set(x.order)
	z.i.Add(x.i, y.i)
	z.i.Mod(z.i, z.order)
	return z
}

func (z *prime) Sub(x, y *prime) *prime {
	z.order.Set(x.order)
	z.i.Sub(x.i, y.i)
	z.i.Mod(z.i, z.order)
	return z
}

func (z *prime) Mul(x, y *prime) *prime {
	z.order.Set(x.order)
	z.i.Mul(x.i, y.i)
	z.i.Mod(z.i, z.order)
	return z
}

func (z *prime) Div(x, y *prime) *prime {
	z.Inv(y)
	z.i.Mul(x.i, z.i)
	z.i.Mod(z.i, z.order)
	return z
}

func (z *prime) Inv(x *prime) *prime {
	if x.i.Sign() == 0 {
		panic("division by zero")
	}
	z.order.Set(x.order)
	z.i.ModInverse(x.i, z.order)
	return z
}

func (x *prime) String() string {
	return x.i.String()
}

func (x *prime) characteristic() int {
	return int(x.order.Int64())
}

func (x *prime) primePower() int { return 1 }
