// Package field implements [finite field] arithmetic.
//
// [finite field]: https://en.wikipedia.org/wiki/Finite_field
package field

import (
	"math/big"

	"github.com/fumin/nag"
)

// A Finite is an element of a finite field GF(q) with order q = p^n where p is a prime number.
type Finite[K nag.Field[K]] interface {
	nag.Field[K]
	// Characteristic returns p in the finite field GF(p^n).
	Characteristic() *big.Int
	// PrimePower returns n in the finite field GF(p^n).
	PrimePower() *big.Int
	// Ith returns the i'th element in the field.
	Ith(i *big.Int) K
}

// An IrreduciblePoly is an irreducible polynomial for the [construction] of a prime field extension.
//
// [construction]: https://en.wikipedia.org/wiki/Finite_field#Explicit_construction
type IrreduciblePoly struct {
	*nag.Polynomial[*prime]
}

// NewIrreduciblePoly returns an irreducible polynomial for the finite field GF(p^n), where p is a prime number and n >= 1.
func NewIrreduciblePoly(p *big.Int, n int) *IrreduciblePoly {
	k := &prime{order: new(big.Int).Set(p), i: big.NewInt(0)}
	xn := nag.NewPolynomial(k, nag.Deglex, nag.PolynomialTerm[*prime]{Coefficient: k.NewOne(), Monomial: make([]nag.Symbol, n)})
	xn.SymbolStringer = func(nag.Symbol) string { return "x" }
	order := new(big.Int).Exp(p, big.NewInt(int64(n)), nil)
	for i := big.NewInt(1); i.Cmp(order) < 0; i.Add(i, big.NewInt(1)) {
		poly := poly0(xn).Set(xn)
		ithPoly(poly, i, p)

		factors := Factor(poly)
		if !(len(factors) == 1 && factors[0].N.Cmp(big.NewInt(1)) == 0) {
			continue
		}
		// Only check primitivity if for small orders, since this involves iterating through every element in the field.
		if order.Cmp(big.NewInt(1024)) < 0 && !primitive(poly) {
			continue
		}
		return &IrreduciblePoly{poly}
	}
	panic("unable to find irreducible polynomial")
}

// Ext returns the i'th element in the finite field GF(p^n).
func (irr *IrreduciblePoly) Ext(i *big.Int) *PrimeExt {
	irrp := irr.Polynomial
	x := &PrimeExt{irr: &IrreduciblePoly{poly0(irrp).Set(irrp)}, poly: poly0(irrp)}

	order := new(big.Int).Exp(x.Characteristic(), x.PrimePower(), nil)
	ith := new(big.Int).Mod(i, order)
	ithPoly(x.poly, ith, x.Characteristic())

	return x
}

// A PrimeExt is an element in the finite field GF(p^n), where p is a prime number and n >= 1.
type PrimeExt struct {
	irr  *IrreduciblePoly
	poly *nag.Polynomial[*prime]
}

// NewZero returns the additive identity 0.
func (x *PrimeExt) NewZero() *PrimeExt {
	y := &PrimeExt{}
	y.irr = &IrreduciblePoly{poly0(x.irr.Polynomial).Set(x.irr.Polynomial)}
	y.poly = poly0(y.irr.Polynomial)
	return y
}

// NewOne returns the multiplicative identity 1.
func (x *PrimeExt) NewOne() *PrimeExt {
	y := &PrimeExt{}
	y.irr = &IrreduciblePoly{poly0(x.irr.Polynomial).Set(x.irr.Polynomial)}
	y.poly = poly1(y.irr.Polynomial)
	return y
}

// Equal reports whether x and y are equal.
func (x *PrimeExt) Equal(y *PrimeExt) bool {
	if !x.irr.Polynomial.Equal(y.irr.Polynomial) {
		return false
	}
	if !x.poly.Equal(y.poly) {
		return false
	}
	return true
}

// Add sets z to the sum x+y and returns z.
func (z *PrimeExt) Add(x, y *PrimeExt) *PrimeExt {
	z.irr.Polynomial.Set(x.irr.Polynomial)
	z.poly.Add(x.poly, y.poly)
	_, remainder := divide(z.poly, z.irr.Polynomial)
	z.poly.Set(remainder)
	return z
}

// Sub sets z to the difference x-y and returns z.
func (z *PrimeExt) Sub(x, y *PrimeExt) *PrimeExt {
	z.irr.Polynomial.Set(x.irr.Polynomial)

	// Compute negative of y.poly.
	k := z.irr.Field()
	neg1 := k.NewZero()
	neg1.Sub(neg1, k.NewOne())
	negY := mulScalar(y.poly, neg1)

	z.poly.Add(x.poly, negY)
	_, remainder := divide(z.poly, z.irr.Polynomial)
	z.poly.Set(remainder)
	return z
}

// Mul sets z to the product x*y and returns z.
func (z *PrimeExt) Mul(x, y *PrimeExt) *PrimeExt {
	z.irr.Polynomial.Set(x.irr.Polynomial)

	xp := x.poly
	if xp == z.poly {
		xp = poly0(xp).Set(xp)
	}
	yp := y.poly
	if yp == z.poly {
		yp = poly0(yp).Set(yp)
	}
	z.poly.Mul(xp, yp)

	_, remainder := divide(z.poly, z.irr.Polynomial)
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
	z.irr.Polynomial.Set(x.irr.Polynomial)
	z.poly.Set(inverse(x.poly, z.irr.Polynomial))
	return z
}

// String returns the integer representation of x.
func (x *PrimeExt) String() string {
	p := x.Characteristic()
	ith := big.NewInt(0)
	for c, w := range x.poly.Terms() {
		incr := new(big.Int).Mul(c.i, new(big.Int).Exp(p, big.NewInt(int64(len(w))), nil))
		ith.Add(ith, incr)
	}
	return ith.String()
}

func (x *PrimeExt) Characteristic() *big.Int {
	return x.irr.LeadingTerm().Coefficient.Characteristic()
}

func (x *PrimeExt) PrimePower() *big.Int {
	return big.NewInt(int64(len(x.irr.LeadingTerm().Monomial)))
}

func (x *PrimeExt) Ith(i *big.Int) *PrimeExt {
	return x.irr.Ext(i)
}

type prime struct {
	// Order is the number of elements in the finite field which must be a prime number.
	order *big.Int
	// I is the integer representation of the element.
	i *big.Int
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

func (x *prime) Characteristic() *big.Int {
	return big.NewInt(0).Set(x.order)
}

func (x *prime) PrimePower() *big.Int { return big.NewInt(1) }

func (x *prime) Ith(i *big.Int) *prime {
	return &prime{order: new(big.Int).Set(x.order), i: new(big.Int).Set(i)}
}
