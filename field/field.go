// Package field implements [finite field] arithmetic.
//
// [finite field]: https://en.wikipedia.org/wiki/Finite_field
package field

import (
	"fmt"
	"math/big"

	"github.com/fumin/nag"
)

// An Extension is a field extension [constructed] from irreducible polynomials.
//
// [constructed]: https://en.wikipedia.org/wiki/Finite_field#Explicit_construction
type Extension[K nag.Field[K]] struct {
	// Irr is the irreducible polynomial modulus.
	Irr *nag.Polynomial[K]
	// Poly is the polynomial representing the element.
	Poly *nag.Polynomial[K]
}

// NewZero returns the additive identity 0.
func (x *Extension[K]) NewZero() *Extension[K] {
	y := &Extension[K]{}
	y.Irr = poly0(x.Irr).Set(x.Irr)
	y.Poly = poly0(y.Irr)
	return y
}

// NewOne returns the multiplicative identity 1.
func (x *Extension[K]) NewOne() *Extension[K] {
	y := &Extension[K]{}
	y.Irr = poly0(x.Irr).Set(x.Irr)
	y.Poly = poly1(y.Irr)
	return y
}

// Equal reports whether x and y are equal.
func (x *Extension[K]) Equal(y *Extension[K]) bool {
	if !x.Irr.Equal(y.Irr) {
		return false
	}
	if !x.Poly.Equal(y.Poly) {
		return false
	}
	return true
}

// Set sets x to y.
func (x *Extension[K]) Set(y *Extension[K]) *Extension[K] {
	x.Irr.Set(y.Irr)
	x.Poly.Set(y.Poly)
	return x
}

// Add sets z to the sum x+y and returns z.
func (z *Extension[K]) Add(x, y *Extension[K]) *Extension[K] {
	z.Irr.Set(x.Irr)
	z.Poly.Add(x.Poly, y.Poly)
	_, remainder := divide(z.Poly, z.Irr)
	z.Poly.Set(remainder)
	return z
}

// Sub sets z to the difference x-y and returns z.
func (z *Extension[K]) Sub(x, y *Extension[K]) *Extension[K] {
	z.Irr.Set(x.Irr)

	// Compute negative of y.Poly.
	k := z.Irr.Field()
	neg1 := k.NewZero()
	neg1.Sub(neg1, k.NewOne())
	negY := mulScalar(y.Poly, neg1)

	z.Poly.Add(x.Poly, negY)
	_, remainder := divide(z.Poly, z.Irr)
	z.Poly.Set(remainder)
	return z
}

// Mul sets z to the product x*y and returns z.
func (z *Extension[K]) Mul(x, y *Extension[K]) *Extension[K] {
	z.Irr.Set(x.Irr)
	z.Poly.Mul(x.Poly, y.Poly)
	_, remainder := divide(z.Poly, z.Irr)
	z.Poly.Set(remainder)
	return z
}

// Div sets z to the quotient x/y and returns z.
func (z *Extension[K]) Div(x, y *Extension[K]) *Extension[K] {
	yInv := z.NewZero().Inv(y)
	z.Mul(x, yInv)
	return z
}

// Inv sets z to 1/x and returns z.
func (z *Extension[K]) Inv(x *Extension[K]) *Extension[K] {
	if x.Poly.Len() == 0 {
		panic("no inverse for 0")
	}
	z.Irr.Set(x.Irr)
	z.Poly.Set(inverse(x.Poly, z.Irr))
	return z
}

// String returns the coefficient representation of x.
func (x *Extension[K]) String() string {
	coeffs := x.Coeffs()
	zero := coeffs[0].NewZero()
	isScalar := true
	for i := 1; i < len(coeffs); i++ {
		if !coeffs[i].Equal(zero) {
			isScalar = false
			break
		}
	}
	if isScalar {
		return coeffs[0].String()
	}
	return fmt.Sprintf("%v", x.Coeffs())
}

// Degree returns the [degree] of the field extension.
//
// [degree]: https://en.wikipedia.org/wiki/Degree_of_a_field_extension
func (x *Extension[K]) Degree() int {
	return len(x.Irr.LeadingTerm().Monomial)
}

// SetCoeffs sets the coefficients of the polynomial representation of x.
func (x *Extension[K]) SetCoeffs(coefficients ...K) *Extension[K] {
	k := x.Irr.Field()
	x.Poly = poly0(x.Poly)
	for deg, coeff := range coefficients {
		c := k.NewZero().Set(coeff)
		w := make([]nag.Symbol, deg)
		x.Poly.Add(x.Poly, nag.NewPolynomial(k, x.Poly.Order(), nag.PolynomialTerm[K]{Coefficient: c, Monomial: w}))
	}
	return x
}

// Coeffs returns the coefficients of x's polynomial representation.
func (x *Extension[K]) Coeffs() []K {
	k := x.Irr.Field()
	coeffs := make([]K, x.Degree())
	for i := range coeffs {
		coeffs[i] = k.NewZero()
	}

	for c, w := range x.Poly.Terms() {
		d := len(w)
		coeffs[d].Set(c)
	}
	return coeffs
}

// A Finite is an element of a finite field GF(q) with order q = p^n
// where p is a prime number, and n is a positive integer.
type Finite[K nag.Field[K]] interface {
	nag.Field[K]
	// SetCoeffs sets the polynomial representation of the receiver.
	SetCoeffs(...*big.Int) K
	// Coeffs returns the polynomial representation.
	Coeffs() []*big.Int
	// Degree returns the degree of the finite field.
	Degree() int
	// Characteristic returns the characteristic the finite field.
	Characteristic() *big.Int
}

// A PrimeExt is an element of a finite field GF(q) with order q = p^n
// where p is a prime number, and n is a positive integer.
type PrimeExt Extension[*prime]

// NewPrimeExt returns an element of the finite field GF(p, len(irr)-1).
// The coefficients of the irreducible polynomial defining the field is irr.
func NewPrimeExt(p *big.Int, irr []*big.Int) *PrimeExt {
	k := &prime{order: new(big.Int).Set(p), i: big.NewInt(0)}
	poly := nag.NewPolynomial(k, nag.Deglex)
	poly.SymbolStringer = func(nag.Symbol) string { return "x" }
	for deg, cIth := range irr {
		c := setIth(k.NewZero(), cIth)
		w := make([]nag.Symbol, deg)
		poly.Add(poly, nag.NewPolynomial(k, poly.Order(), nag.PolynomialTerm[*prime]{Coefficient: c, Monomial: w}))
	}
	x := &PrimeExt{Irr: poly}
	x.Poly = poly0(x.Irr)
	return x
}

// NewPrimeExtDeg returns an element of the finite field GF(p, n).
// The irreducible polynomial is automatically calculated.
func NewPrimeExtDeg(p *big.Int, n int) *PrimeExt {
	k := &prime{order: new(big.Int).Set(p), i: big.NewInt(0)}
	xn := nag.NewPolynomial(k, nag.Deglex, nag.PolynomialTerm[*prime]{Coefficient: k.NewOne(), Monomial: make([]nag.Symbol, n)})
	xn.SymbolStringer = func(nag.Symbol) string { return "x" }
	irr := nag.NewPolynomial(xn.Field(), xn.Order())

	order := new(big.Int).Exp(p, big.NewInt(int64(n)), nil)
	for i := big.NewInt(1); i.Cmp(order) < 0; i.Add(i, big.NewInt(1)) {
		irr.Set(xn)
		ithPoly(irr, i, p)

		factors := Factor(irr)
		if !(len(factors) == 1 && factors[0].N.Cmp(big.NewInt(1)) == 0) {
			continue
		}
		// Only check primitivity if for small orders, since this involves iterating through every element in the field.
		if order.Cmp(big.NewInt(1024)) < 0 && !primitive(irr) {
			continue
		}
		return &PrimeExt{Irr: irr, Poly: poly0(irr)}
	}
	panic("unable to find irreducible polynomial")
}

// NewZero returns the additive identity 0.
func (x *PrimeExt) NewZero() *PrimeExt {
	return (*PrimeExt)((*Extension[*prime])(x).NewZero())
}

// NewOne returns the multiplicative identity 1.
func (x *PrimeExt) NewOne() *PrimeExt {
	return (*PrimeExt)((*Extension[*prime])(x).NewOne())
}

// Equal reports whether x and y are equal.
func (x *PrimeExt) Equal(y *PrimeExt) bool {
	return (*Extension[*prime])(x).Equal((*Extension[*prime])(y))
}

// Set sets x to y.
func (x *PrimeExt) Set(y *PrimeExt) *PrimeExt {
	return (*PrimeExt)((*Extension[*prime])(x).Set((*Extension[*prime])(y)))
}

// Add sets z to the sum x+y and returns z.
func (z *PrimeExt) Add(x, y *PrimeExt) *PrimeExt {
	return (*PrimeExt)((*Extension[*prime])(z).Add((*Extension[*prime])(x), (*Extension[*prime])(y)))
}

// Sub sets z to the difference x-y and returns z.
func (z *PrimeExt) Sub(x, y *PrimeExt) *PrimeExt {
	return (*PrimeExt)((*Extension[*prime])(z).Sub((*Extension[*prime])(x), (*Extension[*prime])(y)))
}

// Mul sets z to the product x*y and returns z.
func (z *PrimeExt) Mul(x, y *PrimeExt) *PrimeExt {
	return (*PrimeExt)((*Extension[*prime])(z).Mul((*Extension[*prime])(x), (*Extension[*prime])(y)))
}

// Div sets z to the quotient x/y and returns z.
func (z *PrimeExt) Div(x, y *PrimeExt) *PrimeExt {
	return (*PrimeExt)((*Extension[*prime])(z).Div((*Extension[*prime])(x), (*Extension[*prime])(y)))
}

// Inv sets x to 1/y and returns x.
func (x *PrimeExt) Inv(y *PrimeExt) *PrimeExt {
	return (*PrimeExt)((*Extension[*prime])(x).Inv((*Extension[*prime])(y)))
}

// String returns the coefficient representation of x.
func (x *PrimeExt) String() string {
	return (*Extension[*prime])(x).String()
}

// SetCoeffs sets the coefficients of the polynomial representation of x.
func (x *PrimeExt) SetCoeffs(cs ...*big.Int) *PrimeExt {
	k := x.Poly.Field()
	coeffs := make([]*prime, len(cs))
	for i := range coeffs {
		coeffs[i] = k.NewZero()
		coeffs[i].i.Set(cs[i])
	}
	return (*PrimeExt)((*Extension[*prime])(x).SetCoeffs(coeffs...))
}

// Coeffs returns the coefficients of x's polynomial representation.
func (x *PrimeExt) Coeffs() []*big.Int {
	coeffs := (*Extension[*prime])(x).Coeffs()
	cs := make([]*big.Int, len(coeffs))
	for i := range cs {
		cs[i] = new(big.Int).Set(coeffs[i].i)
	}
	return cs
}

// Degree returns the [degree] of the field extension.
//
// [degree]: https://en.wikipedia.org/wiki/Degree_of_a_field_extension
func (x *PrimeExt) Degree() int {
	return (*Extension[*prime])(x).Degree()
}

// Characteristic returns the [characteristic] of the field.
//
// [characteristic]: https://en.wikipedia.org/wiki/Characteristic_(algebra)
func (x *PrimeExt) Characteristic() *big.Int {
	return x.Irr.LeadingTerm().Coefficient.Characteristic()
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

func (x *prime) Set(y *prime) *prime {
	x.order.Set(y.order)
	x.i.Set(y.i)
	return x
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

func (x *prime) Degree() int { return 1 }

func (x *prime) SetCoeffs(coefficients ...*big.Int) *prime {
	x.i.Set(coefficients[0])
	return x
}

func (x *prime) Coeffs() []*big.Int {
	return []*big.Int{x.i}
}

// setIth sets x to the i'th element in the field and returns x.
func setIth[K Finite[K]](x K, i *big.Int) K {
	i = new(big.Int).Set(i)
	p := x.Characteristic()
	r := new(big.Int)
	coeffs := make([]*big.Int, 0)
	for i.Sign() != 0 {
		i.QuoRem(i, p, r)
		coeffs = append(coeffs, new(big.Int).Set(r))
	}
	if len(coeffs) == 0 {
		coeffs = append(coeffs, big.NewInt(0))
	}
	return x.SetCoeffs(coeffs...)
}
