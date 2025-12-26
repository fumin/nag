// Package nag implements algorithms in noncommutative algebraic geometry.
// In particular, this package provides functions to perform [polynomial division] and compute [Gröbner bases].
// Applications of this package include simplifying expressions and solving equations containing noncommutative variables.
//
// [polynomial division]: https://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor#Euclidean_division
// [Gröbner bases]: https://en.wikipedia.org/wiki/Gr%C3%B6bner_basis
package nag

import (
	"bytes"
	"cmp"
	"fmt"
	"iter"
	"math"
	"math/big"
	"reflect"
	"slices"
	"strings"

	"github.com/jba/omap"
)

// A Field is an element whose addition and multiplication operations satisfy the [field] axioms.
//
// [field]: https://en.wikipedia.org/wiki/Field_(mathematics)
type Field[T any] interface {
	// NewZero returns the additive identity of the field.
	NewZero() T
	// NewOne returns the multiplicative identity of the field.
	NewOne() T

	// Equal reports whether x and y are equal, where x is the method receiver.
	Equal(y T) bool
	// Add sets z to the sum x+y and returns z, where z is the method receiver.
	Add(x, y T) T
	// Sub sets z to the difference x-y and returns z, where z is the method receiver.
	Sub(x, y T) T
	// Mul sets z to the product x*y and returns z, where z is the method receiver.
	Mul(x, y T) T
	// Div sets z to the quotient x/y and returns z, where z is the method receiver.
	Div(x, y T) T
	// Inv sets z to 1/x and returns z, where z is the method receiver.
	Inv(x T) T

	// String returns the string representation.
	String() string
}

// A Rat represents a quotient of arbitrary precision.
type Rat struct{ *big.Rat }

// NewRat creates a new [Rat] with numerator a and denominator b.
func NewRat(a, b int64) *Rat { return &Rat{big.NewRat(a, b)} }

// NewZero returns the additive identity 0.
func (x *Rat) NewZero() *Rat {
	return &Rat{big.NewRat(0, 1)}
}

// NewOne returns the multiplicative identity 1.
func (x *Rat) NewOne() *Rat {
	return &Rat{big.NewRat(1, 1)}
}

// Add sets z to the sum x+y and returns z.
func (z *Rat) Add(x, y *Rat) *Rat { return &Rat{z.Rat.Add(x.Rat, y.Rat)} }

// Sub sets z to the difference x-y and returns z.
func (z *Rat) Sub(x, y *Rat) *Rat { return &Rat{z.Rat.Sub(x.Rat, y.Rat)} }

// Mul sets z to the product x*y and returns z.
func (z *Rat) Mul(x, y *Rat) *Rat { return &Rat{z.Rat.Mul(x.Rat, y.Rat)} }

// Div sets z to the quotient x/y and returns z. If y == 0, Div panics.
func (z *Rat) Div(x, y *Rat) *Rat { return &Rat{z.Rat.Quo(x.Rat, y.Rat)} }

// Inv sets z to 1/x and returns z. If x == 0, Inv panics.
func (z *Rat) Inv(x *Rat) *Rat { return &Rat{z.Rat.Inv(x.Rat)} }

// Equal reports whether x and y are equal.
func (x *Rat) Equal(y *Rat) bool {
	return x.Rat.Cmp(y.Rat) == 0
}

// String returns a string representation of x in the form "a/b" if b != 1, and in the form "a" if b == 1.
func (x *Rat) String() string {
	return x.RatString()
}

// A Symbol is a variable in a monomial.
// This definition implies that the maximum number of variables supported by this package is 256.
type Symbol = byte

// A [Monomial] is a product of variables chained by multiplication.
//
// [Monomial]: https://en.wikipedia.org/wiki/Monomial
type Monomial []Symbol

// An Order is a [monomial order] for comparing monomials.
// The meaning of the return value is the same as [cmp.Compare].
//
// [monomial order]: https://en.wikipedia.org/wiki/Monomial_order
type Order func(x, y Monomial) int

// [Deglex] compares x, y by first comparing their degrees, and in case of a tie applies the lexicographic order.
//
// [Deglex]: https://en.wikipedia.org/wiki/Monomial_order#Graded_lexicographic_order
func Deglex(x, y Monomial) int {
	if c := cmp.Compare(len(x), len(y)); c != 0 {
		return c
	}
	return lexicographic(x, y)
}

// ElimOrder returns a monomial order that first compares monomials as commutative words lexicographically, and in case of a tie applies noncommutative lexicographic order.
func ElimOrder() Order {
	var xb, yb Monomial
	order := func(x, y Monomial) int {
		xb, yb = xb[:0], yb[:0]
		xb = append(xb, x...)
		yb = append(yb, y...)

		// Compare as commutative monomials.
		slices.SortFunc(xb, func(a, b Symbol) int { return -cmp.Compare(a, b) })
		slices.SortFunc(yb, func(a, b Symbol) int { return -cmp.Compare(a, b) })
		if c := lexicographic(xb, yb); c != 0 {
			return c
		}

		return lexicographic(x, y)
	}
	return order
}

// A PolynomialTerm is a term in a polynomial.
type PolynomialTerm[K Field[K]] struct {
	Coefficient K
	Monomial    Monomial
}

// A Polynomial is a polynomial of noncommutative variables.
type Polynomial[K Field[K]] struct {
	// SymbolStringer specifies how a symbol in a monomial is formated when the polynomial is printed out.
	SymbolStringer func(s Symbol) string

	field K
	order Order
	m     *omap.MapFunc[Monomial, K]
}

// NewPolynomial returns a new polynomial containing the given terms.
func NewPolynomial[K Field[K]](field K, order Order, terms ...PolynomialTerm[K]) *Polynomial[K] {
	x := &Polynomial[K]{
		SymbolStringer: englishSymbolStringer,
		field:          field,
		order:          order,
		m:              omap.NewMapFunc[Monomial, K](order),
	}
	for _, term := range terms {
		x.addTerm(1, term)
	}
	return x
}

// Field returns the field of the coefficients in x.
func (x *Polynomial[K]) Field() K { return x.field }

// Order returns the monomial order employed by x.
func (x *Polynomial[K]) Order() Order { return x.order }

// Len reports the number of terms in x.
func (x *Polynomial[K]) Len() int { return x.m.Len() }

// Terms iterates the terms in a polynomial.
func (x *Polynomial[K]) Terms() iter.Seq2[K, Monomial] {
	return func(yield func(K, Monomial) bool) {
		for w, c := range x.m.Backward() {
			if !yield(c, w) {
				return
			}
		}
	}
}

// Equal reports whether x and y have the same coefficients and monomials.
func (x *Polynomial[K]) Equal(y *Polynomial[K]) bool {
	if x.m.Len() != y.m.Len() {
		return false
	}
	for i := range x.m.Len() {
		xw, xc := x.m.At(x.m.Len() - 1 - i)
		yw, yc := y.m.At(y.m.Len() - 1 - i)
		if !monomialEq(xw, yw) {
			return false
		}
		if !xc.Equal(yc) {
			return false
		}
	}
	return true
}

// Set sets z to x and returns z.
func (z *Polynomial[K]) Set(x *Polynomial[K]) *Polynomial[K] {
	if z == x {
		return z
	}
	z.SymbolStringer = x.SymbolStringer
	z.field = x.field
	z.order = x.order
	z.m = omap.NewMapFunc[Monomial, K](z.order)
	for xw, xc := range x.m.All() {
		w := make(Monomial, len(xw))
		copy(w, xw)
		z.addTerm(1, PolynomialTerm[K]{Coefficient: xc, Monomial: w})
	}
	return z
}

// Add sets z to the sum x+y and returns z.
func (z *Polynomial[K]) Add(x, y *Polynomial[K]) *Polynomial[K] {
	// Set z = x, while handling the case where x or y is z itself.
	if y == z {
		x, y = y, x
	}
	if z != x {
		z.m.Clear()
		for xw, c := range x.m.All() {
			w := make(Monomial, len(xw))
			copy(w, xw)
			z.addTerm(1, PolynomialTerm[K]{Coefficient: c, Monomial: w})
		}
	}

	// Compute z += y.
	for yw, c := range y.m.All() {
		w := make(Monomial, len(yw))
		copy(w, yw)
		z.addTerm(1, PolynomialTerm[K]{Coefficient: c, Monomial: w})
	}

	return z
}

// Mul sets z to the product x*y and returns z.
func (z *Polynomial[K]) Mul(x, y *Polynomial[K]) *Polynomial[K] {
	if z == x {
		panic(fmt.Sprintf("z == x"))
	}
	if z == y {
		panic(fmt.Sprintf("z == y"))
	}

	z.m.Clear()
	for xw, xc := range x.m.Backward() {
		for yw, yc := range y.m.Backward() {
			c := z.field.Mul(xc, yc)
			w := make(Monomial, 0, len(xw)+len(yw))
			w = append(append(w, xw...), yw...)
			z.addTerm(1, PolynomialTerm[K]{Coefficient: c, Monomial: w})
		}
	}

	return z
}

// Pow sets z to the power x^y and returns z.
func (z *Polynomial[K]) Pow(x *Polynomial[K], y int) *Polynomial[K] {
	if z == x {
		panic("z == x")
	}

	z.Set(x)
	buf := NewPolynomial[K](z.field, z.order)
	for range y - 1 {
		buf.Mul(z, x)
		z, buf = buf, z
	}
	if y%2 == 0 {
		z, buf = buf, z
		z.Set(buf)
	}

	return z
}

// LeadingTerm returns the polynomial term of the leading monomial.
// Note that the leading term depends on the monomial order employed by the polynomial.
func (x *Polynomial[K]) LeadingTerm() PolynomialTerm[K] {
	w, ok := x.m.Max()
	if !ok {
		panic("zero polynomial has no terms")
	}
	c, _ := x.m.Get(w)
	return PolynomialTerm[K]{Coefficient: c, Monomial: w}
}

// String returns the string representation of x.
// Symbols in x are formatted using x.SymbolStringer.
func (x *Polynomial[K]) String() string {
	if x.Len() == 0 {
		return "0"
	}
	var b strings.Builder
	for i := range x.m.Len() {
		w, c := x.m.At(x.m.Len() - 1 - i)

		// Print c.
		s := c.String()
		if s[0] != '-' {
			s = "+" + s
		}
		switch {
		case i == 0 && s == "+1" && len(w) != 0:
			s = ""
		case i == 0 && s[0] == '+':
			s = s[1:]
		case s == "+1" && len(w) != 0:
			s = "+"
		case s == "-1" && len(w) != 0:
			s = "-"
		}
		fmt.Fprintf(&b, "%s", s)

		// Print w.
		printMonomial(&b, w, x.SymbolStringer)
	}
	return b.String()
}

func (x *Polynomial[K]) addTerm(sign int, term PolynomialTerm[K]) {
	c, ok := x.m.Get(term.Monomial)
	if !ok {
		c = x.field.NewZero()
	}

	tc := term.Coefficient
	tcv := reflect.ValueOf(tc)
	kind := tcv.Kind()
	if (kind == reflect.Ptr || kind == reflect.Interface) && tcv.IsNil() {
		tc = x.field.NewOne()
	}
	if sign < 0 {
		c.Sub(c, tc)
	} else {
		c.Add(c, tc)
	}

	if c.Equal(x.field.NewZero()) {
		x.m.Delete(term.Monomial)
	} else {
		x.m.Set(term.Monomial, c)
	}
}

func (z *Polynomial[K]) add(sign int, c K, left Monomial, x *Polynomial[K], right Monomial) {
	for xw, xc := range x.m.Backward() {
		c := z.field.Mul(c, xc)
		w := make(Monomial, 0, len(left)+len(xw)+len(right))
		w = append(append(append(w, left...), xw...), right...)
		z.addTerm(sign, PolynomialTerm[K]{Coefficient: c, Monomial: w})
	}
}

func (z *Polynomial[K]) mulScalar(scalar K, x *Polynomial[K]) *Polynomial[K] {
	if z == x {
		for zw, zc := range z.m.All() {
			zc.Mul(scalar, zc)
			z.m.Set(zw, zc)
		}
		return z
	}

	z.m.Clear()
	for xw, xc := range x.m.All() {
		c := z.field.Mul(scalar, xc)
		w := make(Monomial, len(xw))
		copy(w, xw)
		z.addTerm(1, PolynomialTerm[K]{Coefficient: c, Monomial: w})
	}
	return z
}

// A Quotient is the resulting quotient of a polynomial division.
// Let f be a polynomial and g an ideal, the result of f divided by g is:
//
//	f = g*quotient + remainder
//
// For the exact form of a quotient, see Theorem 3.2.1, Xiu Xingqiang.
//
// Xiu, Xingqiang. "Non-commutative Gröbner bases and applications." PhD diss., Universität Passau, 2012.
type Quotient[K Field[K]] struct {
	// Coefficient is c_{ij} in Theorem 3.2.1, Xiu Xingqiang.
	Coefficient K
	// Coefficient is w_{ij} in Theorem 3.2.1, Xiu Xingqiang.
	Left Monomial
	// Coefficient is w'_{ij} in Theorem 3.2.1, Xiu Xingqiang.
	Right Monomial
}

// Divide divides the polynomial f by the ideal g, and returns the quotient and remainder.
// The polynomial f is modified upon return.
// For more details, please see Theorem 3.2.1, Xiu Xingqiang.
//
// Xiu, Xingqiang. "Non-commutative Gröbner bases and applications." PhD diss., Universität Passau, 2012.
func Divide[K Field[K]](quotient [][]Quotient[K], f *Polynomial[K], g []*Polynomial[K]) (outQuotient [][]Quotient[K], remainder *Polynomial[K]) {
	if quotient != nil {
		short := len(g) - len(quotient)
		if short > 0 {
			quotient = append(quotient, make([][]Quotient[K], short)...)
		}
		quotient = quotient[:len(g)]
		for i := range quotient {
			quotient[i] = quotient[i][:0]
		}
	}
	p := NewPolynomial[K](f.field, f.order)
	p.SymbolStringer = f.SymbolStringer
	v := f

	for v.m.Len() != 0 {
		lmv := v.LeadingTerm()
		ltv := lmv.Monomial

		// Find basis where ltv = left * ltg * right.
		basis, leftEnd, lmg := -1, -1, PolynomialTerm[K]{}
		for i, gi := range g {
			if gi == nil {
				continue
			}
			lmg = gi.LeadingTerm()
			ltg := lmg.Monomial

			if leftEnd = monomialIndex(ltv, ltg); leftEnd != -1 {
				basis = i
				break
			}
		}

		if basis == -1 {
			p.addTerm(1, lmv)
			v.addTerm(-1, lmv)
		} else {
			q := Quotient[K]{
				Coefficient: f.field.NewZero().Div(lmv.Coefficient, lmg.Coefficient),
				Left:        ltv[:leftEnd],
				Right:       ltv[leftEnd+len(lmg.Monomial):],
			}
			if quotient != nil {
				quotient[basis] = append(quotient[basis], q)
			}
			v.add(-1, q.Coefficient, q.Left, g[basis], q.Right)
		}
	}

	return quotient, p
}

// Buchberger returns the Gröbner basis of the ideal g, using the Buchberger algorithm.
// For noncummutative algebras, a Gröbner basis may not be finite and therefore only a partial basis can be computed.
// In this case, upon reaching maxIter without a complete basis, Buchberger sets complete to false.
// For more details, please see Theorem 4.2.24, Xiu Xingqiang.
//
// Xiu, Xingqiang. "Non-commutative Gröbner bases and applications." PhD diss., Universität Passau, 2012.
func Buchberger[K Field[K]](g []*Polynomial[K], maxIter int) (basis []*Polynomial[K], complete bool) {
	// Buffers.
	r0 := g[0].field.NewZero()
	buf := &Monomial{}

	g = interreduce(g)
	// t tracks unwanted polynomials in g.
	t := make([]*Polynomial[K], len(g))
	deleteUnwanted := func() {
		for i := range t {
			if t[i] != nil {
				g[i] = nil
			}
		}
	}
	restoreUnwanted := func() {
		for i := range t {
			if t[i] != nil {
				g[i] = t[i]
			}
		}
	}
	// b is the set of obstructions.
	var b []obstruction[K]
	for l := 1; l <= len(g); l++ {
		gl := g[:l]
		b = addObstructions(b, gl, buf)
	}

	for range maxIter {
		if len(b) == 0 {
			complete = true
			break
		}

		// Take oij from b.
		oij := b[0]
		for i := range len(b) - 1 {
			b[i] = b[i+1]
		}
		b = b[:len(b)-1]

		// Compute S-polynomial of oij.
		s := sPolynomial[K](oij, g, r0)
		deleteUnwanted()
		_, sP := Divide(nil, s, g)
		restoreUnwanted()
		if sP.m.Len() == 0 {
			continue
		}

		// Add sP to g and add new obstructions.
		g = append(g, sP)
		t = append(t, nil)
		b = addObstructions(b, g, buf)

		// Set gi as unwanted if ltgi is a multiple of ltgs.
		ltgs := g[len(g)-1].LeadingTerm().Monomial
		for i := range len(g) - 1 {
			ltgi := g[i].LeadingTerm().Monomial
			if monomialIndex(ltgi, ltgs) != -1 {
				t[i] = g[i]
			}
		}
	}

	// Remove unwanted polynomials.
	deleteUnwanted()
	g = slices.DeleteFunc(g, func(gi *Polynomial[K]) bool { return gi == nil })
	g = interreduce(g)
	// Make basis monic.
	for i := range g {
		lc := g[i].LeadingTerm().Coefficient
		g[i].mulScalar(r0.Inv(lc), g[i])
	}
	slices.SortFunc(g, polynomialCmp[K])
	return g, complete
}

// BuchbergerHomogeneous returns the Gröbner basis of the [homogeneous] ideal g, using the Buchberger algorithm.
// All basis polynomials with degree less or equal than maxDeg are returned.
// For more details, please see Theorem 4.3.16, Xiu Xingqiang.
//
// Xiu, Xingqiang. "Non-commutative Gröbner bases and applications." PhD diss., Universität Passau, 2012.
//
// [homogeneous]: https://en.wikipedia.org/wiki/Homogeneous_polynomial
func BuchbergerHomogeneous[K Field[K]](g []*Polynomial[K], maxDeg int) (basis []*Polynomial[K], complete bool) {
	// Check that g is a homogeneous ideal.
	for _, f := range g {
		if !homogeneous(f) {
			panic(fmt.Sprintf("non-homogeneous polynomial %v", f))
		}
	}

	// Buffers.
	r0 := g[0].field.NewZero()
	m0 := &Monomial{}
	p0 := NewPolynomial(g[0].field, g[0].order)

	// Make a copy of g since we will be deleting from it.
	newG := make([]*Polynomial[K], len(g))
	copy(newG, g)
	g = newG
	// b is the set of obstructions.
	var b []obstruction[K]
	// gd and bd are subsets of g and b with degree d.
	var gd []*Polynomial[K]
	var bd []obstruction[K]
	var numDeleted int

	for {
		if len(g) == 0 && len(b) == 0 {
			if numDeleted == 0 {
				complete = true
			}
			break
		}
		var stMax bool
		g, gd, b, bd, stMax = smallestDegreeSet(g, gd, b, bd, maxDeg)
		if !stMax {
			break
		}

		for len(gd) > 0 {
			gL := gd[len(gd)-1]
			gd = gd[:len(gd)-1]
			_, gP := Divide(nil, p0.Set(gL), basis)
			if gP.m.Len() == 0 {
				continue
			}

			basis = append(basis, gP)
			b = addObstructions(b, basis, m0)
			var nDel int
			b, nDel = deleteHighDegObs(b, basis, maxDeg, r0)
			numDeleted += nDel
		}

		for len(bd) > 0 {
			o := bd[len(bd)-1]
			bd = bd[:len(bd)-1]
			_, sP := Divide(nil, p0.Set(o.sPolynomial), basis)
			if sP.m.Len() == 0 {
				continue
			}

			basis = append(basis, sP)
			b = addObstructions(b, basis, m0)
			var nDel int
			b, nDel = deleteHighDegObs(b, basis, maxDeg, r0)
			numDeleted += nDel
		}
	}

	basis = interreduce(basis)
	// Make basis monic.
	for i := range basis {
		lc := basis[i].LeadingTerm().Coefficient
		basis[i].mulScalar(r0.Inv(lc), basis[i])
	}
	slices.SortFunc(basis, polynomialCmp[K])
	return basis, complete
}

func interreduce[K Field[K]](g []*Polynomial[K]) []*Polynomial[K] {
	i, s := 0, len(g)
	for i != s {
		gi := g[i]
		if gi == nil {
			i++
			continue
		}
		f := NewPolynomial(g[0].field, g[0].order).Set(gi)
		g[i] = nil
		_, giP := Divide(nil, f, g)

		switch {
		case giP.m.Len() == 0:
			g[i] = nil
			i++
		case !giP.Equal(gi):
			g[i] = giP
			i = 0
		default:
			g[i] = gi
			i++
		}
	}
	return slices.DeleteFunc(g, func(x *Polynomial[K]) bool { return x == nil })
}

func smallestDegreeSet[K Field[K]](g, gd []*Polynomial[K], b, bd []obstruction[K], maxDeg int) ([]*Polynomial[K], []*Polynomial[K], []obstruction[K], []obstruction[K], bool) {
	// Compute the smallest degree.
	var d int = math.MaxInt
	for _, f := range g {
		d = min(d, len(f.LeadingTerm().Monomial))
	}
	for _, o := range b {
		d = min(d, len(o.sPolynomial.LeadingTerm().Monomial))
	}
	if d > maxDeg {
		return nil, nil, nil, nil, false
	}

	// Extract from g.
	gd = gd[:0]
	for i, f := range g {
		if len(f.LeadingTerm().Monomial) == d {
			gd = append(gd, f)
			g[i] = nil
		}
	}
	g = slices.DeleteFunc(g, func(f *Polynomial[K]) bool { return f == nil })

	//Extract from b.
	bd = bd[:0]
	for i, o := range b {
		if len(o.sPolynomial.LeadingTerm().Monomial) == d {
			bd = append(bd, o)
			b[i].sPolynomial = nil
		}
	}
	b = slices.DeleteFunc(b, func(o obstruction[K]) bool { return o.sPolynomial == nil })

	return g, gd, b, bd, true
}

func deleteHighDegObs[K Field[K]](b []obstruction[K], basis []*Polynomial[K], maxDeg int, r0 K) ([]obstruction[K], int) {
	var numDeleted int
	for i := len(b) - 1; i >= 0; i-- {
		if b[i].sPolynomial != nil {
			break
		}
		b[i].sPolynomial = sPolynomial(b[i], basis, r0)
		if b[i].sPolynomial.Len() == 0 {
			b[i].removed = true
			continue
		}
		if len(b[i].sPolynomial.LeadingTerm().Monomial) > maxDeg {
			b[i].removed = true
			numDeleted++
		}
	}
	b = slices.DeleteFunc(b, func(o obstruction[K]) bool { return o.removed })
	return b, numDeleted
}

type obstruction[K Field[K]] struct {
	i           int
	j           int
	iLeft       Monomial
	iRight      Monomial
	jLeft       Monomial
	jRight      Monomial
	sPolynomial *Polynomial[K]

	removed bool
}

func sPolynomial[K Field[K]](o obstruction[K], g []*Polynomial[K], buf K) *Polynomial[K] {
	gi, gj := g[o.i], g[o.j]
	lcgi := gi.LeadingTerm().Coefficient
	lcgj := gj.LeadingTerm().Coefficient
	s := NewPolynomial(gi.field, gi.order)
	s.SymbolStringer = gi.SymbolStringer
	s.add(1, buf.Inv(lcgi), o.iLeft, gi, o.iRight)
	s.add(-1, buf.Inv(lcgj), o.jLeft, gj, o.jRight)
	return s
}

func delRemoved[K Field[K]](sPObs, obs []obstruction[K]) ([]obstruction[K], []obstruction[K]) {
	spPrevLen := len(sPObs)
	sPObs = slices.DeleteFunc(sPObs, func(o obstruction[K]) bool { return o.removed })
	obs = obs[:len(obs)-(spPrevLen-len(sPObs))]
	return sPObs, obs
}

func addObstructions[K Field[K]](obs []obstruction[K], g []*Polynomial[K], buf *Monomial) []obstruction[K] {
	// Add sPObs.
	prevLen := len(obs)
	obs = overlapObstruction(obs, g)
	sPObs := obs[prevLen:]

	// Remove from sPObs using step 4b Theorem 4.2.22.
	remove4b(sPObs, g[0].order)
	sPObs, obs = delRemoved(sPObs, obs)

	// Remove from sPObs using step 4c Theorem 4.2.22.
	b := obs[:prevLen]
	ltgs := g[len(g)-1].LeadingTerm().Monomial
	remove4c(sPObs, b, ltgs)
	sPObs, obs = delRemoved(sPObs, obs)

	// Remove from b using step 4d Theorem 4.2.22.
	remove4d(b, sPObs, g, buf)
	obs = slices.DeleteFunc(obs, func(o obstruction[K]) bool { return o.removed })

	return obs
}

func remove4b[K Field[K]](sPObs []obstruction[K], order Order) {
	for k := range sPObs {
		oi := sPObs[k]
		for l := range sPObs {
			if l == k {
				continue
			}
			oj := sPObs[l]
			if oj.removed {
				continue
			}

			// Check that oi encloses oj.
			w, ok := cutSuffix(oi.jLeft, oj.jLeft)
			if !ok {
				continue
			}
			wP, ok := cutPrefix(oi.jRight, oj.jRight)
			if !ok {
				continue
			}
			wwp1 := ((len(w) == 0) && (len(wP) == 0))
			wiLtwj := (order(oi.iLeft, oj.iLeft) == 1)

			// Remove oi if possible.
			var removed bool
			i, j := oi.i, oj.i
			switch {
			case i > j:
				removed = true
			case i <= j && !wwp1:
				removed = true
			case i == j && wwp1 && wiLtwj:
				removed = true
			}
			if removed {
				sPObs[k].removed = removed
				break
			}
		}
	}
}

func remove4c[K Field[K]](sPObs, b []obstruction[K], ltgs Monomial) {
	for k := range sPObs {
		oj := sPObs[k]
		for l := range b {
			oi := b[l]
			if oi.j != oj.i {
				continue
			}

			w, ok := cutSuffix(oj.iLeft, oi.jLeft)
			if !ok {
				continue
			}
			wP, ok := cutPrefix(oj.iRight, oi.jRight)
			if !ok {
				continue
			}

			// Check whether w * oi.iLeft is a multiple of oj.jLeft * ltgs
			multiple := len(oj.jLeft)+len(ltgs) <= len(w)+len(oi.iLeft)
			// Check whether oi.iRight * wP is a multiple of ltgs * oj.jRight
			multipleP := len(ltgs)+len(oj.jRight) <= len(oi.iRight)+len(wP)
			if multiple || multipleP {
				sPObs[k].removed = true
				break
			}
		}
	}
}

func remove4d[K Field[K]](b, sPObs []obstruction[K], g []*Polynomial[K], buf *Monomial) {
	sP := len(g) - 1
	ltgs := g[sP].LeadingTerm().Monomial
	for k := range b {
		oij := b[k]

		ltgi := g[oij.i].LeadingTerm().Monomial
		ltgj := g[oij.j].LeadingTerm().Monomial
		*buf = (*buf)[:0]
		*buf = append(append(append(*buf, oij.jLeft...), ltgj...), oij.jRight...)
		wjw := *buf
		var wjwStart int
		for wjwStart < len(wjw) {
			// Find ws and wsP such that, ws * ltgs * wsP == wjw.
			sEnd := monomialIndex(wjw[wjwStart:], ltgs)
			if sEnd == -1 {
				break
			}
			ws, wsP := wjw[:wjwStart+sEnd], wjw[wjwStart+sEnd+len(ltgs):]
			wjwStart += (sEnd + 1)

			// Check o(i, s).
			ois := obstruction[K]{i: oij.i, j: sP, iLeft: oij.iLeft, iRight: oij.iRight, jLeft: ws, jRight: wsP}
			iNoOverlap := !hasOverlap(ois, ltgi, ltgs)
			ois = shrink(ois)
			iInS := slices.ContainsFunc(sPObs, func(o obstruction[K]) bool { return obsEq(o, ois) })
			if !(iNoOverlap || iInS) {
				continue
			}

			// Check o(j, s).
			ojs := obstruction[K]{i: oij.j, j: sP, iLeft: oij.jLeft, iRight: oij.jRight, jLeft: ws, jRight: wsP}
			jNoOverlap := !hasOverlap(ojs, ltgj, ltgs)
			ojs = shrink(ojs)
			jInS := slices.ContainsFunc(sPObs, func(o obstruction[K]) bool { return obsEq(o, ojs) })
			if !(jNoOverlap || jInS) {
				continue
			}

			b[k].removed = true
			break
		}
	}
}

func hasOverlap[K Field[K]](o obstruction[K], im, jm Monomial) bool {
	if len(o.iLeft)+len(im) <= len(o.jLeft) {
		return false
	}
	if len(im)+len(o.iRight) <= len(o.jRight) {
		return false
	}
	return true
}

func shrink[K Field[K]](o obstruction[K]) obstruction[K] {
	// Shrink left.
	leftLen := min(len(o.iLeft), len(o.jLeft))
	leftStart := leftLen
	for i := range leftLen {
		if o.iLeft[i] != o.jLeft[i] {
			leftStart = i
			break
		}
	}
	o.iLeft = o.iLeft[leftStart:]
	o.jLeft = o.jLeft[leftStart:]

	// Shrink right.
	rightLen := min(len(o.iRight), len(o.jRight))
	rightEnd := rightLen
	for i := range rightLen {
		if o.iRight[len(o.iRight)-1-i] != o.jRight[len(o.jRight)-1-i] {
			rightEnd = i
			break
		}
	}
	o.iRight = o.iRight[:len(o.iRight)-rightEnd]
	o.jRight = o.jRight[:len(o.jRight)-rightEnd]

	return o
}

func overlapObstruction[K Field[K]](obs []obstruction[K], g []*Polynomial[K]) []obstruction[K] {
	sP := len(g) - 1
	for i := range sP {
		obs = leftObstruction(obs, i, sP, g)
		obs = rightObstruction(obs, i, sP, g)
		obs = centerObstruction(obs, i, sP, g)
	}
	obs = rightObstruction(obs, sP, sP, g)
	return obs
}

func leftObstruction[K Field[K]](obs []obstruction[K], i, j int, g []*Polynomial[K]) []obstruction[K] {
	ltgi := g[i].LeadingTerm().Monomial
	ltgj := g[j].LeadingTerm().Monomial
	var iEnd, jStart int
	switch {
	case i == j:
		iEnd, jStart = len(ltgi)-1, 1
	case len(ltgi) < len(ltgj):
		iEnd, jStart = len(ltgi), len(ltgj)-len(ltgi)
	default:
		iEnd, jStart = len(ltgj), 0
	}

	for jStart < len(ltgj) {
		iOverlap := ltgi[:iEnd]
		jOverlap := ltgj[jStart:]
		if monomialEq(iOverlap, jOverlap) {
			o := obstruction[K]{i: i, j: j}
			o.iLeft = ltgj[:jStart]
			o.jRight = ltgi[iEnd:]
			obs = append(obs, o)
		}

		iEnd--
		jStart++
	}

	return obs
}

func rightObstruction[K Field[K]](obs []obstruction[K], i, j int, g []*Polynomial[K]) []obstruction[K] {
	ltgi := g[i].LeadingTerm().Monomial
	ltgj := g[j].LeadingTerm().Monomial
	var iStart, jEnd int
	switch {
	case i == j:
		iStart, jEnd = 1, len(ltgj)-1
	case len(ltgi) < len(ltgj):
		iStart, jEnd = 0, len(ltgi)
	default:
		iStart, jEnd = len(ltgi)-len(ltgj), len(ltgj)
	}

	for iStart < len(ltgi) {
		iOverlap := ltgi[iStart:]
		jOverlap := ltgj[:jEnd]
		if monomialEq(iOverlap, jOverlap) {
			o := obstruction[K]{i: i, j: j}
			o.iRight = ltgj[jEnd:]
			o.jLeft = ltgi[:iStart]
			obs = append(obs, o)
		}

		iStart++
		jEnd--
	}

	return obs
}

func centerObstruction[K Field[K]](obs []obstruction[K], i, j int, g []*Polynomial[K]) []obstruction[K] {
	ltgi := g[i].LeadingTerm().Monomial
	ltgj := g[j].LeadingTerm().Monomial

	if len(ltgi) < len(ltgj) {
		for jStart := 1; jStart <= len(ltgj)-len(ltgi)-1; jStart++ {
			jEnd := jStart + len(ltgi)
			jOverlap := ltgj[jStart:jEnd]
			if monomialEq(ltgi, jOverlap) {
				o := obstruction[K]{i: i, j: j}
				o.iLeft = ltgj[:jStart]
				o.iRight = ltgj[jEnd:]
				obs = append(obs, o)
			}
		}
	} else {
		for iStart := 1; iStart <= len(ltgi)-len(ltgj)-1; iStart++ {
			iEnd := iStart + len(ltgj)
			iOverlap := ltgi[iStart:iEnd]
			if monomialEq(iOverlap, ltgj) {
				o := obstruction[K]{i: i, j: j}
				o.jLeft = ltgi[:iStart]
				o.jRight = ltgi[iEnd:]
				obs = append(obs, o)
			}
		}
	}

	return obs
}

func homogeneous[K Field[K]](x *Polynomial[K]) bool {
	for i := range x.m.Len() - 1 {
		im, _ := x.m.At(i)
		i1m, _ := x.m.At(i + 1)
		if len(im) != len(i1m) {
			return false
		}
	}
	return true
}

func englishSymbolStringer(s Symbol) string {
	return string((s % 26) - 1 + 'a')
}

func printSymbol(b *strings.Builder, s Symbol, power int, symbolStringer func(Symbol) string) {
	v := symbolStringer(s)
	switch {
	case power == 1:
		fmt.Fprintf(b, "%s", v)
	default:
		fmt.Fprintf(b, "%s^%d", v, power)
	}
}

func printMonomial(b *strings.Builder, w Monomial, ss func(Symbol) string) {
	if len(w) == 0 {
		return
	}
	prev, pow := w[0], 1
	for _, s := range w[1:] {
		if s != prev {
			printSymbol(b, prev, pow, ss)
			prev, pow = s, 1
		} else {
			pow++
		}
	}
	printSymbol(b, prev, pow, ss)
}

func monomialEq(x, y Monomial) bool { return slices.Equal(x, y) }

func monomialIndex(x, y Monomial) int {
	return bytes.Index(x, y)
}

func cutSuffix(x, y Monomial) (Monomial, bool) {
	return bytes.CutSuffix(x, y)
}

func cutPrefix(x, y Monomial) (Monomial, bool) {
	return bytes.CutPrefix(x, y)
}

func lexicographic(x, y Monomial) int {
	for i := range x {
		if !(i < len(y)) {
			return 1
		}
		if c := cmp.Compare(x[i], y[i]); c != 0 {
			return c
		}
	}
	if len(x) == len(y) {
		return 0
	}
	return -1
}

func obsEq[K Field[K]](x, y obstruction[K]) bool {
	if x.i != y.i {
		return false
	}
	if x.j != y.j {
		return false
	}
	if !monomialEq(x.iLeft, y.iLeft) {
		return false
	}
	if !monomialEq(x.iRight, y.iRight) {
		return false
	}
	if !monomialEq(x.jLeft, y.jLeft) {
		return false
	}
	if !monomialEq(x.jRight, y.jRight) {
		return false
	}
	return true
}

func polynomialCmp[K Field[K]](x, y *Polynomial[K]) int {
	// Compare monomials.
	for i := range x.m.Len() {
		if i >= y.m.Len() {
			return 1
		}
		xw, _ := x.m.At(x.m.Len() - 1 - i)
		yw, _ := y.m.At(y.m.Len() - 1 - i)
		if wo := x.order(xw, yw); wo != 0 {
			return wo
		}
	}
	if x.m.Len() < y.m.Len() {
		return -1
	}

	// Compare coefficients.
	for i := range x.m.Len() {
		_, xc := x.m.At(x.m.Len() - 1 - i)
		_, yc := y.m.At(y.m.Len() - 1 - i)
		if co := cmp.Compare(xc.String(), yc.String()); co != 0 {
			return co
		}
	}
	return 0
}
