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
	"math/big"
	"slices"
	"strings"

	"github.com/jba/omap"
)

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
type PolynomialTerm struct {
	Coefficient *big.Rat
	Monomial    Monomial
}

// A Polynomial is a polynomial of noncommutative variables.
type Polynomial struct {
	// SymbolStringer specifies how a symbol in a monomial is formated when the polynomial is printed out.
	SymbolStringer func(s Symbol) string

	order Order
	m     *omap.MapFunc[Monomial, *big.Rat]

	// Buffers
	r0 *big.Rat
}

// NewPolynomial returns a new polynomial containing the given terms.
func NewPolynomial(order Order, terms ...PolynomialTerm) *Polynomial {
	x := &Polynomial{
		SymbolStringer: englishSymbolStringer,
		order:          order,
		m:              omap.NewMapFunc[Monomial, *big.Rat](order),
		r0:             big.NewRat(0, 1),
	}
	for _, term := range terms {
		if term.Coefficient == nil {
			term.Coefficient = big.NewRat(1, 1)
		}
		x.addTerm(1, term)
	}
	return x
}

// Terms iterates the terms in a polynomial.
func (x *Polynomial) Terms() iter.Seq2[*big.Rat, Monomial] {
	return func(yield func(*big.Rat, Monomial) bool) {
		for w, c := range x.m.Backward() {
			if !yield(c, w) {
				return
			}
		}
	}
}

// Cmp compares two polynomials by first comparing their monomials consecutively in descending monomial order.
// In the case of a tie, Cmp compares their coefficients again consecutively.
func (x *Polynomial) Cmp(y *Polynomial) int {
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
		if co := xc.Cmp(yc); co != 0 {
			return co
		}
	}
	return 0
}

// Set sets z to x and returns z.
func (z *Polynomial) Set(x *Polynomial) *Polynomial {
	z.SymbolStringer = x.SymbolStringer
	z.order = x.order
	z.m = omap.NewMapFunc[Monomial, *big.Rat](z.order)
	for xw, xc := range x.m.All() {
		w := make(Monomial, len(xw))
		copy(w, xw)
		z.addTerm(1, PolynomialTerm{Coefficient: xc, Monomial: w})
	}
	return z
}

// Add sets z to the sum x+y and returns z.
func (z *Polynomial) Add(x, y *Polynomial) *Polynomial {
	// Set z = x, while handling the case where x or y is z itself.
	if y == z {
		x, y = y, x
	}
	if z != x {
		z.m.Clear()
		for xw, c := range x.m.All() {
			w := make(Monomial, len(xw))
			copy(w, xw)
			z.addTerm(1, PolynomialTerm{Coefficient: c, Monomial: w})
		}
	}

	// Compute z += y.
	for yw, c := range y.m.All() {
		w := make(Monomial, len(yw))
		copy(w, yw)
		z.addTerm(1, PolynomialTerm{Coefficient: c, Monomial: w})
	}

	return z
}

// Mul sets z to the product x*y and returns z.
func (z *Polynomial) Mul(x, y *Polynomial) *Polynomial {
	if z == x {
		panic(fmt.Sprintf("z == x"))
	}
	if z == y {
		panic(fmt.Sprintf("z == y"))
	}

	z.m.Clear()
	for xw, xc := range x.m.Backward() {
		for yw, yc := range y.m.Backward() {
			c := z.r0.Mul(xc, yc)
			w := make(Monomial, 0, len(xw)+len(yw))
			w = append(append(w, xw...), yw...)
			z.addTerm(1, PolynomialTerm{Coefficient: c, Monomial: w})
		}
	}

	return z
}

// Pow sets z to the power x^y and returns z.
func (z *Polynomial) Pow(x *Polynomial, y int) *Polynomial {
	if z == x {
		panic("z == x")
	}

	z.Set(x)
	buf := NewPolynomial(z.order)
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
func (x *Polynomial) LeadingTerm() PolynomialTerm {
	w, ok := x.m.Max()
	if !ok {
		return PolynomialTerm{}
	}
	c, _ := x.m.Get(w)
	return PolynomialTerm{Coefficient: c, Monomial: w}
}

// String returns the string representation of x.
// Symbols in x are formatted using x.SymbolStringer.
func (x *Polynomial) String() string {
	var b strings.Builder
	var i int = -1
	for w, c := range x.m.Backward() {
		// Print c.
		i++
		switch {
		case len(w) != 0 && c.Cmp(big.NewRat(1, 1)) == 0:
			if i > 0 {
				fmt.Fprintf(&b, "+")
			}
		case len(w) != 0 && c.Cmp(big.NewRat(-1, 1)) == 0:
			fmt.Fprintf(&b, "-")
		default:
			if c.Sign() == 1 && i > 0 {
				fmt.Fprintf(&b, "+")
			}
			fmt.Fprintf(&b, "%s", c.RatString())
		}

		printMonomial(&b, w, x.SymbolStringer)
	}
	return b.String()
}

func (x *Polynomial) addTerm(sign int, term PolynomialTerm) {
	c, ok := x.m.Get(term.Monomial)
	if !ok {
		c = big.NewRat(0, 1)
	}

	if sign < 0 {
		c.Sub(c, term.Coefficient)
	} else {
		c.Add(c, term.Coefficient)
	}

	if c.Sign() == 0 {
		x.m.Delete(term.Monomial)
	} else {
		x.m.Set(term.Monomial, c)
	}
}

func (z *Polynomial) add(sign int, c *big.Rat, left Monomial, x *Polynomial, right Monomial) {
	for xw, xc := range x.m.Backward() {
		c := z.r0.Mul(c, xc)
		w := make(Monomial, 0, len(left)+len(xw)+len(right))
		w = append(append(append(w, left...), xw...), right...)
		z.addTerm(sign, PolynomialTerm{Coefficient: c, Monomial: w})
	}
}

func (z *Polynomial) mulScalar(scalar *big.Rat, x *Polynomial) *Polynomial {
	if z == x {
		for zw, zc := range z.m.All() {
			zc.Mul(scalar, zc)
			z.m.Set(zw, zc)
		}
		return z
	}

	z.m.Clear()
	for xw, xc := range x.m.All() {
		c := z.r0.Mul(scalar, xc)
		w := make(Monomial, len(xw))
		copy(w, xw)
		z.addTerm(1, PolynomialTerm{Coefficient: c, Monomial: w})
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
type Quotient struct {
	// Coefficient is c_{ij} in Theorem 3.2.1, Xiu Xingqiang.
	Coefficient *big.Rat
	// Coefficient is w_{ij} in Theorem 3.2.1, Xiu Xingqiang.
	Left Monomial
	// Coefficient is w'_{ij} in Theorem 3.2.1, Xiu Xingqiang.
	Right Monomial
}

// Divide divides the polynomial f by the ideal g, and returns the remainder polynomial and quotient.
// The polynomial f is modified upon return.
// For more details, please see Theorem 3.2.1, Xiu Xingqiang.
//
// Xiu, Xingqiang. "Non-commutative Gröbner bases and applications." PhD diss., Universität Passau, 2012.
func Divide(quotient [][]Quotient, f *Polynomial, g []*Polynomial) (remainder *Polynomial, outQuotient [][]Quotient) {
	if quotient != nil {
		short := len(g) - len(quotient)
		if short > 0 {
			quotient = append(quotient, make([][]Quotient, short)...)
		}
		quotient = quotient[:len(g)]
		for i := range quotient {
			quotient[i] = quotient[i][:0]
		}
	}
	p := NewPolynomial(f.order)
	p.SymbolStringer = f.SymbolStringer
	v := f

	for v.m.Len() != 0 {
		lmv := v.LeadingTerm()
		ltv := lmv.Monomial

		// Find basis where ltv = left * ltg * right.
		basis, leftEnd, lmg := -1, -1, PolynomialTerm{}
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
			q := Quotient{
				Coefficient: big.NewRat(0, 1).Quo(lmv.Coefficient, lmg.Coefficient),
				Left:        ltv[:leftEnd],
				Right:       ltv[leftEnd+len(lmg.Monomial):],
			}
			if quotient != nil {
				quotient[basis] = append(quotient[basis], q)
			}
			v.add(-1, q.Coefficient, q.Left, g[basis], q.Right)
		}
	}

	return p, quotient
}

// Buchberger returns the Gröbner basis of an ideal g, using the Buchberger algorithm.
// In the noncummutative case, a Gröbner basis may not be finite and therefore may need to be truncated.
// Buchberger truncates the returned basis when maxiter is reached, and returns complete as false.
// For more details, please see Theorem 4.2.24, Xiu Xingqiang.
//
// Xiu, Xingqiang. "Non-commutative Gröbner bases and applications." PhD diss., Universität Passau, 2012.
func Buchberger(g []*Polynomial, maxiter int) (outG []*Polynomial, complete bool) {
	// Buffers.
	r0 := big.NewRat(0, 1)
	buf := &Monomial{}

	g = interreduce(g)
	// t tracks unwanted polynomials in g.
	t := make([]*Polynomial, len(g))
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
	var b []obstruction
	for l := 1; l <= len(g); l++ {
		gl := g[:l]
		b = addObstructions(b, gl, buf)
	}

	for range maxiter {
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
		gi, gj := g[oij.i], g[oij.j]
		lcgi := gi.LeadingTerm().Coefficient
		lcgj := gj.LeadingTerm().Coefficient
		s := NewPolynomial(g[0].order)
		s.SymbolStringer = g[0].SymbolStringer
		s.add(1, r0.Inv(lcgi), oij.iLeft, gi, oij.iRight)
		s.add(-1, r0.Inv(lcgj), oij.jLeft, gj, oij.jRight)
		deleteUnwanted()
		sP, _ := Divide(nil, s, g)
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
	g = slices.DeleteFunc(g, func(gi *Polynomial) bool { return gi == nil })
	// Make basis monic.
	g = interreduce(g)
	for i := range g {
		lc := g[i].LeadingTerm().Coefficient
		g[i].mulScalar(r0.Inv(lc), g[i])
	}
	slices.SortFunc(g, func(a, b *Polynomial) int { return a.Cmp(b) })
	return g, complete
}

func interreduce(g []*Polynomial) []*Polynomial {
	i, s := 0, len(g)
	for i != s {
		gi := g[i]
		if gi == nil {
			i++
			continue
		}
		f := NewPolynomial(Deglex).Set(gi)
		g[i] = nil
		giP, _ := Divide(nil, f, g)

		switch {
		case giP.m.Len() == 0:
			g[i] = nil
			i++
		case giP.Cmp(gi) != 0:
			g[i] = giP
			i = 0
		default:
			g[i] = gi
			i++
		}
	}
	return slices.DeleteFunc(g, func(x *Polynomial) bool { return x == nil })
}

type obstruction struct {
	i      int
	j      int
	iLeft  Monomial
	iRight Monomial
	jLeft  Monomial
	jRight Monomial

	removed bool
}

func delRemoved(sPObs, obs []obstruction) ([]obstruction, []obstruction) {
	spPrevLen := len(sPObs)
	sPObs = slices.DeleteFunc(sPObs, func(o obstruction) bool { return o.removed })
	obs = obs[:len(obs)-(spPrevLen-len(sPObs))]
	return sPObs, obs
}

func addObstructions(obs []obstruction, g []*Polynomial, buf *Monomial) []obstruction {
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
	obs = slices.DeleteFunc(obs, func(o obstruction) bool { return o.removed })

	return obs
}

func remove4b(sPObs []obstruction, order Order) {
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

func remove4c(sPObs, b []obstruction, ltgs Monomial) {
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

func remove4d(b, sPObs []obstruction, g []*Polynomial, buf *Monomial) {
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
			ois := obstruction{i: oij.i, j: sP, iLeft: oij.iLeft, iRight: oij.iRight, jLeft: ws, jRight: wsP}
			iNoOverlap := !hasOverlap(ois, ltgi, ltgs)
			ois = shrink(ois)
			iInS := slices.ContainsFunc(sPObs, func(o obstruction) bool { return obsEq(o, ois) })
			if !(iNoOverlap || iInS) {
				continue
			}

			// Check o(j, s).
			ojs := obstruction{i: oij.j, j: sP, iLeft: oij.jLeft, iRight: oij.jRight, jLeft: ws, jRight: wsP}
			jNoOverlap := !hasOverlap(ojs, ltgj, ltgs)
			ojs = shrink(ojs)
			jInS := slices.ContainsFunc(sPObs, func(o obstruction) bool { return obsEq(o, ojs) })
			if !(jNoOverlap || jInS) {
				continue
			}

			b[k].removed = true
			break
		}
	}
}

func hasOverlap(o obstruction, im, jm Monomial) bool {
	if len(o.iLeft)+len(im) <= len(o.jLeft) {
		return false
	}
	if len(im)+len(o.iRight) <= len(o.jRight) {
		return false
	}
	return true
}

func shrink(o obstruction) obstruction {
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

func overlapObstruction(obs []obstruction, g []*Polynomial) []obstruction {
	sP := len(g) - 1
	for i := range sP {
		obs = leftObstruction(obs, i, sP, g)
		obs = rightObstruction(obs, i, sP, g)
		obs = centerObstruction(obs, i, sP, g)
	}
	obs = rightObstruction(obs, sP, sP, g)
	return obs
}

func leftObstruction(obs []obstruction, i, j int, g []*Polynomial) []obstruction {
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
			o := obstruction{i: i, j: j}
			o.iLeft = ltgj[:jStart]
			o.jRight = ltgi[iEnd:]
			obs = append(obs, o)
		}

		iEnd--
		jStart++
	}

	return obs
}

func rightObstruction(obs []obstruction, i, j int, g []*Polynomial) []obstruction {
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
			o := obstruction{i: i, j: j}
			o.iRight = ltgj[jEnd:]
			o.jLeft = ltgi[:iStart]
			obs = append(obs, o)
		}

		iStart++
		jEnd--
	}

	return obs
}

func centerObstruction(obs []obstruction, i, j int, g []*Polynomial) []obstruction {
	ltgi := g[i].LeadingTerm().Monomial
	ltgj := g[j].LeadingTerm().Monomial

	if len(ltgi) < len(ltgj) {
		for jStart := 1; jStart <= len(ltgj)-len(ltgi)-1; jStart++ {
			jEnd := jStart + len(ltgi)
			jOverlap := ltgj[jStart:jEnd]
			if monomialEq(ltgi, jOverlap) {
				o := obstruction{i: i, j: j}
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
				o := obstruction{i: i, j: j}
				o.jLeft = ltgi[:iStart]
				o.jRight = ltgi[iEnd:]
				obs = append(obs, o)
			}
		}
	}

	return obs
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

func obsEq(x, y obstruction) bool {
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
