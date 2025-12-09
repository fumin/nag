package nag

import (
	"bytes"
	"cmp"
	"fmt"
	"math/big"
	"slices"
	"strings"
	"unsafe"

	"github.com/jba/omap"
)

type Symbol byte

type Monomial []Symbol

type Ordering func(a, b Monomial) int

func Deglex(a, b Monomial) int {
	if c := cmp.Compare(len(a), len(b)); c != 0 {
		return c
	}
	return lexicographic(a, b)
}

type PolynomialTerm struct {
	Coefficient *big.Rat
	Monomial    Monomial
}

type Polynomial struct {
	SymbolStringer func(s Symbol) string

	ordering Ordering
	m        *omap.MapFunc[Monomial, *big.Rat]

	// Buffers
	r0 *big.Rat
}

func NewPolynomial(ordering Ordering, terms ...PolynomialTerm) *Polynomial {
	x := &Polynomial{
		SymbolStringer: englishSymbolStringer,
		ordering:       ordering,
		m:              omap.NewMapFunc[Monomial, *big.Rat](ordering),
		r0:             big.NewRat(0, 1),
	}
	for _, term := range terms {
		x.addTerm(1, term)
	}
	return x
}

func (x *Polynomial) Cmp(y *Polynomial) int {
	for i := range x.m.Len() {
		if !(i < y.m.Len()) {
			return 1
		}
		xw, xc := x.m.At(i)
		yw, yc := y.m.At(i)
		if wo := x.ordering(xw, yw); wo != 0 {
			return wo
		}
		if co := xc.Cmp(yc); co != 0 {
			return co
		}
	}
	if x.m.Len() == y.m.Len() {
		return 0
	}
	return -1
}

func (z *Polynomial) Set(x *Polynomial) *Polynomial {
	z.SymbolStringer = x.SymbolStringer
	z.ordering = x.ordering
	z.m = omap.NewMapFunc[Monomial, *big.Rat](z.ordering)
	for xw, xc := range x.m.All() {
		w := make(Monomial, len(xw))
		copy(w, xw)
		z.addTerm(1, PolynomialTerm{Coefficient: xc, Monomial: w})
	}
	return z
}

func (z *Polynomial) Add(x, y *Polynomial) *Polynomial {
	// Set z = x, while handling the case where x or y is z itself.
	if y == z {
		x, y = y, x
	}
	if z != x {
		z.m.Clear()
		for w, c := range x.m.All() {
			z.addTerm(1, PolynomialTerm{Coefficient: c, Monomial: w})
		}
	}

	// Compute z += y.
	for w, c := range y.m.All() {
		z.addTerm(1, PolynomialTerm{Coefficient: c, Monomial: w})
	}

	return z
}

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

func (x *Polynomial) LeadingTerm() PolynomialTerm {
	w, ok := x.m.Max()
	if !ok {
		return PolynomialTerm{}
	}
	c, _ := x.m.Get(w)
	return PolynomialTerm{Coefficient: c, Monomial: w}
}

func (x *Polynomial) String() string {
	var b strings.Builder
	var i int = -1
	for w, c := range x.m.Backward() {
		// Print c.
		i++
		if c.Sign() == 1 && i > 0 {
			fmt.Fprintf(&b, "+")
		}
		fmt.Fprintf(&b, "%v", c)

		// Print w.
		for _, s := range w {
			b.WriteString(x.SymbolStringer(s))
		}
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

type Quotient struct {
	Coefficient *big.Rat
	Left        Monomial
	Right       Monomial
}

// Theorem 3.2.1
// f is modified.
func Divide(quotient [][]Quotient, f *Polynomial, g []*Polynomial) (*Polynomial, [][]Quotient) {
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
	p := NewPolynomial(f.ordering)
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

func interreduce(g []*Polynomial) []*Polynomial {
	i, s := 0, len(g)
	for i != s {
		gi := g[i]
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
	remove4b(sPObs, g[0].ordering)
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

func remove4b(sPObs []obstruction, ordering Ordering) {
	for k := range sPObs {
		for l := range sPObs {
			if l == k {
				continue
			}
			oi, oj := sPObs[k], sPObs[l]

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
			wiLtwj := (ordering(oi.iLeft, oj.iLeft) == 1)

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

func monomialEq(x, y Monomial) bool { return slices.Equal(x, y) }

func monomialIndex(x, y Monomial) int {
	xb := *(*[]byte)(unsafe.Pointer(&x))
	yb := *(*[]byte)(unsafe.Pointer(&y))
	return bytes.Index(xb, yb)
}

func cutSuffix(x, y Monomial) (Monomial, bool) {
	xb := *(*[]byte)(unsafe.Pointer(&x))
	yb := *(*[]byte)(unsafe.Pointer(&y))
	beforeb, ok := bytes.CutSuffix(xb, yb)
	before := *(*Monomial)(unsafe.Pointer(&beforeb))
	return before, ok
}

func cutPrefix(x, y Monomial) (Monomial, bool) {
	xb := *(*[]byte)(unsafe.Pointer(&x))
	yb := *(*[]byte)(unsafe.Pointer(&y))
	afterb, ok := bytes.CutPrefix(xb, yb)
	after := *(*Monomial)(unsafe.Pointer(&afterb))
	return after, ok
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
