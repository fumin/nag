package field

import (
	"cmp"
	"math/rand"
	"slices"

	"github.com/pkg/errors"

	"github.com/fumin/nag"
)

// A finite is an element of a finite field GF(q) with order q = p^n where p is a prime number.
type finite[K nag.Field[K]] interface {
	nag.Field[K]
	// characteristic returns p in the finite field GF(p^n).
	characteristic() int
	// primePower returns n in the finite field GF(p^n).
	primePower() int
}

func primitive[K finite[K]](f *nag.Polynomial[K]) bool {
	m := len(f.LeadingTerm().Monomial)
	k := f.Field()
	p := k.characteristic()
	n := expi(p, m) - 1
	for i := 1; i < n; i++ {
		g := nag.NewPolynomial(k, f.Order(), nag.PolynomialTerm[K]{Coefficient: k.NewOne(), Monomial: make([]nag.Symbol, i)})
		g.SymbolStringer = func(x nag.Symbol) string { return "x" }
		g = sub(g, poly1(g))
		_, r := divide(g, f)
		if r.Len() == 0 {
			return false
		}
	}
	return true
}

type factor[K finite[K]] struct {
	A *nag.Polynomial[K]
	N int
}

func factorize[K finite[K]](f *nag.Polynomial[K]) []factor[K] {
	factors := make([]factor[K], 0)
	sf := squareFree(f)
	for _, sfi := range sf {
		dd := distinctDegree(sfi.A)
		for _, ddi := range dd {
			ed := equalDegree(ddi.A, ddi.N)
			for _, edi := range ed {
				factors = append(factors, factor[K]{A: edi, N: sfi.N})
			}
		}
	}
	return factors
}

func equalDegree[K finite[K]](f *nag.Polynomial[K], d int) []*nag.Polynomial[K] {
	n := len(f.LeadingTerm().Monomial)
	r := n / d
	k := f.Field()
	// q is the order of the field k.
	q := expi(k.characteristic(), k.primePower())
	// numFq is the number of GF(q) polynomials with degree less than n.
	numFq := expi(q, n)

	factors := []*nag.Polynomial[K]{f}
	for len(factors) < r {
		hi := rand.Intn(numFq-q) + q
		h := ithPoly(poly0(f), hi, q)
		g, _, _, _, _ := gcd(h, f)
		if isOne(g) {
			qd := (expi(q, d) - 1) / 2
			hqd := sub(pow(h, qd), poly1(h))
			_, g = divide(hqd, f)
		}

		checked := make(map[string]struct{})
		for {
			// Find u in factors where deg(u) > d.
			var ui int = -1
			for i, f := range factors {
				if _, ok := checked[f.String()]; ok {
					continue
				}
				if len(f.LeadingTerm().Monomial) > d {
					ui = i
					break
				}
			}
			if ui == -1 {
				break
			}
			u := factors[ui]

			gcdu, _, _, _, _ := gcd(g, u)
			if !isOne(gcdu) && !gcdu.Equal(u) {
				factors = slices.Delete(factors, ui, ui+1)
				factors = append(factors, gcdu)
				ugcd, _ := divide(u, gcdu)
				factors = append(factors, ugcd)
			}
			checked[u.String()] = struct{}{}
		}
	}

	return factors
}

func distinctDegree[K finite[K]](f *nag.Polynomial[K]) []factor[K] {
	degF := len(f.LeadingTerm().Monomial)
	k := f.Field()
	// q is the order of the field k.
	q := expi(k.characteristic(), k.primePower())
	neg1 := k.NewZero()
	neg1.Sub(k.NewZero(), k.NewOne())

	fs := poly0(f).Set(f)
	s := make([]factor[K], 0)
	for i := 1; i <= degF/2+1; i++ {
		qi := expi(q, i)
		xqi := nag.NewPolynomial(k, fs.Order(),
			nag.PolynomialTerm[K]{Coefficient: k.NewOne(), Monomial: make([]nag.Symbol, qi)},
			nag.PolynomialTerm[K]{Coefficient: neg1, Monomial: make([]nag.Symbol, 1)})
		_, xqi = divide(xqi, fs)
		g, _, _, _, _ := gcd(fs, xqi)
		if !isOne(g) {
			s = append(s, factor[K]{A: g, N: i})
			fs, _ = divide(fs, g)
		}
	}
	if !isOne(fs) {
		s = append(s, factor[K]{A: fs, N: degF})
	}
	if len(s) == 0 {
		s = append(s, factor[K]{A: f, N: 1})
	}
	return s
}

func squareFree[K finite[K]](f *nag.Polynomial[K]) []factor[K] {
	c, _, _, _, _ := gcd(f, differentiate(f))
	w, _ := divide(f, c)

	r := make([]factor[K], 0)
	var i int = 1
	for !isOne(w) {
		y, _, _, _, _ := gcd(w, c)
		fac, _ := divide(w, y)
		if !isOne(fac) {
			r = append(r, factor[K]{A: fac, N: i})
		}
		w = y
		c, _ = divide(c, y)
		i++
	}

	if !isOne(c) {
		// Let the field be GF(q), where q = p^e.
		// Compute cRoot = c^{1/p}.
		k := f.Field()
		p, e := k.characteristic(), k.primePower()
		cRoot := poly0(c)
		for cc, cw := range c.Terms() {
			// x^{ap} -> x^a
			w := make([]nag.Symbol, len(cw)/p)
			copy(w, cw)

			// c -> c^{p^(e-1)}
			pe1 := expi(p, e-1)
			rc := exp(cc, pe1)

			cRoot.Add(cRoot, nag.NewPolynomial(k, c.Order(), nag.PolynomialTerm[K]{Coefficient: rc, Monomial: w}))
		}

		// Do factorization on c^{1/p}.
		sf := squareFree(cRoot)
		for i := range sf {
			sf[i].N *= p
		}
		r = append(r, sf...)
	}
	return r
}

func differentiate[K nag.Field[K]](a *nag.Polynomial[K]) *nag.Polynomial[K] {
	k := a.Field()
	aP := poly0(a)
	for ac, aw := range a.Terms() {
		if len(aw) == 0 {
			continue
		}
		w := make([]nag.Symbol, len(aw)-1)
		copy(w, aw)

		c := k.NewZero()
		for range len(aw) {
			c.Add(c, ac)
		}

		aP.Add(aP, nag.NewPolynomial(k, a.Order(), nag.PolynomialTerm[K]{Coefficient: c, Monomial: w}))
	}
	return aP
}

func inverse[K nag.Field[K]](a, p *nag.Polynomial[K]) *nag.Polynomial[K] {
	r, _, v, _, _ := gcd(p, a)
	if !isOne(r) {
		return nil
	}
	return v
}

func gcd[K nag.Field[K]](a, b *nag.Polynomial[K]) (g, u, v, a1, b1 *nag.Polynomial[K]) {
	r0, r1 := a, b
	s0, s1 := poly1(a), poly0(a)
	t0, t1 := poly0(a), poly1(a)

	var n1i int = 1
	for r1.Len() != 0 {
		n1i *= -1
		q, _ := divide(r0, r1)
		r2 := sub(r0, mul(q, r1))
		s2 := sub(s0, mul(q, s1))
		t2 := sub(t0, mul(q, t1))

		r0, r1 = r1, r2
		s0, s1 = s1, s2
		t0, t1 = t1, t2
	}

	k := a.Field()
	var n1, nn1 K
	if n1i == 1 {
		n1 = k.NewOne()
		nn1 = k.NewZero().Sub(k.NewZero(), k.NewOne())
	} else {
		n1 = k.NewZero().Sub(k.NewZero(), k.NewOne())
		nn1 = k.NewOne()
	}
	a1 = mulScalar(t1, n1)
	b1 = mulScalar(s1, nn1)

	// Make g monic.
	c := r0.LeadingTerm().Coefficient
	g = mulScalar(r0, k.NewZero().Inv(c))
	u = mulScalar(s0, k.NewZero().Inv(c))
	v = mulScalar(t0, k.NewZero().Inv(c))
	a1 = mulScalar(a1, c)
	b1 = mulScalar(b1, c)
	return g, u, v, a1, b1
}

func divide[K nag.Field[K]](a, b *nag.Polynomial[K]) (*nag.Polynomial[K], *nag.Polynomial[K]) {
	k := a.Field()
	a2 := poly0(a).Set(a)
	quotient := make([][]nag.Quotient[K], 0)
	quotient, r := nag.Divide(quotient, a2, []*nag.Polynomial[K]{b})

	q := poly0(a)
	for i := range quotient {
		for j := range quotient[i] {
			c := nag.NewPolynomial(k, q.Order(), nag.PolynomialTerm[K]{Coefficient: quotient[i][j].Coefficient})
			left := nag.NewPolynomial(k, q.Order(), nag.PolynomialTerm[K]{Coefficient: k.NewOne(), Monomial: quotient[i][j].Left})
			right := nag.NewPolynomial(k, q.Order(), nag.PolynomialTerm[K]{Coefficient: k.NewOne(), Monomial: quotient[i][j].Right})
			cwgw := mul[K](c, left, right)
			q.Add(q, cwgw)
		}
	}
	return q, r
}

func ithPoly[K nag.Field[K]](poly *nag.Polynomial[K], ith, base int) *nag.Polynomial[K] {
	k := poly.Field()
	var pow int = -1
	for ith != 0 {
		pow++
		var r int
		ith, r = ith/base, ith%base

		c := k.NewZero()
		for range r {
			c.Add(c, k.NewOne())
		}
		w := make([]nag.Symbol, pow)
		poly.Add(poly, nag.NewPolynomial(k, poly.Order(), nag.PolynomialTerm[K]{Coefficient: c, Monomial: w}))
	}
	return poly
}

func isOne[K nag.Field[K]](a *nag.Polynomial[K]) bool {
	if a.Len() != 1 {
		return false
	}
	lt := a.LeadingTerm()
	if len(lt.Monomial) != 0 {
		return false
	}
	if !lt.Coefficient.Equal(a.Field().NewOne()) {
		return false
	}
	return true
}

func pow[K nag.Field[K]](a *nag.Polynomial[K], n int) *nag.Polynomial[K] {
	y := poly1(a)
	for range n {
		y = mul(y, a)
	}
	return y
}

func sub[K nag.Field[K]](x, y *nag.Polynomial[K]) *nag.Polynomial[K] {
	k := x.Field()
	neg1 := k.NewZero()
	neg1.Sub(neg1, k.NewOne())
	negY := mulScalar(y, neg1)
	return poly0(x).Add(x, negY)
}

func mulScalar[K nag.Field[K]](a *nag.Polynomial[K], b K) *nag.Polynomial[K] {
	bP := nag.NewPolynomial(a.Field(), a.Order(), nag.PolynomialTerm[K]{Coefficient: b})
	ab := poly0(a).Mul(a, bP)
	ab.SymbolStringer = a.SymbolStringer
	return ab
}

func mul[K nag.Field[K]](x *nag.Polynomial[K], y ...*nag.Polynomial[K]) *nag.Polynomial[K] {
	z := x
	for i := range y {
		z = poly0(z).Mul(z, y[i])
		z.SymbolStringer = x.SymbolStringer
	}
	return z
}

func expi(x, n int) int {
	return int(exp(nag.NewRat(int64(x), 1), n).Num().Int64())
}

func exp[K nag.Field[K]](x K, n int) K {
	switch {
	case n < 0:
		return exp(x.NewOne().Inv(x), -n)
	case n == 0:
		return x.NewOne()
	case n%2 == 0:
		return exp(x.NewOne().Mul(x, x), n/2)
	default:
		return x.NewOne().Mul(x, exp(x.NewOne().Mul(x, x), (n-1)/2))
	}
}

func poly1[K nag.Field[K]](x *nag.Polynomial[K]) *nag.Polynomial[K] {
	y := nag.NewPolynomial(x.Field(), x.Order(), nag.PolynomialTerm[K]{Coefficient: x.Field().NewOne()})
	y.SymbolStringer = x.SymbolStringer
	return y
}

func poly0[K nag.Field[K]](x *nag.Polynomial[K]) *nag.Polynomial[K] {
	y := nag.NewPolynomial(x.Field(), x.Order())
	y.SymbolStringer = x.SymbolStringer
	return y
}

func polynomialCmp[K nag.Field[K]](x, y *nag.Polynomial[K]) int {
	xTerms := make([]nag.PolynomialTerm[K], 0)
	for c, w := range x.Terms() {
		xTerms = append(xTerms, nag.PolynomialTerm[K]{Coefficient: c, Monomial: w})
	}
	yTerms := make([]nag.PolynomialTerm[K], 0)
	for c, w := range y.Terms() {
		yTerms = append(yTerms, nag.PolynomialTerm[K]{Coefficient: c, Monomial: w})
	}

	// Compare monomials.
	for i := range xTerms {
		if i >= len(yTerms) {
			return 1
		}
		xw := xTerms[i].Monomial
		yw := yTerms[i].Monomial
		if wo := x.Order()(xw, yw); wo != 0 {
			return wo
		}
	}
	if len(xTerms) < len(yTerms) {
		return -1
	}

	// Compare coefficients.
	for i := range xTerms {
		xc := xTerms[i].Coefficient
		yc := yTerms[i].Coefficient
		if co := cmp.Compare(xc.String(), yc.String()); co != 0 {
			return co
		}
	}
	return 0
}

func parse(order int, s string) (*nag.Polynomial[*prime], error) {
	variables := map[string]nag.Symbol{"x": 0}
	rp, err := nag.Parse(variables, nag.Deglex, s)
	if err != nil {
		return nil, errors.Wrap(err, "")
	}

	// Cast coefficients from rationals to GF(order).
	field := newPrime(order, 0)
	p := nag.NewPolynomial[*prime](field, rp.Order())
	p.SymbolStringer = rp.SymbolStringer
	for rc, w := range rp.Terms() {
		c := newPrime(order, int(rc.Num().Int64()))
		p.Add(p, nag.NewPolynomial(p.Field(), p.Order(), nag.PolynomialTerm[*prime]{Coefficient: c, Monomial: w}))
	}
	return p, nil
}
