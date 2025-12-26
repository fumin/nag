package field

import (
	"github.com/pkg/errors"

	"github.com/fumin/nag"
)

type factor[K Finite[K]] struct {
	A   *nag.Polynomial[K]
	Pow int
}

func squareFree[K Finite[K]](f *nag.Polynomial[K]) []factor[K] {
	c := gcd(f, differentiate(f))
	w, _ := divide(f, c)

	r := make([]factor[K], 0)
	var i int = 1
	for !isOne(w) {
		y := gcd(w, c)
		fac, _ := divide(w, y)
		if !isOne(fac) {
			r = append(r, factor[K]{A: fac, Pow: i})
		}
		w = y
		c, _ = divide(c, y)
		i++
	}

	if !isOne(c) {
		// Let the field be GF(q), where q = p^e.
		// Compute cRoot = c^{1/p}.
		lt := c.LeadingTerm().Coefficient
		p, e := lt.Characteristic(), lt.PrimePower()
		k := f.Field()
		cRoot := nag.NewPolynomial(k, c.Order())
		cRoot.SymbolStringer = c.SymbolStringer
		for cc, cw := range c.Terms() {
			// x^{ap} -> x^a
			w := make([]nag.Symbol, len(cw)/p)
			copy(w, cw)

			// c -> c^{p^(e-1)}
			pe1 := int(exp(nag.Rationals, nag.NewRat(int64(p), 1), e-1).Num().Int64())
			rc := exp(k, cc, pe1)

			cRoot.Add(cRoot, nag.NewPolynomial(k, c.Order(), nag.PolynomialTerm[K]{Coefficient: rc, Monomial: w}))
		}

		// Do factorization on c^{1/p}.
		sf := squareFree(cRoot)
		for i := range sf {
			sf[i].Pow *= p
		}
		r = append(r, sf...)
	}
	return r
}

func differentiate[K nag.Field[K]](a *nag.Polynomial[K]) *nag.Polynomial[K] {
	k := a.Field()
	aP := nag.NewPolynomial[K](k, a.Order())
	aP.SymbolStringer = a.SymbolStringer
	for ac, aw := range a.Terms() {
		if len(aw) == 0 {
			continue
		}
		w := make([]nag.Symbol, len(aw)-1)
		copy(w, aw)

		pow := k.NewZero()
		for range len(aw) {
			pow.Add(pow, k.NewOne())
		}
		c := k.NewZero().Mul(ac, pow)
		aP.Add(aP, nag.NewPolynomial(k, a.Order(), nag.PolynomialTerm[K]{Coefficient: c, Monomial: w}))
	}
	return aP
}

func gcd[K nag.Field[K]](a, b *nag.Polynomial[K]) *nag.Polynomial[K] {
	r0, r1 := a, b
	for r1.Len() != 0 {
		_, r2 := divide(r0, r1)
		r0, r1 = r1, r2
	}

	// Make result monic.
	c := r0.LeadingTerm().Coefficient
	r0 = mulScalar(r0, r0.Field().NewZero().Inv(c))
	return r0
}

func divide[K nag.Field[K]](a, b *nag.Polynomial[K]) (*nag.Polynomial[K], *nag.Polynomial[K]) {
	k := a.Field()
	a2 := nag.NewPolynomial(k, a.Order()).Set(a)
	quotient := make([][]nag.Quotient[K], 0)
	quotient, r := nag.Divide(quotient, a2, []*nag.Polynomial[K]{b})

	q := nag.NewPolynomial(k, a.Order())
	q.SymbolStringer = r.SymbolStringer
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

func mulScalar[K nag.Field[K]](a *nag.Polynomial[K], b K) *nag.Polynomial[K] {
	bP := nag.NewPolynomial(a.Field(), a.Order(), nag.PolynomialTerm[K]{Coefficient: b})
	ab := nag.NewPolynomial(a.Field(), a.Order()).Mul(a, bP)
	ab.SymbolStringer = a.SymbolStringer
	return ab
}

func mul[K nag.Field[K]](x *nag.Polynomial[K], y ...*nag.Polynomial[K]) *nag.Polynomial[K] {
	z := x
	for i := range y {
		z = nag.NewPolynomial(z.Field(), z.Order()).Mul(z, y[i])
		z.SymbolStringer = x.SymbolStringer
	}
	return z
}

func exp[K nag.Field[K]](k nag.FieldSet[K], x K, n int) K {
	switch {
	case n < 0:
		return exp(k, k.NewOne().Inv(x), -n)
	case n == 0:
		return k.NewOne()
	case n%2 == 0:
		return exp(k, k.NewOne().Mul(x, x), n/2)
	default:
		return k.NewOne().Mul(x, exp(k, k.NewOne().Mul(x, x), (n-1)/2))
	}
}

func parse(order int, s string) (*nag.Polynomial[*Prime], error) {
	variables := map[string]nag.Symbol{"x": 1}
	rp, err := nag.Parse(variables, nag.Deglex, s)
	if err != nil {
		return nil, errors.Wrap(err, "")
	}

	// Cast coefficients from rationals to GF(order).
	field := NewPrimes(order)
	p := nag.NewPolynomial[*Prime](field, rp.Order())
	p.SymbolStringer = rp.SymbolStringer
	for rc, w := range rp.Terms() {
		c := NewPrime(order, int(rc.Num().Int64()))
		p.Add(p, nag.NewPolynomial(p.Field(), p.Order(), nag.PolynomialTerm[*Prime]{Coefficient: c, Monomial: w}))
	}
	return p, nil
}
