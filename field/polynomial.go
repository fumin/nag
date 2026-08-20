package field

import (
	"cmp"
	"crypto/rand"
	"fmt"
	"math/big"
	"slices"

	"github.com/fumin/nag"
	"github.com/pkg/errors"
)

func primitive[K Finite[K]](f *nag.Polynomial[K]) bool {
	m := len(f.LeadingTerm().Monomial)
	k := f.Field()
	p := k.Characteristic()
	q := new(big.Int).Exp(p, big.NewInt(int64(m)), nil)
	n := int(new(big.Int).Sub(q, big.NewInt(1)).Int64())
	// Let a be the primitive element we are looking for, then the n-1 non-zero elements are:
	// a, a^2, ..., a^(n-2), a^(n-1) == 1.
	// Therefore, we need only scan through a to a^(n-2).
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

// An IrrFactor is a irreducible polynomial P raised to the power of N.
type IrrFactor[K Finite[K]] struct {
	// P is the irreducible polynomial.
	P *nag.Polynomial[K]
	// N is the exponent.
	N *big.Int
}

// Factor factors f into irreducible polynomials.
func Factor[K Finite[K]](f *nag.Polynomial[K]) []IrrFactor[K] {
	factors := make([]IrrFactor[K], 0)
	sf := squareFree(f)
	for _, sfi := range sf {
		dd := distinctDegree(sfi.P)
		for _, ddi := range dd {
			ed := equalDegree(ddi.P, ddi.N)
			for _, edi := range ed {
				factors = append(factors, IrrFactor[K]{P: edi, N: sfi.N})
			}
		}
	}
	return factors
}

// Factoring Polynomials Over Finite Fields: A Survey, Joachim von zur Gathen, Daniel Panario.
func equalDegree[K Finite[K]](f *nag.Polynomial[K], d *big.Int) []*nag.Polynomial[K] {
	n := big.NewInt(int64(len(f.LeadingTerm().Monomial)))
	r := new(big.Int).Div(n, d)
	k := f.Field()
	// q is the order of the field k.
	q := new(big.Int).Exp(k.Characteristic(), big.NewInt(int64(k.Degree())), nil)

	factors := []*nag.Polynomial[K]{f}
	for uint64(len(factors)) < r.Uint64() {
		h := randPoly(poly0(f), n)
		g, _, _, _, _ := gcd(h, f)
		if isOne(g) {
			qd := new(big.Int).Exp(q, d, nil)
			qdExp := new(big.Int).Div(new(big.Int).Sub(qd, big.NewInt(1)), big.NewInt(2))
			hqd := sub(nag.Pow(newPolynomialMod(h, f), qdExp).p, poly1(h))
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
				if deg(f).Cmp(d) > 0 {
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

func distinctDegree[K Finite[K]](f *nag.Polynomial[K]) []IrrFactor[K] {
	k := f.Field()
	// q is the order of the field k.
	q := new(big.Int).Exp(k.Characteristic(), big.NewInt(int64(k.Degree())), nil)
	neg1 := k.NewZero()
	neg1.Sub(k.NewZero(), k.NewOne())

	fs := poly0(f).Set(f)
	s := make([]IrrFactor[K], 0)
	for i := big.NewInt(1); deg(fs).Cmp(new(big.Int).Mul(big.NewInt(2), i)) >= 0; i = new(big.Int).Add(i, big.NewInt(1)) {
		qi := new(big.Int).Exp(q, i, nil)
		x := newPolynomialMod(nag.NewPolynomial(k, fs.Order(), nag.PolynomialTerm[K]{Coefficient: k.NewOne(), Monomial: make(nag.Monomial, 1)}), fs)
		xqi := nag.Pow(x, qi).p
		xqi.Add(xqi, nag.NewPolynomial(k, fs.Order(), nag.PolynomialTerm[K]{Coefficient: neg1, Monomial: make(nag.Monomial, 1)}))
		g, _, _, _, _ := gcd(fs, xqi)
		if !isOne(g) {
			s = append(s, IrrFactor[K]{P: g, N: new(big.Int).Set(i)})
			fs, _ = divide(fs, g)
		}
	}
	if !isOne(fs) {
		s = append(s, IrrFactor[K]{P: fs, N: deg(fs)})
	}
	if len(s) == 0 {
		s = append(s, IrrFactor[K]{P: f, N: big.NewInt(1)})
	}
	return s
}

func squareFree[K Finite[K]](f *nag.Polynomial[K]) []IrrFactor[K] {
	c, _, _, _, _ := gcd(f, differentiate(f))
	w, _ := divide(f, c)

	r := make([]IrrFactor[K], 0)
	i := big.NewInt(1)
	for !isOne(w) {
		y, _, _, _, _ := gcd(w, c)
		fac, _ := divide(w, y)
		if !isOne(fac) {
			r = append(r, IrrFactor[K]{P: fac, N: new(big.Int).Set(i)})
		}
		w = y
		c, _ = divide(c, y)
		i.Add(i, big.NewInt(1))
	}

	if !isOne(c) {
		// Let the field be GF(q), where q = p^e.
		// Compute cRoot = c^{1/p}.
		k := f.Field()
		p, e := k.Characteristic(), big.NewInt(int64(k.Degree()))
		cRoot := poly0(c)
		for cc, cw := range c.Terms() {
			// x^{ap} -> x^a
			wLen := new(big.Int).Div(big.NewInt(int64(len(cw))), p)
			w := make(nag.Monomial, wLen.Int64())
			copy(w, cw)

			// c -> c^{p^(e-1)}
			e1 := new(big.Int).Sub(e, big.NewInt(1))
			pe1 := new(big.Int).Exp(p, e1, nil)
			rc := nag.Pow(cc.NewZero().Set(cc), pe1)

			cRoot.Add(cRoot, nag.NewPolynomial(k, c.Order(), nag.PolynomialTerm[K]{Coefficient: rc, Monomial: w}))
		}

		// Do factorization on c^{1/p}.
		sf := squareFree(cRoot)
		for i := range sf {
			sf[i].N.Mul(sf[i].N, p)
		}
		r = append(r, sf...)
	}
	return r
}

func differentiate[K nag.Field[K]](a *nag.Polynomial[K]) *nag.Polynomial[K] {
	k := a.Field()
	aP := poly0(a)
	for ac, aw := range a.Terms() {
		deg := int64(len(aw))
		if deg == 0 {
			continue
		}
		w := make([]nag.Symbol, deg-1)
		copy(w, aw)

		c := mulInt(ac.NewZero().Set(ac), big.NewInt(deg))
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

// gcd returns the greatest common divisor of a and b.
// gcd(a, b) = g = u*a + v*b
// a = a1*g
// b = b1*g
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

func ithPoly[K Finite[K]](poly *nag.Polynomial[K], ith, base *big.Int) *nag.Polynomial[K] {
	k := poly.Field()
	ith = new(big.Int).Set(ith)
	xPow := -1
	r := new(big.Int)
	for ith.Sign() != 0 {
		xPow++
		ith.QuoRem(ith, base, r)

		c := setIth(k.NewZero(), r)
		w := make([]nag.Symbol, xPow)
		poly.Add(poly, nag.NewPolynomial(k, poly.Order(), nag.PolynomialTerm[K]{Coefficient: c, Monomial: w}))
	}
	return poly
}

func deg[K nag.Field[K]](a *nag.Polynomial[K]) *big.Int {
	w := a.LeadingTerm().Monomial
	return big.NewInt(int64(len(w)))
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

func randPoly[K Finite[K]](poly *nag.Polynomial[K], n *big.Int) *nag.Polynomial[K] {
	k := poly.Field()
	// q is the order of the field k.
	q := new(big.Int).Exp(k.Characteristic(), big.NewInt(int64(k.Degree())), nil)

	for i := big.NewInt(0); i.Cmp(n) < 0; i.Add(i, big.NewInt(1)) {
		cint, _ := rand.Int(rand.Reader, q)
		if cint.Sign() == 0 {
			continue
		}

		c := setIth(k.NewZero(), cint)
		w := make(nag.Monomial, int(i.Uint64()))
		poly.Add(poly, nag.NewPolynomial(k, poly.Order(), nag.PolynomialTerm[K]{Coefficient: c, Monomial: w}))
	}

	return poly
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

func mulInt[K nag.Field[K]](x K, n *big.Int) K {
	return nag.Pow(&fieldAddGroup[K]{x}, n).field
}

// A polynomialMod is an element of the quotient ring K[x]/(mod), used to efficiently compute p^n mod mod
// for astronomically large n via [pow].
type polynomialMod[K nag.Field[K]] struct {
	p   *nag.Polynomial[K]
	mod *nag.Polynomial[K]
}

func newPolynomialMod[K nag.Field[K]](p, mod *nag.Polynomial[K]) *polynomialMod[K] {
	z := &polynomialMod[K]{p: nag.NewPolynomial(p.Field(), p.Order()).Set(p), mod: nag.NewPolynomial(mod.Field(), mod.Order()).Set(mod)}
	_, z.p = nag.Divide[K](nil, z.p, []*nag.Polynomial[K]{z.mod})
	return z
}

func (x *polynomialMod[K]) NewOne() *polynomialMod[K] {
	z := &polynomialMod[K]{p: nag.NewPolynomial(x.p.Field(), x.p.Order(), nag.PolynomialTerm[K]{Coefficient: x.p.Field().NewOne()}), mod: nag.NewPolynomial(x.mod.Field(), x.mod.Order()).Set(x.mod)}
	z.p.SymbolStringer = x.p.SymbolStringer
	return z
}

func (z *polynomialMod[K]) Equal(x *polynomialMod[K]) bool {
	return z.p.Equal(x.p) && z.mod.Equal(x.mod)
}

func (x *polynomialMod[K]) Set(y *polynomialMod[K]) *polynomialMod[K] {
	x.p.Set(y.p)
	x.mod.Set(y.mod)
	return x
}

func (z *polynomialMod[K]) Mul(x, y *polynomialMod[K]) *polynomialMod[K] {
	z.p.Mul(x.p, y.p)
	z.mod.Set(x.mod)
	_, z.p = nag.Divide[K](nil, z.p, []*nag.Polynomial[K]{z.mod})
	return z
}

func (z *polynomialMod[K]) Inv(x *polynomialMod[K]) *polynomialMod[K] {
	g, u, _, _, _ := gcd(x.p, x.mod)
	if !isOne(g) {
		panic("inverse does not exist")
	}
	z.mod.Set(x.mod)
	z.p.Set(u)
	return z
}

func (x *polynomialMod[K]) String() string {
	return fmt.Sprintf("%s/%s", x.p, x.mod)
}

type fieldAddGroup[K nag.Field[K]] struct {
	field K
}

func (x *fieldAddGroup[K]) NewOne() *fieldAddGroup[K] {
	return &fieldAddGroup[K]{field: x.field.NewZero()}
}

func (z *fieldAddGroup[K]) Equal(x *fieldAddGroup[K]) bool {
	return z.field.Equal(x.field)
}

func (x *fieldAddGroup[K]) Set(y *fieldAddGroup[K]) *fieldAddGroup[K] {
	x.field.Set(y.field)
	return x
}

func (z *fieldAddGroup[K]) Mul(x, y *fieldAddGroup[K]) *fieldAddGroup[K] {
	return &fieldAddGroup[K]{z.field.Add(x.field, y.field)}
}

func (z *fieldAddGroup[K]) Inv(x *fieldAddGroup[K]) *fieldAddGroup[K] {
	return &fieldAddGroup[K]{z.field.Sub(z.field.NewZero(), x.field)}
}

func (x *fieldAddGroup[K]) String() string {
	return x.field.String()
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

// Parse parses input and returns the polynomial it represents.
func Parse[K Finite[K]](variables map[string]nag.Symbol, field K, input string) (*nag.Polynomial[K], error) {
	rp, err := nag.Parse(variables, nag.Deglex, input)
	if err != nil {
		return nil, errors.Wrap(err, "")
	}

	// Cast coefficients from rationals to GF(order).
	p := nag.NewPolynomial[K](field, rp.Order())
	p.SymbolStringer = rp.SymbolStringer
	for rc, w := range rp.Terms() {
		c := field.NewZero().Div(setIth(field.NewZero(), rc.Num()), setIth(field.NewZero(), rc.Denom()))
		p.Add(p, nag.NewPolynomial(p.Field(), p.Order(), nag.PolynomialTerm[K]{Coefficient: c, Monomial: w}))
	}
	return p, nil
}
