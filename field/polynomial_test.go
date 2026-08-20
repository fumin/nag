package field

import (
	"fmt"
	"math/big"
	"slices"
	"testing"

	"github.com/fumin/nag"
)

func TestPrimitive(t *testing.T) {
	tests := []struct {
		order *big.Int
		f     string
		ok    bool
	}{
		{order: big.NewInt(2), f: "x^6+x+1", ok: true},
		{order: big.NewInt(2), f: "x^6+x^3+1", ok: false},
		{order: big.NewInt(2), f: "x^6+x^4+x^2+x+1", ok: false},
		{order: big.NewInt(2), f: "x^6+x^4+x^3+x+1", ok: true},
	}

	for testI, test := range tests {
		t.Run(fmt.Sprintf("%d", testI), func(t *testing.T) {
			t.Parallel()
			f := parseMust(primeField(test.order), test.f)
			ok := primitive(f)
			if ok != test.ok {
				t.Errorf("primitive(%v): got %v want %v", f, ok, test.ok)
			}
		})
	}
}

func TestFactor(t *testing.T) {
	tests := []struct {
		order  *big.Int
		a      string
		factor []factorStr
	}{
		{
			order: big.NewInt(2),
			a:     "x^2+x",
			factor: []factorStr{
				{a: "x", n: big.NewInt(1)},
				{a: "x+1", n: big.NewInt(1)},
			},
		},
		{
			order: big.NewInt(2),
			a:     "x^2+1",
			factor: []factorStr{
				{a: "x+1", n: big.NewInt(2)},
			},
		},
		{
			order: big.NewInt(2),
			a:     "x^2+x+1",
			factor: []factorStr{
				{a: "x^2+x+1", n: big.NewInt(1)},
			},
		},
		{
			order: big.NewInt(7),
			a:     "(x+2)^6(x+5)^6(x^5+x+4)^6(x^2+x+3)^4(x^2+2x+5)^4",
			factor: []factorStr{
				{a: "x+2", n: big.NewInt(6)},
				{a: "x+5", n: big.NewInt(6)},
				{a: "x^2+x+3", n: big.NewInt(4)},
				{a: "x^2+2x+5", n: big.NewInt(4)},
				{a: "x^5+x+4", n: big.NewInt(6)},
			},
		},
		{
			order: int10("21888242871839275222246405745257275088548364400416034343698204186575808495617"),
			a:     "x^2+2x-8",
			factor: []factorStr{
				{a: "x-2", n: big.NewInt(1)},
				{a: "x+4", n: big.NewInt(1)},
			},
		},
		// The below tests are from test_gf_factor in
		// https://github.com/sympy/sympy/blob/16fa855354eb7bcabd3fe10993841e03b1382692/sympy/polys/tests/test_galoistools.py#L630
		{
			order:  big.NewInt(11),
			a:      "1",
			factor: []factorStr{},
		},
		{
			order: big.NewInt(11),
			a:     "x+1",
			factor: []factorStr{
				{a: "x+1", n: big.NewInt(1)},
			},
		},
		{
			order: big.NewInt(2),
			a:     "x^4+x",
			factor: []factorStr{
				{a: "x", n: big.NewInt(1)},
				{a: "x+1", n: big.NewInt(1)},
				{a: "x^2+x+1", n: big.NewInt(1)},
			},
		},
		{
			order: big.NewInt(11),
			a:     "x^6-3x^5+x^4-3x^3-x^2-3x+1",
			factor: []factorStr{
				{a: "x+1", n: big.NewInt(1)},
				{a: "x^2+5x+3", n: big.NewInt(1)},
				{a: "x^3+2x^2+3x+4", n: big.NewInt(1)},
			},
		},
		{
			order: big.NewInt(11),
			a:     "x^3+5x^2+8x+4",
			factor: []factorStr{
				{a: "x+1", n: big.NewInt(1)},
				{a: "x+2", n: big.NewInt(2)},
			},
		},
		{
			order: big.NewInt(11),
			a:     "x^9+x^8+10x^7+x^6+10x^4+10x^3+10x^2",
			factor: []factorStr{
				{a: "x", n: big.NewInt(2)},
				{a: "x^2+9x+5", n: big.NewInt(1)},
				{a: "x^5+3x^4+8x^2+5x+2", n: big.NewInt(1)},
			},
		},
		{
			order: big.NewInt(11),
			a:     "x^32+1",
			factor: []factorStr{
				{a: "x^16+3x^8+10", n: big.NewInt(1)},
				{a: "x^16+8x^8+10", n: big.NewInt(1)},
			},
		},
		{
			order: big.NewInt(11),
			a:     "(8x^32+5)(1/8)",
			factor: []factorStr{
				{a: "x+3", n: big.NewInt(1)},
				{a: "x+8", n: big.NewInt(1)},
				{a: "x^2+9", n: big.NewInt(1)},
				{a: "x^2+2x+2", n: big.NewInt(1)},
				{a: "x^2+9x+2", n: big.NewInt(1)},
				{a: "x^4+5x^2+7", n: big.NewInt(1)},
				{a: "x^4+6x^2+7", n: big.NewInt(1)},
				{a: "x^8+x^4+6", n: big.NewInt(1)},
				{a: "x^8+10x^4+6", n: big.NewInt(1)},
			},
		},
		{
			order: big.NewInt(11),
			a:     "(8x^63+5)(1/8)",
			factor: []factorStr{
				{a: "x+7", n: big.NewInt(1)},
				{a: "x^2+4x+5", n: big.NewInt(1)},
				{a: "x^3+6x^2+8x+2", n: big.NewInt(1)},
				{a: "x^3+9x^2+9x+2", n: big.NewInt(1)},
				{a: "x^6+9x^3+4", n: big.NewInt(1)},
				{a: "x^6+2x^5+8x^3+4x^2+6x+4", n: big.NewInt(1)},
				{a: "x^6+2x^5+6x^4+8x^2+4x+4", n: big.NewInt(1)},
				{a: "x^6+5x^5+6x^4+8x^2+6x+4", n: big.NewInt(1)},
				{a: "x^6+2x^5+3x^4+8x^3+6x+4", n: big.NewInt(1)},
				{a: "x^6+10x^5+10x^4+x^3+4x^2+9x+4", n: big.NewInt(1)},
				{a: "x^6+10x^5+4x^4+7x^3+10x^2+7x+4", n: big.NewInt(1)},
				{a: "x^6+3x^5+3x^4+x^3+6x^2+8x+4", n: big.NewInt(1)},
				{a: "x^6+6x^5+2x^4+7x^3+9x^2+8x+4", n: big.NewInt(1)},
			},
		},
		{
			// Gathen polynomial: x^n+x+1 (mod p > 2^n * pi).
			order: big.NewInt(102953),
			a:     "x^15+x+1",
			factor: []factorStr{
				{a: "x^2+22730x+68144", n: big.NewInt(1)},
				{a: "x^4+81553x^3+77449x^2+86810x+4724", n: big.NewInt(1)},
				{a: "x^4+86276x^3+56779x^2+14859x+31575", n: big.NewInt(1)},
				{a: "x^5+15347x^4+95022x^3+84569x^2+94508x+92335", n: big.NewInt(1)},
			},
		},
		{
			// Shoup polynomial: a_0 x^n + a_1 x^(n-1) + ... + a_n (mod p > 2^(n-2) * pi),
			// where a_n = a_(n-1)^2+1, a_0 = 1.
			order: big.NewInt(53),
			a:     "x^6+2x^5+5x^4+26x^3+41x^2+39x+38",
			factor: []factorStr{
				{a: "x^2+44x+26", n: big.NewInt(1)},
				{a: "x^4+11x^3+25x^2+18x+30", n: big.NewInt(1)},
			},
		},
	}

	for testI, test := range tests {
		t.Run(fmt.Sprintf("%d", testI), func(t *testing.T) {
			t.Parallel()
			a := parseMust(primeField(test.order), test.a)
			tFactors := make([]IrrFactor[*prime], len(test.factor))
			for i, f := range test.factor {
				tFactors[i] = IrrFactor[*prime]{P: parseMust(primeField(test.order), f.a), N: f.n}
			}

			factors := Factor(a)
			slices.SortFunc(factors, func(a, b IrrFactor[*prime]) int { return polynomialCmp(a.P, b.P) })
			if len(factors) != len(tFactors) {
				t.Fatalf("%v", factors)
			}
			for i, f := range factors {
				tf := tFactors[i]
				if !(f.P.Equal(tf.P) && f.N.Cmp(tf.N) == 0) {
					t.Errorf("%d got %v want %v", i, f, tf)
				}
			}
			// Check that the reconstruction from factors equals a.
			af := nag.NewPolynomial(a.Field(), a.Order(), nag.PolynomialTerm[*prime]{Coefficient: a.Field().NewOne()})
			buf := nag.NewPolynomial(a.Field(), a.Order())
			for _, f := range factors {
				for range f.N.Int64() {
					af.Set(buf.Mul(af, f.P))
				}
			}
			if !af.Equal(a) {
				t.Errorf("got %v want %v", af, a)
			}
		})
	}
}

func TestEqualDegree(t *testing.T) {
	tests := []struct {
		order  *big.Int
		a      factorStr
		factor []string
	}{
		{
			order:  big.NewInt(2),
			a:      factorStr{a: "x(x+1)", n: big.NewInt(1)},
			factor: []string{"x", "x+1"},
		},
		{
			order:  big.NewInt(3),
			a:      factorStr{a: "x(x+2)", n: big.NewInt(1)},
			factor: []string{"x", "x+2"},
		},
		{
			order:  big.NewInt(2),
			a:      factorStr{a: "(x^5+x^2+1)(x^5+x^3+1)(x^5+x^4+x^2+x+1)", n: big.NewInt(5)},
			factor: []string{"x^5+x^2+1", "x^5+x^3+1", "x^5+x^4+x^2+x+1"},
		},
		// The below test is from test_gf_edf in
		// https://github.com/sympy/sympy/blob/16fa855354eb7bcabd3fe10993841e03b1382692/sympy/polys/tests/test_galoistools.py#L615
		{
			order:  big.NewInt(3),
			a:      factorStr{a: "x^4+x^3+x+2", n: big.NewInt(2)},
			factor: []string{"x^2+1", "x^2+x+2"},
		},
		// The below test is from test_issue_23174 in
		// https://github.com/sympy/sympy/blob/16fa855354eb7bcabd3fe10993841e03b1382692/sympy/polys/tests/test_galoistools.py#L623
		{
			order:  big.NewInt(2),
			a:      factorStr{a: "x^16+x^15+x^14+x^13+x^12+x^11+x^10+x^9+x^8+x^7+x^6+x^5+x^4+x^3+x^2+x+1", n: big.NewInt(8)},
			factor: []string{"x^8+x^5+x^4+x^3+1", "x^8+x^7+x^6+x^4+x^2+x+1"},
		},
	}

	for testI, test := range tests {
		t.Run(fmt.Sprintf("%d", testI), func(t *testing.T) {
			t.Parallel()
			a := IrrFactor[*prime]{P: parseMust(primeField(test.order), test.a.a), N: test.a.n}
			tfactors := make([]*nag.Polynomial[*prime], len(test.factor))
			for i, f := range test.factor {
				tfactors[i] = parseMust(primeField(test.order), f)
			}

			factors := equalDegree(a.P, a.N)
			slices.SortFunc(factors, polynomialCmp)
			if len(factors) != len(tfactors) {
				t.Fatalf("%v", factors)
			}
			for i, f := range factors {
				if !f.Equal(tfactors[i]) {
					t.Errorf("got %v want %v", f, tfactors[i])
				}
			}
		})
	}
}

func TestDistinctDegree(t *testing.T) {
	tests := []struct {
		order  *big.Int
		a      string
		factor []factorStr
	}{
		{
			order: big.NewInt(3),
			a:     "x(x+2)(x^2+x+2)",
			factor: []factorStr{
				{a: "x(x+2)", n: big.NewInt(1)},
				{a: "x^2+x+2", n: big.NewInt(2)},
			},
		},
		// The below tests are from test_gf_ddf in
		// https://github.com/sympy/sympy/blob/16fa855354eb7bcabd3fe10993841e03b1382692/sympy/polys/tests/test_galoistools.py#L572
		{
			order: big.NewInt(11),
			a:     "x^15-1",
			factor: []factorStr{
				{a: "x^5+10", n: big.NewInt(1)},
				{a: "x^10+x^5+1", n: big.NewInt(2)},
			},
		},
		{
			order: big.NewInt(2),
			a:     "x^63+1",
			factor: []factorStr{
				{a: "x+1", n: big.NewInt(1)},
				{a: "x^2+x+1", n: big.NewInt(2)},
				{a: "x^6+x^5+x^4+x^3+x^2+x+1", n: big.NewInt(3)},
				{a: "x^54+x^53+x^51+x^50+x^48+x^46+x^45+x^43+x^42+x^33+x^32+x^30+x^29+x^27+x^25+x^24+x^22+x^21+x^12+x^11+x^9+x^8+x^6+x^4+x^3+x+1", n: big.NewInt(6)},
			},
		},
		{
			order: big.NewInt(3),
			a:     "x^6-x^5+x^4+x^3-x",
			factor: []factorStr{
				{a: "x^2+x", n: big.NewInt(1)},
				{a: "x^4+x^3+x+2", n: big.NewInt(2)},
			},
		},
		{
			order: big.NewInt(809),
			a:     "x^10+2x^9+5x^8+26x^7+677x^6+436x^5+791x^4+325x^3+456x^2+24x+577",
			factor: []factorStr{
				{a: "x+701", n: big.NewInt(1)},
				{a: "x^9+110x^8+559x^7+532x^6+694x^5+151x^4+110x^3+70x^2+735x+122", n: big.NewInt(9)},
			},
		},
		{
			order: big.NewInt(102953),
			a:     "x^15+x+1",
			factor: []factorStr{
				{a: "x^2+22730x+68144", n: big.NewInt(2)},
				{a: "x^8+64876x^7+83977x^6+10787x^5+12561x^4+68608x^3+52650x^2+88001x+84356", n: big.NewInt(4)},
				{a: "x^5+15347x^4+95022x^3+84569x^2+94508x+92335", n: big.NewInt(5)},
			},
		},
	}

	for testI, test := range tests {
		t.Run(fmt.Sprintf("%d", testI), func(t *testing.T) {
			t.Parallel()
			a := parseMust(primeField(test.order), test.a)
			tfactors := make([]IrrFactor[*prime], len(test.factor))
			for i, f := range test.factor {
				tfactors[i] = IrrFactor[*prime]{P: parseMust(primeField(test.order), f.a), N: f.n}
			}

			factors := distinctDegree(a)
			if len(factors) != len(tfactors) {
				t.Fatalf("%v", factors)
			}
			for i, f := range factors {
				tf := tfactors[i]
				if !(f.P.Equal(tf.P) && f.N.Cmp(tf.N) == 0) {
					t.Errorf("got %v want %v", f, tf)
				}
			}
		})
	}
}

func TestSquareFree(t *testing.T) {
	tests := []struct {
		order  *big.Int
		a      string
		factor []factorStr
	}{
		{
			order: big.NewInt(3),
			a:     "x^11+2x^9+2x^8+x^6+x^5+2x^3+2x^2+1",
			factor: []factorStr{
				{a: "x+1", n: big.NewInt(1)},
				{a: "x+2", n: big.NewInt(4)},
				{a: "x^2+1", n: big.NewInt(3)},
			},
		},
		{
			order:  big.NewInt(5),
			a:      "x^6+x^4+x^3-x^2-2x-1",
			factor: []factorStr{{a: "x^3+3x+3", n: big.NewInt(2)}},
		},
		{
			order: big.NewInt(13),
			a:     "x^7+3x^6+5x^5+7x^4+7x^3+5x^2+3x+1",
			factor: []factorStr{
				{a: "x^2+1", n: big.NewInt(2)},
				{a: "x+1", n: big.NewInt(3)},
			},
		},
		{
			order: big.NewInt(7),
			a:     "(x+2)^6(x+5)^6",
			factor: []factorStr{
				{a: "x^2+3", n: big.NewInt(6)},
			},
		},
	}

	for testI, test := range tests {
		t.Run(fmt.Sprintf("%d", testI), func(t *testing.T) {
			t.Parallel()
			a := parseMust(primeField(test.order), test.a)
			testFactors := make([]IrrFactor[*prime], len(test.factor))
			for i, f := range test.factor {
				testFactors[i] = IrrFactor[*prime]{P: parseMust(primeField(test.order), f.a), N: f.n}
			}

			factors := squareFree(a)
			if len(factors) != len(testFactors) {
				t.Fatalf("%v", factors)
			}
			for i, f := range factors {
				tf := testFactors[i]
				if !(f.P.Equal(tf.P) && f.N.Cmp(tf.N) == 0) {
					t.Errorf("%d got %v want %v", i, f, tf)
				}
			}
			// Check that the reconstruction from factors equals a.
			af := nag.NewPolynomial(a.Field(), a.Order(), nag.PolynomialTerm[*prime]{Coefficient: a.Field().NewOne()})
			buf := nag.NewPolynomial(a.Field(), a.Order())
			for _, f := range factors {
				for range f.N.Int64() {
					af.Set(buf.Mul(af, f.P))
				}
			}
			if !af.Equal(a) {
				t.Errorf("got %v want %v", af, a)
			}
		})
	}
}

func TestDifferentiate(t *testing.T) {
	tests := []struct {
		a  *nag.Polynomial[*prime]
		aP *nag.Polynomial[*prime]
	}{
		{
			a:  parseMust(primeField(big.NewInt(13)), "3x^7+8x^5-2x^2+5"),
			aP: parseMust(primeField(big.NewInt(13)), "-5x^6+x^4-4x"),
		},
	}

	for testI, test := range tests {
		t.Run(fmt.Sprintf("%d", testI), func(t *testing.T) {
			t.Parallel()
			aP := differentiate(test.a)
			if !aP.Equal(test.aP) {
				t.Errorf("got %v want %v", aP, test.aP)
			}
		})
	}
}

func TestInverse(t *testing.T) {
	tests := []struct {
		order *big.Int
		a     string
		p     string
		inv   string
	}{
		{
			order: big.NewInt(2),
			a:     "x^6+x^4+x+1",
			p:     "x^8+x^4+x^3+x+1",
			inv:   "x^7+x^6+x^3+x",
		},
	}

	for testI, test := range tests {
		t.Run(fmt.Sprintf("%d", testI), func(t *testing.T) {
			t.Parallel()
			a := parseMust(primeField(test.order), test.a)
			p := parseMust(primeField(test.order), test.p)
			tinv := parseMust(primeField(test.order), test.inv)
			inv := inverse(a, p)
			if !inv.Equal(tinv) {
				t.Errorf("got %v want %v", inv, tinv)
			}
		})
	}
}

func TestGcd(t *testing.T) {
	tests := []struct {
		order *big.Int
		a     string
		b     string
		g     string
		u     string
		v     string
		a1    string
		b1    string
	}{
		{
			order: big.NewInt(101),
			a:     "x^2+7x+6",
			b:     "x^2-5x-6",
			g:     "x+1",
			u:     "59",
			v:     "42",
			a1:    "x+6",
			b1:    "x-6",
		},
	}

	for testI, test := range tests {
		t.Run(fmt.Sprintf("%d", testI), func(t *testing.T) {
			t.Parallel()
			a := parseMust(primeField(test.order), test.a)
			b := parseMust(primeField(test.order), test.b)
			tg := parseMust(primeField(test.order), test.g)
			tu := parseMust(primeField(test.order), test.u)
			tv := parseMust(primeField(test.order), test.v)
			ta1 := parseMust(primeField(test.order), test.a1)
			tb1 := parseMust(primeField(test.order), test.b1)
			g, u, v, a1, b1 := gcd(a, b)
			if !g.Equal(tg) {
				t.Errorf("got %v want %v", g, tg)
			}
			if !u.Equal(tu) {
				t.Errorf("got %v want %v", u, tu)
			}
			if !v.Equal(tv) {
				t.Errorf("got %v want %v", v, tv)
			}
			if !a1.Equal(ta1) {
				t.Errorf("got %v want %v", a1, ta1)
			}
			if !b1.Equal(tb1) {
				t.Errorf("got %v want %v", b1, tb1)
			}
			// Check Bezout's identity.
			aubv := poly0(g).Add(mul(a, u), mul(b, v))
			if !aubv.Equal(tg) {
				t.Errorf("got %v want %v", aubv, tg)
			}
			if ga1 := mul(g, a1); !ga1.Equal(a) {
				t.Errorf("got %v want %v", ga1, a)
			}
			if gb1 := mul(g, b1); !gb1.Equal(b) {
				t.Errorf("got %v want %v", gb1, b)
			}
		})
	}
}

func TestIthPoly(t *testing.T) {
	tests := []struct {
		ith  *big.Int
		base *big.Int
		poly string
	}{
		{
			ith:  big.NewInt(3),
			base: big.NewInt(2),
			poly: "x+1",
		},
	}

	for testI, test := range tests {
		t.Run(fmt.Sprintf("%d", testI), func(t *testing.T) {
			t.Parallel()
			tpoly := parseMust(primeField(test.base), test.poly)
			poly := ithPoly(poly0(tpoly), test.ith, test.base)
			if !poly.Equal(tpoly) {
				t.Errorf("%v %v", poly, tpoly)
			}
		})
	}
}

func TestNewPolynomialMod(t *testing.T) {
	tests := []struct {
		order int64
		mod   string
		p     string
		e     string
	}{
		{order: 2, mod: "x^2+x+1", p: "x^4+x^3+1", e: "x"},
		{order: 3, mod: "x^3+2x+1", p: "(x^2+x+2)(2x^2+1)", e: "x^2+x"},
	}
	for i, test := range tests {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			mod := parseMust(primeField(big.NewInt(test.order)), test.mod)
			p := parseMust(primeField(big.NewInt(test.order)), test.p)
			e := parseMust(primeField(big.NewInt(test.order)), test.e)

			z := newPolynomialMod(p, mod)
			if !z.mod.Equal(mod) {
				t.Errorf("%v", z.mod)
			}
			if !z.p.Equal(e) {
				t.Errorf("%v", z.p)
			}
		})
	}
}

func TestPolynomialModMul(t *testing.T) {
	tests := []struct {
		order int64
		mod   string
		x     string
		y     string
		mul   string
	}{
		{order: 2, mod: "x^2+x+1", x: "x", y: "0", mul: "0"},
		{order: 2, mod: "x^2+x+1", x: "x", y: "x", mul: "x+1"},
		{order: 2, mod: "x^2+x+1", x: "x", y: "x+1", mul: "1"},
		{order: 3, mod: "x^3+2x+1", x: "x^2+x+2", y: "2x^2+1", mul: "x^2+x"},
	}
	for i, test := range tests {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			mod := parseMust(primeField(big.NewInt(test.order)), test.mod)
			x := newPolynomialMod(parseMust(primeField(big.NewInt(test.order)), test.x), mod)
			y := newPolynomialMod(parseMust(primeField(big.NewInt(test.order)), test.y), mod)
			mul := newPolynomialMod(parseMust(primeField(big.NewInt(test.order)), test.mul), mod)

			if z := x.NewOne().Mul(x, y); !z.Equal(mul) {
				t.Errorf("%v", z)
			}
			// Check multiply by one is a no-op.
			if z1 := x.NewOne().Mul(mul, x.NewOne()); !z1.Equal(mul) {
				t.Errorf("%v", z1)
			}
		})
	}
}

func TestPolynomialModInv(t *testing.T) {
	tests := []struct {
		order int64
		mod   string
		x     string
		inv   string
	}{
		{order: 3, mod: "x^3+2x+1", x: "x^2+1", inv: "2x^2+x+2"},
		{order: 5, mod: "x^3+3x+2", x: "x+2", inv: "-2x^2-x+1"},
	}
	for i, test := range tests {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			mod := parseMust(primeField(big.NewInt(test.order)), test.mod)
			x := newPolynomialMod(parseMust(primeField(big.NewInt(test.order)), test.x), mod)
			inv := newPolynomialMod(parseMust(primeField(big.NewInt(test.order)), test.inv), mod)

			if z := x.NewOne().Inv(x); !z.Equal(inv) {
				t.Errorf("%v", z)
			}
		})
	}
}

func int10(s string) *big.Int {
	i, ok := new(big.Int).SetString(s, 10)
	if !ok {
		panic("SetString error")
	}
	return i
}

func primeField(order *big.Int) *prime {
	return &prime{order: new(big.Int).Set(order), i: big.NewInt(0)}
}

type factorStr struct {
	a string
	n *big.Int
}

func parseMust[K Finite[K]](field K, s string) *nag.Polynomial[K] {
	variables := map[string]nag.Symbol{"x": 0}
	p, err := Parse(variables, field, s)
	if err != nil {
		panic(err)
	}
	return p
}
