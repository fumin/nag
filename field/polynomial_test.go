package field

import (
	"fmt"
	"testing"

	"github.com/fumin/nag"
)

func TestSquareFree(t *testing.T) {
	type factorStr struct {
		a   string
		pow int
	}
	tests := []struct {
		order  int
		a      string
		factor []factorStr
	}{
		{
			order: 3,
			a:     "x^11+2x^9+2x^8+x^6+x^5+2x^3+2x^2+1",
			factor: []factorStr{
				{a: "x+1", pow: 1},
				{a: "x+2", pow: 4},
				{a: "x^2+1", pow: 3},
			},
		},
		{
			order:  5,
			a:      "x^6+x^4+x^3-x^2-2x-1",
			factor: []factorStr{{a: "x^3+3x+3", pow: 2}},
		},
		{
			order: 13,
			a:     "x^7+3x^6+5x^5+7x^4+7x^3+5x^2+3x+1",
			factor: []factorStr{
				{a: "x^2+1", pow: 2},
				{a: "x+1", pow: 3},
			},
		},
	}

	for testI, test := range tests {
		t.Run(fmt.Sprintf("%d", testI), func(t *testing.T) {
			t.Parallel()
			a := parseMust(test.order, test.a)
			testFactors := make([]factor[*Prime], len(test.factor))
			for i, f := range test.factor {
				testFactors[i] = factor[*Prime]{A: parseMust(test.order, f.a), Pow: f.pow}
			}

			factors := squareFree(a)
			if len(factors) != len(testFactors) {
				t.Fatalf("%v", factors)
			}
			for i, f := range factors {
				tf := testFactors[i]
				if !(f.A.Equal(tf.A) && f.Pow == tf.Pow) {
					t.Errorf("%d got %v want %v", i, f, tf)
				}
			}
			// Check that the reconstruction from factors equals a.
			af := nag.NewPolynomial(a.Field(), a.Order(), nag.PolynomialTerm[*Prime]{Coefficient: a.Field().NewOne()})
			buf := nag.NewPolynomial(a.Field(), a.Order())
			for _, f := range factors {
				for range f.Pow {
					af.Set(buf.Mul(af, f.A))
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
		a  *nag.Polynomial[*Prime]
		aP *nag.Polynomial[*Prime]
	}{
		{
			a:  parseMust(13, "3x^7+8x^5-2x^2+5"),
			aP: parseMust(13, "-5x^6+x^4-4x"),
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

func TestGcd(t *testing.T) {
	tests := []struct {
		order int
		a     string
		b     string
		gcd   string
	}{
		{
			order: 101,
			a:     "x^2+7x+6",
			b:     "x^2-5x-6",
			gcd:   "x+1",
		},
		{
			order: 5,
			a:     "x^3-2x^2+x-2",
			b:     "x^4-2x^3+2x+1",
			gcd:   "x^2+x-1",
		},
		{
			order: 7,
			a:     "x^5+2x^4+2x^3+3x^2+6x+6",
			b:     "x^3+5x^2+x+6",
			gcd:   "x^2+2x+2",
		},
		{
			order: 101,
			a:     "-39x^80-39x^79-4x^77-4x^76-9x^68-16x^67-7x^66+3x^64+3x^63+42x^62+4x^61-38x^60+45x^55-x^54-46x^53-8x^49-20x^48-12x^47+36x^42+36x^41-33x^38-33x^37+34x^36-16x^35-50x^34-42x^28-42x^27-30x^25-30x^24+39x^22+39x^21-27x^20-27x^19+5x^17-28x^16-33x^15+12x^13+17x^12+5x^11-48x^10-30x^9+4x^8+31x^7-39x^6+17x^5-4x^4-4x^3-2x^2-5x-3",
			b:     "24x^80+24x^79+18x^77+18x^76+6x^72+6x^71-46x^69+45x^68-10x^67+26x^64+26x^63+13x^62-18x^61+17x^60+7x^59-41x^58+x^56+x^55-22x^54+46x^53-33x^52-25x^47-25x^46-42x^45-42x^44+25x^43+44x^42+x^41+26x^40+44x^39-33x^33-33x^32-17x^31+x^30+23x^29+5x^28+40x^27-32x^26+41x^25+30x^24+18x^23+27x^20+27x^19+43x^18+38x^17-5x^16-31x^14-31x^13+35x^12+37x^11+2x^10+14x^8+14x^7+4x^4+4x^3+2x^2+5x+3",
			gcd:   "x^20+x^19+26x^17+26x^16+8x^8+8x^7-41x^4-41x^3+30x^2-26x+45",
		},
	}

	for testI, test := range tests {
		t.Run(fmt.Sprintf("%d", testI), func(t *testing.T) {
			t.Parallel()
			a := parseMust(test.order, test.a)
			b := parseMust(test.order, test.b)
			testGcd := parseMust(test.order, test.gcd)
			c := gcd(a, b)
			if !c.Equal(testGcd) {
				t.Errorf("got %v want %v", c, testGcd)
			}
		})
	}
}

func parseMust(order int, s string) *nag.Polynomial[*Prime] {
	p, err := parse(order, s)
	if err != nil {
		panic(err)
	}
	return p
}
