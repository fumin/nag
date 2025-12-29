package nag

import (
	"fmt"
	"testing"
)

func TestParse(t *testing.T) {
	tests := []struct {
		variables map[string]Symbol
		order     Order
		input     string
		p         *Polynomial[*Rat]
	}{
		{
			variables: map[string]Symbol{"a": 1, "b": 2},
			order:     Deglex,
			input:     "ba^3",
			p: NewPolynomial(
				NewRat(0, 1), Deglex,
				PolynomialTerm[*Rat]{Coefficient: NewRat(1, 1), Monomial: Monomial{2, 1, 1, 1}},
			),
		},
		{
			variables: map[string]Symbol{"a": 1, `{b^{\dagger}}`: 2},
			order:     Deglex,
			input:     `-{b^{\dagger}}^1a^3`,
			p: NewPolynomial(
				NewRat(0, 1), Deglex,
				PolynomialTerm[*Rat]{Coefficient: NewRat(-1, 1), Monomial: Monomial{2, 1, 1, 1}},
			),
		},
		{
			variables: map[string]Symbol{"a": 1, "b": 2},
			order:     Deglex,
			input:     "(-ba)^3a-a(3b-5b^2)",
			p: NewPolynomial(
				NewRat(0, 1), Deglex,
				PolynomialTerm[*Rat]{Coefficient: NewRat(-1, 1), Monomial: Monomial{2, 1, 2, 1, 2, 1, 1}},
				PolynomialTerm[*Rat]{Coefficient: NewRat(-3, 1), Monomial: Monomial{1, 2}},
				PolynomialTerm[*Rat]{Coefficient: NewRat(5, 1), Monomial: Monomial{1, 2, 2}},
			),
		},
		{
			variables: map[string]Symbol{"a": 1, "b": 2},
			order:     Deglex,
			input:     "(a-b)^3",
			p: NewPolynomial(
				NewRat(0, 1), Deglex,
				PolynomialTerm[*Rat]{Coefficient: NewRat(1, 1), Monomial: Monomial{1, 1, 1}},
				PolynomialTerm[*Rat]{Coefficient: NewRat(-1, 1), Monomial: Monomial{1, 1, 2}},
				PolynomialTerm[*Rat]{Coefficient: NewRat(-1, 1), Monomial: Monomial{1, 2, 1}},
				PolynomialTerm[*Rat]{Coefficient: NewRat(1, 1), Monomial: Monomial{1, 2, 2}},
				PolynomialTerm[*Rat]{Coefficient: NewRat(-1, 1), Monomial: Monomial{2, 1, 1}},
				PolynomialTerm[*Rat]{Coefficient: NewRat(1, 1), Monomial: Monomial{2, 1, 2}},
				PolynomialTerm[*Rat]{Coefficient: NewRat(1, 1), Monomial: Monomial{2, 2, 1}},
				PolynomialTerm[*Rat]{Coefficient: NewRat(-1, 1), Monomial: Monomial{2, 2, 2}},
			),
		},
		{
			variables: map[string]Symbol{"a": 1, "b": 2, "c": 3},
			order:     Deglex,
			input:     "-12/5a^3((a+cc)b)^2a+7/3ca-3/2b",
			p: NewPolynomial(
				NewRat(0, 1), Deglex,
				PolynomialTerm[*Rat]{Coefficient: NewRat(-12, 5), Monomial: Monomial{1, 1, 1, 3, 3, 2, 3, 3, 2, 1}},
				PolynomialTerm[*Rat]{Coefficient: NewRat(-12, 5), Monomial: Monomial{1, 1, 1, 1, 2, 3, 3, 2, 1}},
				PolynomialTerm[*Rat]{Coefficient: NewRat(-12, 5), Monomial: Monomial{1, 1, 1, 3, 3, 2, 1, 2, 1}},
				PolynomialTerm[*Rat]{Coefficient: NewRat(-12, 5), Monomial: Monomial{1, 1, 1, 1, 2, 1, 2, 1}},
				PolynomialTerm[*Rat]{Coefficient: NewRat(7, 3), Monomial: Monomial{3, 1}},
				PolynomialTerm[*Rat]{Coefficient: NewRat(-3, 2), Monomial: Monomial{2}},
			),
		},
		{
			variables: map[string]Symbol{"a": 1, "b": 2, "c": 3},
			order:     Deglex,
			input:     "5/3b(a+b)^2c+9a",
			p: NewPolynomial(
				NewRat(0, 1), Deglex,
				PolynomialTerm[*Rat]{Coefficient: NewRat(5, 3), Monomial: Monomial{2, 1, 1, 3}},
				PolynomialTerm[*Rat]{Coefficient: NewRat(5, 3), Monomial: Monomial{2, 1, 2, 3}},
				PolynomialTerm[*Rat]{Coefficient: NewRat(5, 3), Monomial: Monomial{2, 2, 1, 3}},
				PolynomialTerm[*Rat]{Coefficient: NewRat(5, 3), Monomial: Monomial{2, 2, 2, 3}},
				PolynomialTerm[*Rat]{Coefficient: NewRat(9, 1), Monomial: Monomial{1}},
			),
		},
	}

	for i, test := range tests {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			t.Parallel()
			p, err := Parse(test.variables, test.order, test.input)
			if err != nil {
				t.Fatalf("%+v", err)
			}
			if !p.Equal(test.p) {
				t.Errorf("%v %v", p, test.p)
			}
		})
	}
}
