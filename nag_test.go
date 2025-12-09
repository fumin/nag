// Tests in this file come from the following references:
//
// Clemens Hofstadler. Certifying operator identities and ideal membership of noncommutative polynomials. Master’s thesis. Johannes Kepler University Linz, Austria, 2020.
// Green, Edward L. "Noncommutative Gröbner bases, and projective resolutions." Computational Methods for Representations of Groups and Algebras: Euroconference in Essen (Germany), April 1–5, 1977. Basel: Birkhäuser Basel, 1999. 29-60.
// Kreuzer, Martin, and Xingqiang Xiu. "Non-commutative Gebauer-Moeller criteria." arXiv preprint arXiv:1302.3805 (2013).
// Mora, Teo. "An introduction to commutative and noncommutative Gröbner bases." Theoretical Computer Science 134.1 (1994): 131-173.
// Xiu, Xingqiang. "Non-commutative Gröbner bases and applications." PhD diss., Universität Passau, 2012.
// Bergman software by Jörgen Backelin, https://servus.math.su.se/bergman
// NCAlgebra Mathematica library, version 6.0.3, https://mathweb.ucsd.edu/~ncalg
package nag

import (
	"flag"
	"fmt"
	"log"
	"math/big"
	"slices"
	"testing"
)

func TestBuchberger(t *testing.T) {
	tests := []struct {
		ideal    []*Polynomial
		maxiter  int
		basis    []*Polynomial
		complete bool
		long     bool
	}{
		// Example 5.12, Mora.
		{
			ideal: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 2, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 1, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2}},
				),
			},
			maxiter: 10,
			basis: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1, 2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1, 2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 1, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2}},
				),
			},
			complete: true,
		},
		// Example 4.1.15, Xiu Xingqiang.
		{
			ideal: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{4, 6, 5, 6, 5, 6, 5}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{4, 6, 5, 6, 5}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{4}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{3}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{5, 6, 5, 6, 5, 6, 2}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{5, 6, 5, 6, 2}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1}},
				),
			},
			maxiter: 4,
			basis: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{4, 6, 1}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{3, 6, 2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{4, 6, 5, 6, 1}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{3, 6, 5, 6, 2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{4, 6, 5, 6, 5, 6, 1}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{3, 6, 5, 6, 5, 6, 2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{4, 6, 5, 6, 5, 6, 5}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{4, 6, 5, 6, 5}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{4}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{3}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{5, 6, 5, 6, 5, 6, 2}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{5, 6, 5, 6, 2}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1}},
				),
			},
			complete: true,
		},
		// Section 6.5 Polynomials and Rules, NCAlgebra.
		{
			ideal: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 2, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2}},
				),
			},
			maxiter: 13,
			basis: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1, 2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 1, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 1, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2}},
				),
			},
			complete: true,
		},
		// Section 6.6 Polynomials and Rules, NCAlgebra.
		{
			ideal: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 1, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-2, 1), Monomial: Monomial{1, 2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 1, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-2, 1), Monomial: Monomial{2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1}},
				),
			},
			maxiter: 20,
			basis: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 1}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 1}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 2), Monomial: Monomial{1}},
				),
			},
			complete: true,
		},
		// Section 1.2.3 Non-commutative algebras, Bergman manual.
		{
			ideal: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2, 2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 2}},
				),
			},
			maxiter: 9,
			basis: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1, 1}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 1, 1}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 1, 1}},
				),
			},
			complete: true,
		},
		// Section 6.9.1 Lex Order, NCAlgebra.
		{
			ideal: []*Polynomial{
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2, 3}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{3, 4, 1}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{3, 2, 1, 1}},
				),
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{3, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{}},
				),
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 3}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{}},
				),
			},
			maxiter: 7,
			basis: []*Polynomial{
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 3}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{}},
				),
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{3, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{}},
				),
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{4}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1, 2, 3, 3}},
				),
			},
			complete: true,
		},
		// Section 2.5 Homogenisation, Bergman manual.
		{
			ideal: []*Polynomial{
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{5, 5}},
					PolynomialTerm{Coefficient: big.NewRat(-2, 1), Monomial: Monomial{}},
				),
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{4, 4}},
					PolynomialTerm{Coefficient: big.NewRat(-3, 1), Monomial: Monomial{}},
				),
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{3, 3}},
					PolynomialTerm{Coefficient: big.NewRat(-5, 1), Monomial: Monomial{}},
				),
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{5}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{4}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{3}},
				),
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{5, 4}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{4, 5}},
				),
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{5, 3}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{3, 5}},
				),
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{5, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2, 5}},
				),
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{4, 3}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{3, 4}},
				),
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{4, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2, 4}},
				),
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{3, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2, 3}},
				),
			},
			maxiter: 28,
			basis: []*Polynomial{
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 2, 2, 2, 2, 2, 2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-40, 1), Monomial: Monomial{2, 2, 2, 2, 2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(352, 1), Monomial: Monomial{2, 2, 2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-960, 1), Monomial: Monomial{2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(576, 1), Monomial: Monomial{}},
				),
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(576, 576), Monomial: Monomial{3}},
					PolynomialTerm{Coefficient: big.NewRat(-5, 576), Monomial: Monomial{2, 2, 2, 2, 2, 2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(194, 576), Monomial: Monomial{2, 2, 2, 2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1520, 576), Monomial: Monomial{2, 2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(2544, 576), Monomial: Monomial{2}},
				),
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(96, 96), Monomial: Monomial{4}},
					PolynomialTerm{Coefficient: big.NewRat(1, 96), Monomial: Monomial{2, 2, 2, 2, 2, 2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-37, 96), Monomial: Monomial{2, 2, 2, 2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(244, 96), Monomial: Monomial{2, 2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-360, 96), Monomial: Monomial{2}},
				),
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(576, 576), Monomial: Monomial{5}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 576), Monomial: Monomial{2, 2, 2, 2, 2, 2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(28, 576), Monomial: Monomial{2, 2, 2, 2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(56, 576), Monomial: Monomial{2, 2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-960, 576), Monomial: Monomial{2}},
				),
			},
			complete: true,
		},
		// G1, Example 4.2.26 Xiu Xingqiang, expected basis from bergman.
		{
			ideal: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "(ababab^2ab^2)^2-1"),
			},
			maxiter: 689,
			basis: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abab-b^2ab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababab^2-bab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2ab-b^2ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ab^2-babab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ababab-bababab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababab^2abab-babab^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababab^2ab-bab^2abababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abababab^2-b^2abababab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2abab^2abab-babab^2abab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2ab^2ab^2ab-ab^2abababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2ab^2aba-b^2ab^2abababababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2abab^2ababa-babab^2abababababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababababababab^2-bab^2ababababab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ababababab^2ab-b^2ababababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2abab^2abab^2abab-abab^2ababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2abababababab^2a-ab^2ab^2ab^2ab^2ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2abab^2abab^2ab^2-ababababababab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2ababababababa-ab^2abab^2abab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2abab^2ab^2a-b^2ababababababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ababababababab^2aba-abab^2abab^2abab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababababababab^2ab-bab^2abababababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2abababababab-babab^2abab^2abababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abababababababab^2-b^2abababababababab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2abab^2abab^2abab^2abab-ab^2abababababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2abab^2abab^2abab^2ab-abababababababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2abab^2abab^2aba-b^2abababababababababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ababababababababab-bab^2abababababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababababababababab^2ab-bababababababababab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2abababababababababa-ab^2abab^2abab^2abab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abababababababababab^2a-abab^2abab^2abab^2abab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababababababababab^2-bab^2abab^2abab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababababababababababab-babababababababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababababababababababab-bababababababababababab^2a"),
			},
			complete: true,
		},
		// G2, Example 4.2.26 Xiu Xingqiang, expected basis from bergman.
		{
			ideal: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "(ababab^2)^3-1"),
			},
			maxiter: 2395,
			basis: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2ababa-abab^2ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2aba-babab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2abab-bab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abab^2ab^2ab-abab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ababab^2aba-ab^2ab^2abab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2abab^2-ababab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ababab^2a-ab^2abab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abab^2a-b^2ababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2ab^2a-bab^2ababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababab^2ab-b^2ab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ababab^2-b^2abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2abab^2ab^2-ab^2ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababababab^2ab-bab^2ababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ababababab^2-b^2ababababab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababababab^2ab^2ab-bab^2ab^2ababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ababababab^2-b^2ababababab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2ab^2ab^2abab-babab^2ab^2ab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ab^2ab^2ab^2ab-bab^2ab^2ab^2ab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2ababab^2-b^2ababab^2ab^2ab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2abab^2abab^2a-bab^2abababababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2abab^2abab^2abab^2-ab^2abababababab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abab^2abab^2abab^2ab^2-abab^2abababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2abababababab^2aba-ab^2abab^2abab^2abab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abababab^2ab^2ab-b^2ababab^2ab^2ab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2ab^2ab^2ab^2aba-babab^2abab^2abab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abababababab^2ab-b^2abab^2abab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ab^2ab^2ababab^2-bab^2ab^2abababab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abababababababab^2aba-ab^2ab^2ababababababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ababababababab^2-abababababababab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2abababababababab^2a-ab^2ababababababab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ababababababab^2a-b^2abababababababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababababababab^2ab^2a-bab^2abababababababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababababababab^2ab-b^2ab^2ababababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abababababababab^2-b^2ababababababab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ababababababab^2ab^2-ab^2abababababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2abababababab^2ab-bab^2ab^2ab^2ab^2abababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2ab^2abababab^2ab^2a-bab^2ab^2abababababababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2abababababab^2-b^2ab^2ab^2ab^2abababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2ab^2ab^2ab^2ab^2abab-ababab^2ab^2ab^2abababab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2ab^2ab^2ab^2ab^2ab^2ab^2ab-abab^2abababab^2ab^2ab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2ab^2ab^2ab^2ab^2ab^2ababa-abab^2ab^2ab^2ab^2ab^2ab^2ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2ab^2ab^2ab^2ababa-b^2ababab^2ab^2ab^2abababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababababababababab^2ab-bab^2abababababababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2ab^2ab^2ab^2ab^2ab^2aba-babab^2ab^2ab^2ab^2ab^2ab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abababababababababab^2-b^2abababababababababab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ababab^2ab^2ab^2abababab^2aba-ab^2ab^2ab^2ab^2ab^2ab^2ab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2ab^2ab^2ababababa-bab^2ab^2ab^2ab^2ababababababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abababababababababab^2-b^2abababababababababab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ababababababababab^2aba-b^2ab^2ab^2ab^2ab^2ab^2ababababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ab^2ab^2ab^2ab^2ab^2abababab-ab^2ab^2ab^2ab^2abababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ab^2ab^2ab^2ab^2abababababa-abab^2ababababababababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2ab^2ab^2abababab^2ab^2-ababab^2ab^2ab^2abababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2ab^2abababababababa-ab^2ab^2ab^2ab^2ab^2ab^2abababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2ab^2abababab^2ab^2a-b^2ababab^2ab^2ab^2abababababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ababab^2ab^2ab^2abababababab^2a-ab^2ab^2ab^2ab^2ab^2abababab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ababababababababababab^2ab-ab^2abababab^2ab^2ab^2abababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2abababababababab^2-b^2ab^2ab^2ab^2ab^2ab^2ab^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababababababababababab^2aba-b^2ab^2abababab^2ab^2ab^2abababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ababababababababababab^2a-b^2abababab^2ab^2ab^2abababab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2abababab^2ab^2ab^2abababab^2a-ab^2ababababababababababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2ab^2ab^2ababababababab-abab^2ab^2ab^2ab^2ab^2ab^2ab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ababababababababab^2ab-bab^2abababab^2ab^2ab^2abababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababab^2ab^2ab^2abababab^2ab^2-bab^2ababababababababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abababab^2ab^2ab^2ababababab-b^2ab^2ab^2ababababababababab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ab^2ab^2abababababababab^2-bab^2ab^2ababababababababab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababababababababababababababab-bababababababababababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2ababababababababab^2-b^2abababab^2ab^2ab^2abababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababab^2ab^2ab^2abababababab^2-b^2ab^2ab^2ab^2ababababababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababababababababababababab-babababababababababababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ababababababababababababab^2aba-ab^2ab^2abababababababababababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2abababababababababababab^2-ababababababababababababab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2abababababababababababab^2aba-abababab^2ab^2ab^2abababababababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ababababababababababababab^2a-ab^2abababababababababababab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2ab^2ab^2ab^2ab^2ab^2ab^2ab-b^2ab^2ab^2ab^2ababababababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2abababababababababab-b^2ab^2ab^2ab^2ab^2ab^2ab^2ab^2ab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ababababababababababab^2-bab^2ab^2ab^2ab^2ab^2ab^2ab^2ab^2ab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abababababababababababab^2a-b^2ababababababababababababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ababababababababababababab-bababababababababababababab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababab^2ab^2ab^2ababababababab-b^2abababab^2ab^2ab^2abababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababababababababababab^2ab^2a-bab^2ababababababababababababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababababababababababababab^2ab-b^2ab^2abababababababababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2ab^2ab^2ab^2ab^2ab^2ab^2ab^2ab^2-b^2ab^2ab^2ababababababababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abababababababababababab^2ab-bababab^2ab^2ab^2ababababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ababababababababababababab^2-b^2abababababababababababab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ab^2ab^2ababababababababab^2-babab^2ab^2ab^2ababababababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababab^2ab^2ab^2ababababababababa-bab^2abababababababababababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababababababababababababab^2ab^2-b^2ab^2abababababababababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2abababababababababababab^2ab^2-ab^2ababababababababababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abababab^2ab^2ab^2abababababababab-ab^2abababababababababababab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2abababababababababababab-babab^2ab^2ab^2ababababababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ab^2ab^2abababababababababab-b^2ab^2ab^2ababababababababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abababababababababababababab^2aba-abab^2ab^2ab^2ababababababababababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2ab^2ab^2ababababababababababab-abababababababababababababab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2ab^2abababababababababababa-b^2abababababababababababababab^2ab"),
			},
			complete: true,
			long:     true,
		},
		// G3, Example 4.2.26 Xiu Xingqiang, expected basis from bergman.
		{
			ideal: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "(abab^2)^2-1"),
			},
			maxiter: 381,
			basis: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab-a^2ba^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2aba-b^2a^2b"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2a-ba^2b^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abab^2-a^2b^2a^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2a^2ba^2-ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ba^2b^2a^2-abab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababa^2b-ba^2baba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2b^2a^2b-b^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2ba^2b^2-bab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "aba^2b^2ab-bab^2a^2ba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2babab^2-b^2ababa^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2ba^2bab-baba^2ba^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2a^2bab^2-b^2aba^2b^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2aba^2bab^2-a^2bababa^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "baba^2ba^2ba-aba^2babab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ba^2b^2ab^2a^2-abababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ba^2bababa^2-aba^2bab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2a^2b^2-b^2a^2b^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2a^2ba-ba^2bababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababababab-bababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2b^2ab^2a^2b-b^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2a^2b^2ab-a^2babababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abababab^2-a^2b^2ab^2a^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2aba^2babab-aba^2ba^2ba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2a^2bababab-a^2ba^2ba^2ba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2a^2b^2ab^2-ababababa^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2a^2babab-aba^2ba^2ba^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2a^2baba^2-a^2bababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "baba^2bababa-ab^2a^2babab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ba^2babababa-ab^2a^2b^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ba^2ba^2ba^2ba-a^2bababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2a^2b^2ab^2a-b^2abababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2a^2bababa-b^2aba^2ba^2b"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "aba^2bababab-b^2abababa^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "aba^2ba^2ba^2b-babababa^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2babababa^2-ba^2ba^2ba^2b"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2ba^2ba^2ba^2-bab^2a^2b^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ababababa^2-ab^2a^2b^2ab^2"),
			},
			complete: true,
		},
		// G4, Example 4.2.26 Xiu Xingqiang, expected basis from bergman.
		{
			ideal: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "(aba^2b^2)^2-1"),
			},
			maxiter: 342,
			basis: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2a^2b-b^2aba^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "aba^2b^2-bab^2a^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2b^2ab-b^2a^2ba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2bab^2-ba^2b^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab-bab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "aba^2bab-baba^2ba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2-b^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2baba^2b-ba^2baba^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "aba^2ba^2ba-bab^2ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2aba^2-b^2aba^2ba^2b"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababa^2b-b^2abababa^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababa^2ba-bab^2ab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2bababab^2-ba^2babab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2a^2-b^2abababa^2b"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababa^2-b^2a^2ba^2ba^2b"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2aba^2ba^2b-b^2a^2bababa^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ababab-bab^2ab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababab^2a^2-b^2a^2bababa^2b"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababab^2-b^2abababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababababa-ba^2babababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababababab^2a-bababababa^2b"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababababab-babababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2babababa^2b-ba^2babababa^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2a-b^2ababababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2bababababa^2-b^2ab^2abababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ababa-b^2a^2babababab^2"),
			},
			complete: true,
		},
		// G5, Example 4.2.26 Xiu Xingqiang, expected basis from bergman.
		{
			ideal: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^5-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "(abab^2)^2-1"),
			},
			maxiter: 183,
			basis: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^5-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab-ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab-bababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abab^2-ab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab-babab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^3-ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^4-abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3abab-babab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^3ab-abab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4abab-babab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^4ab^2-abab^3aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^4ab-ab^2ab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^3ab^2-b^2ab^3ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^3aba-b^3ab^3ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^3ab^2a-b^2ab^3ab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^3ab^3-abab^4aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^3ab^2a-ab^2ab^3ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^3ab^3a-abab^3ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^3ab^2-b^2ab^3aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^3ab^3ab-bab^3ab^3aba"),
			},
			complete: true,
		},
		// G6, Example 4.2.26 Xiu Xingqiang, expected basis from bergman.
		{
			ideal: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^5-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "(ababab^4)^2-1"),
			},
			maxiter: 3095,
			basis: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^5-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^4abab-ab^4aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^4ab-abab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^3abab-babab^3aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2ab-bab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ababa-bab^2ab^3ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^3ab-bab^3ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2aba-bab^3ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^4ab-ab^4ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^3aba-ab^4ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4abab^4-abab^4aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4abab^3a-ab^2ab^4ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^4aba-ab^4abab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3abab^4a-abab^4ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^4ab^4-ababab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^4ab^3a-ab^2abab^4"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^4ab^2ab-ab^3aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ab^4a-abab^2ab^4"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ab^3ab-ababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ab^2aba-ababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^4ab-abab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^3aba-abab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4ab^3ab-b^4ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4ab^2aba-b^4ab^3ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4abab^3-b^2ab^4aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4abab^2a-b^3ab^4ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^4ab-b^4abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3abab^4-bab^4ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^4aba-b^4abab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^4a-bab^4ab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^3ab-bab^3abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^4ab^3-b^2abab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^4ab^2a-b^3abab^4"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^3ab^4-bab^2ab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^3ab^2ab-babab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^4a-bab^3ab^4"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^3ab-bab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4abababab-ab^3ab^3aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babababab^4-abab^3ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3abab^3ab-bab^3abab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^3ab-bab^3ab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4abab^2abab-ab^2ab^3aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ababab^4-ab^4ab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ababab^2ab-ab^3ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ababab^3-ab^2ab^2ab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3abababab-ab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ababab^4-abab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2abab^4-abab^3ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4ababab^3-b^3ababab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4ababab^2a-b^3ab^2abab^4"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ababab^4-b^4ababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^3abab^3-b^3abab^3ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abab^4-b^4abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababab^4a-b^4abab^2ab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababababab^2ab-bab^2ababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababababababa-bab^2ababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2abab^4-ab^4ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2abab^2ab-ab^3ababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2abababa-ab^3ababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4abab^2ab^4-ab^3ab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4abab^2ab^3a-ab^2ababab^4"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3abab^2ab^3-ab^3ab^2ab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3abab^2abab-ab^3abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ababab^2ab-ab^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2abab^2ab^4-ababababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^3abab^3-b^3abab^3ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2abab^4-b^4ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3abab^2ab^4-b^4ab^2abab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^3ab^2ab^2a-bab^2ab^2abab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2abab-babab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababab^2ab^4-bab^2ababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2ab^2ab^4-ab^3ababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2ab^2ab^2ab-ab^3ab^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2ab^2ababa-ab^3ab^2abab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^3abab^2aba-ab^3abab^2abab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3abab^2abab^3a-ab^2ab^3abab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^3abab^2ab^2a-ab^2abab^2abab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ab^2abab^3a-abab^2abab^2ab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2ab^4-abababab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2abab^2aba-ab^2ab^3ab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4ab^2ab^2ab^3-b^3ab^2ab^2ab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4ab^2ab^2abab-b^3ab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^3abab^2ab-b^3abab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ab^2ab^4-b^4ab^2ab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2abab^2aba-b^4ab^2ab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3abab^2abab^3-b^2ab^3abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3abab^2abab^2a-b^3ab^3abab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^3abab^2ab^2-b^2abab^2abab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^3abab^2aba-b^3abab^2abab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2abab^3-bab^2abab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ababab^2ab-bab^2ababab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2abab^3a-b^2ab^3abab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2abab^2ab-bab^2abab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2ab^3a-b^2ab^2ab^2abab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^3ab^2ab^4-ab^3abab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2ab^3ab^4-ab^3ab^2abab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2ab^3ab^2ab-ab^3ab^3abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ab^2ab^2ab^3a-ab^2ababab^2ab^4"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ab^2ab^2ab^2ab-ab^2abababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ab^2ab^2ababa-ab^2abababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^3ab^2ab^4-ab^2ab^3ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ababababab^3-ab^2ab^2ababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4ab^2ab^3ab^3-b^3ab^3ab^2ab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4ab^2ab^3ab^2a-b^2ab^2abab^2abab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^3ab^2ab^4-b^4ab^2ab^3ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ababababab^3-b^3ababababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^3ab^2ab^4a-b^3abab^2abab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2abab^2ab-bab^2abababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2abababa-bab^2abababababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abab^2abab^2-b^2abab^2abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^3ab^2ab^3ab-ab^3abababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^3ab^2ab^2aba-ab^3abababababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2abab^2abab^3-ab^3abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3abab^2abab^2ab^3-ab^2ab^2abab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3abab^2abab^2ab^2a-ab^2ab^3ab^2ab^4"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ab^2ab^3ab^4-ab^2abababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ab^2ab^3ab^3a-ab^3ab^2ab^2abab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ab^2ab^3ab^2ab-ab^3ab^2ab^3ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^3ab^2ab^3ab-ab^2ab^3ab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2abababab^3-ab^2ab^2abab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^3ab^2ab^3ab-b^4ab^2ab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ab^2abab^2ab-bab^3ab^2ab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2abababab^3-b^3abababab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3abababab^2ab^3-b^3ab^2abababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2ab^2ab-bab^2ab^2ab^2ab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^3ab^2ab^3ab^2a-ab^3ab^2ab^3ab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2ab^3ab^2ab^3a-ab^2ab^3ab^2ab^3ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2ab^2ab^2ababab-ab^3ab^2ab^2ab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2ababab^3-ab^2ab^2ab^2ab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2abababab^2ab-ab^2ab^2abababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2abababab^2ab^3a-ab^2abababababab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ab^3ab^2ab^3-b^2ab^3ab^2ab^3ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ab^2ab^2abab^2-b^3ab^2ab^2ababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ab^2ababab^3-b^3ab^2ab^2ab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2abababab^2ab-b^3abababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3abababababab^2a-b^3ab^2abababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^3ab^2ab^3ab^3-b^2ab^2ab^3ab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^3ab^2ab^3ab^2a-b^3ab^2ab^3ab^2ab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abababababab-bab^2ab^2abababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababababab^3a-bab^2abababab^2ab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2ab^2abab^2ab^3-ab^2ab^2abababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2ababab^2ab^4-ab^3ab^2ab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^3ab^2ab^2abab^3-ab^3abababab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2ab^2abababab^2a-ab^2ab^2abababab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2abababab^2ab^2a-ab^3ab^2ab^2abababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ab^2ab^2ab^2abab^2-ab^2ab^2ababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2ab^2abab^2a-ab^3ab^3ab^2ab^2ab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2ab^2ababab-ab^2ab^2ab^2ab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^3ab^2ab^2abab^2-b^3ab^2ab^2abab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ab^2ab^2ab^2ab^2-b^2ab^2ab^2ab^2ab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ab^2abababab^2-b^2ab^2abababab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ab^2ababababa-b^3ab^2abababab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ababab^2ab^4-b^4ab^2ababab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ababab^2ab^3a-b^3ab^2ab^2ab^2ab^2ab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abababab^2ab^3-b^2ab^2ab^2abababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abababab^2ab^2a-b^3ab^2ab^2abababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2ab^2ab^2ab^2ab^3a-ab^3ab^2ababab^2ab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2ab^2ab^2ab^2ababa-ab^2ab^2ab^2ab^2ababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2ab^2ababababab^2-ab^2ab^2ab^2ab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2abababababab^3-ab^3ab^2ab^2ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ab^2ab^2ab^2abab-b^2ab^2ab^2ab^2ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2abababababab^2-b^2ab^2ab^2ababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2ababab^2a-b^3ab^2ab^2ab^2ab^2abab"),
			},
			complete: true,
			long:     true,
		},
		// G7, Example 4.2.26 Xiu Xingqiang, expected basis from bergman.
		{
			ideal: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^5-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "(abab^2ab^4)^2-1"),
			},
			maxiter: 5323,
			basis: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^5-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^4abab-ab^4aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^4ab-abab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^3abab-babab^3aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2ab-bab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ababa-bab^2ab^3ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^3ab-bab^3ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2aba-bab^3ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^4ab-ab^4ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^3aba-ab^4ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4abab^4-abab^4aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4abab^3a-ab^2ab^4ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^4aba-ab^4abab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3abab^4a-abab^4ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^4ab^4-ababab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^4ab^3a-ab^2abab^4"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^4ab^2ab-ab^3aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ab^4a-abab^2ab^4"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ab^3ab-ababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ab^2aba-ababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^4ab-abab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^3aba-abab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4ab^3ab-b^4ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4ab^2aba-b^4ab^3ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4abab^3-b^2ab^4aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4abab^2a-b^3ab^4ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^4ab-b^4abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3abab^4-bab^4ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^4aba-b^4abab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^4a-bab^4ab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^3ab-bab^3abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^4ab^3-b^2abab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^4ab^2a-b^3abab^4"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^3ab^4-bab^2ab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^3ab^2ab-babab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^4a-bab^3ab^4"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^3ab-bab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4abababab-ab^3ab^3aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babababab^4-abab^3ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3abab^3ab-bab^3abab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^3ab-bab^3ab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4abab^2abab-ab^2ab^3aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ababab^4-ab^4ab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ababab^2ab-ab^3ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ababab^3-ab^2ab^2ab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3abababab-ab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ababab^4-abab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2abab^4-abab^3ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4ababab^3-b^3ababab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4ababab^2a-b^3ab^2abab^4"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ababab^4-b^4ababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^3abab^3-b^3abab^3ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abab^4-b^4abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababab^4a-b^4abab^2ab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababababab^2ab-bab^2ababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababababababa-bab^2ababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2abab^4-ab^4ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2abab^2ab-ab^3ababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2abababa-ab^3ababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4abab^2ab^4-ab^3ab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4abab^2ab^3a-ab^2ababab^4"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3abab^2ab^3-ab^3ab^2ab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3abab^2abab-ab^3abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ababab^2ab-ab^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2abab^2ab^4-ababababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^3abab^3-b^3abab^3ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2abab^4-b^4ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3abab^2ab^4-b^4ab^2abab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^3ab^2ab^2a-bab^2ab^2abab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2abab-babab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababab^2ab^4-bab^2ababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2ab^2ab^4-ab^3ababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2ab^2ab^2ab-ab^3ab^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2ab^2ababa-ab^3ab^2abab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^3abab^2aba-ab^3abab^2abab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3abab^2abab^3a-ab^2ab^3abab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^3abab^2ab^2a-ab^2abab^2abab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ab^2abab^3a-abab^2abab^2ab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2ab^4-abababab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2abab^2aba-ab^2ab^3ab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4ab^2ab^2ab^3-b^3ab^2ab^2ab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4ab^2ab^2abab-b^3ab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^3abab^2ab-b^3abab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ab^2ab^4-b^4ab^2ab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2abab^2aba-b^4ab^2ab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3abab^2abab^3-b^2ab^3abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3abab^2abab^2a-b^3ab^3abab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^3abab^2ab^2-b^2abab^2abab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^3abab^2aba-b^3abab^2abab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2abab^3-bab^2abab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ababab^2ab-bab^2ababab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2abab^3a-b^2ab^3abab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2abab^2ab-bab^2abab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2ab^3a-b^2ab^2ab^2abab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^3ab^2ab^4-ab^3abab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2ab^3ab^4-ab^3ab^2abab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2ab^3ab^2ab-ab^3ab^3abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ab^2ab^2ab^3a-ab^2ababab^2ab^4"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ab^2ab^2ab^2ab-ab^2abababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ab^2ab^2ababa-ab^2abababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^3ab^2ab^4-ab^2ab^3ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ababababab^3-ab^2ab^2ababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4ab^2ab^3ab^3-b^3ab^3ab^2ab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4ab^2ab^3ab^2a-b^2ab^2abab^2abab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^3ab^2ab^4-b^4ab^2ab^3ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ababababab^3-b^3ababababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^3ab^2ab^4a-b^3abab^2abab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2abab^2ab-bab^2abababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2abababa-bab^2abababababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abab^2abab^2-b^2abab^2abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^3ab^2ab^3ab-ab^3abababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^3ab^2ab^2aba-ab^3abababababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2abab^2abab^3-ab^3abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3abab^2abab^2ab^3-ab^2ab^2abab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3abab^2abab^2ab^2a-ab^2ab^3ab^2ab^4"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ab^2ab^3ab^4-ab^2abababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ab^2ab^3ab^3a-ab^3ab^2ab^2abab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ab^2ab^3ab^2ab-ab^3ab^2ab^3ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^3ab^2ab^3ab-ab^2ab^3ab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2abababab^3-ab^2ab^2abab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^3ab^2ab^3ab-b^4ab^2ab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ab^2abab^2ab-bab^3ab^2ab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2abababab^3-b^3abababab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3abababab^2ab^3-b^3ab^2abababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2ab^2ab-bab^2ab^2ab^2ab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^3ab^2ab^3ab^2a-ab^3ab^2ab^3ab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2ab^3ab^2ab^3a-ab^2ab^3ab^2ab^3ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2ab^2ab^2ababab-ab^3ab^2ab^2ab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2ababab^3-ab^2ab^2ab^2ab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2abababab^2ab-ab^2ab^2abababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2abababab^2ab^3a-ab^2abababababab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ab^3ab^2ab^3-b^2ab^3ab^2ab^3ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ab^2ab^2abab^2-b^3ab^2ab^2ababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ab^2ababab^3-b^3ab^2ab^2ab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2abababab^2ab-b^3abababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3abababababab^2a-b^3ab^2abababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^3ab^2ab^3ab^3-b^2ab^2ab^3ab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^3ab^2ab^3ab^2a-b^3ab^2ab^3ab^2ab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abababababab-bab^2ab^2abababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababababab^3a-bab^2abababab^2ab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2ab^2abab^2ab^3-ab^2ab^2abababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2ababab^2ab^4-ab^3ab^2ab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^3ab^2ab^2abab^3-ab^3abababab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2ab^2abababab^2a-ab^2ab^2abababab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2abababab^2ab^2a-ab^3ab^2ab^2abababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ab^2ab^2ab^2abab^2-ab^2ab^2ababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2ab^2abab^2a-ab^3ab^3ab^2ab^2ab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2ab^2ababab-ab^2ab^2ab^2ab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^3ab^2ab^2abab^2-b^3ab^2ab^2abab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ab^2ab^2ab^2ab^2-b^2ab^2ab^2ab^2ab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ab^2abababab^2-b^2ab^2abababab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ab^2ababababa-b^3ab^2abababab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ababab^2ab^4-b^4ab^2ababab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ababab^2ab^3a-b^3ab^2ab^2ab^2ab^2ab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abababab^2ab^3-b^2ab^2ab^2abababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abababab^2ab^2a-b^3ab^2ab^2abababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2ab^2ab^2ab^2ab^3a-ab^3ab^2ababab^2ab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2ab^2ab^2ab^2ababa-ab^2ab^2ab^2ab^2ababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2ab^2ababababab^2-ab^2ab^2ab^2ab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2abababababab^3-ab^3ab^2ab^2ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ab^2ab^2ab^2abab-b^2ab^2ab^2ab^2ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2abababababab^2-b^2ab^2ab^2ababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2ababab^2a-b^3ab^2ab^2ab^2ab^2abab"),
			},
			complete: true,
			long:     true,
		},
		// G8, Example 4.2.26 Xiu Xingqiang, expected basis from bergman.
		{
			ideal: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "(ababab^3)^2-1"),
			},
			maxiter: 309,
			basis: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3abab-ab^3aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^3ab-abab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ab-b^3ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3abab^2-b^2ab^3aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ababa-b^3ab^3ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^3ab-b^3abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^3-bab^3ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^3ab^2-b^2abab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^3aba-b^3abab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^3-bab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab-babab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^3a-bab^3ab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ab-bab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababababa-bab^2ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^3aba-ab^3abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3abab^3a-abab^3ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ab^3a-ababab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2aba-abababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^3ab-b^3ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2ab-bab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2ab-babababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ababab^3-ab^3ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3abababab-ab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babababab^3-abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^3a-bab^2ababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab-bab^2ab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abab^2-b^2abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababab^3-b^3ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3abab^2ab^2a-abababab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ababab^2a-ab^2ab^2ab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2abababab-ab^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ababab-b^2abababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababab^2-b^2ab^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2ab^2ab^3-ab^2ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2abababa-ab^2abababab"),
			},
			complete: true,
		},
		// G9, Example 4.2.26 Xiu Xingqiang, expected basis from bergman.
		{
			ideal: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "(abab^2)^2-1"),
			},
			maxiter: 30,
			basis: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab-b^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2-bab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab-bababa"),
			},
			complete: true,
		},
		// G10, Example 4.2.26 Xiu Xingqiang, expected basis from bergman.
		{
			ideal: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "(ababab^2)^2-1"),
			},
			maxiter: 115,
			basis: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2aba-ab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abab^2a-abab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2a-ababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2abab-ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2ab-abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab-b^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2-bab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababa-b^2ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2-babab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2aba-b^2abab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2a-bab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababab-babababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ababab^2-ab^2ab^2a"),
			},
			complete: true,
		},
		// G11, Example 4.2.26 Xiu Xingqiang, expected basis from bergman.
		{
			ideal: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "(abababab^2)^2-1"),
			},
			maxiter: 325,
			basis: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ababab-ab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2abab-ab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bababab^2ab-abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababa-b^2ab^2ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ababa-b^2ab^2abab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2aba-b^2abab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababab^2a-bab^2ab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababababab-bababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ab^2aba-ab^2ababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2abab^2a-abab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abab^2ab^2a-ababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abababab^2-ab^2ab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2a-abababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab-b^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abab^2-bab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2ab^2-babab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2ab^2-bababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2ab-bab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ab^2abab-babab^2ab^2ababa"),
			},
			complete: true,
		},
		// G12, Example 4.2.26 Xiu Xingqiang, expected basis from bergman.
		{
			ideal: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "(ababab^2abab^2)^2-1"),
			},
			maxiter: 1256,
			basis: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2abab^2aba-ab^2abab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abab^2ab^2aba-ab^2ababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abab^2abab^2a-abab^2abab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2abab^2a-abab^2ababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2abab^2ab^2a-ababab^2abab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2abab^2abab-ab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ababab^2ab-abab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2abab^2ab-abab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abab^2ab-b^2abab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2ab^2ab-b^2ababab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2abab^2-bab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2ababa-b^2ab^2abab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababab^2aba-b^2abab^2ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2abab^2-bab^2ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2ab^2-babab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2aba-b^2abab^2abab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ababab^2a-bab^2ab^2abab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abababab-babababab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2abab^2a-bab^2abab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababababab^2ab-bab^2ababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abab^2ababab^2-ab^2ab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ababab^2abab^2-ab^2abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababababab^2-b^2ababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababab^2abab^2-b^2abab^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2ab^2ab^2ab-bab^2ab^2ab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2abab^2-b^2abab^2ab^2ab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abab^2ab^2ab^2ababa-ab^2ab^2abababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ababab^2ab^2ab^2aba-ab^2abababab^2ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2ababab^2a-abab^2ab^2abababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2abababab^2ab-abab^2ab^2ab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2abababab^2ab^2ab-ababab^2ab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abababab^2aba-b^2abab^2ab^2ab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababab^2ab^2ab^2ab-b^2abababab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababab^2ab^2aba-b^2ababab^2ab^2ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababab^2ababab-bab^2abababab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2ab^2ababab^2-bab^2ab^2abababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2abababab^2a-bab^2ab^2ab^2ababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ab^2ab^2abab^2-bab^2abababab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ababababababab^2-ababab^2ab^2ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababababababab^2a-b^2ababab^2ab^2ababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2ab^2abababab-babababab^2ab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ababab^2ab^2ababab^2a-ab^2ababababababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababab^2ab^2ababab^2-b^2ababababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababababababababab-babababababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abababab^2ab^2abababa-ab^2ab^2ab^2ab^2ab^2ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2ab^2ab^2ab^2ab-abababab^2ab^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2ab^2ab^2aba-b^2abababab^2ab^2ababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2abababab^2-b^2abababab^2ab^2ab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababab^2ab^2ab^2ab^2-b^2ab^2ab^2ab^2abababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababababababababab-bababababababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababab^2ab^2abababab^2-bab^2ab^2ab^2ab^2ab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abababab^2ab^2ab^2ababa-abababab^2ab^2ab^2ababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babababab^2ab^2ab^2ababab-abababab^2ab^2ab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2ab^2abababa-bab^2ab^2ab^2abababababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abababababababab-bababab^2ab^2abababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababab^2ab^2ab^2abab-bababab^2ab^2ab^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababab^2ab^2ab^2abababa-b^2abababab^2ab^2ab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababab^2ab^2ababababab-b^2ab^2ababababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ab^2ab^2ab^2ab^2ababab-ab^2ab^2ab^2ababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2ab^2ab^2ababab^2-abab^2ab^2abababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2ababababababa-ab^2ab^2ab^2ab^2ab^2ababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ababababababab-babab^2ab^2abababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ab^2abababababab^2-b^2ab^2ab^2abababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2abababababab-bab^2ab^2ab^2ab^2abababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ababab^2ab^2ababababababa-ababab^2ab^2ababababababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bababab^2ab^2ababababababab-ababab^2ab^2ababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababab^2ab^2abababababab-babab^2ab^2abababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ab^2abababababababa-b^2ababab^2ab^2abababababab"),
			},
			complete: true,
		},
		// G13, Example 4.2.26 Xiu Xingqiang, expected basis from bergman.
		{
			ideal: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "(ababababab^2ab^2)^2-1"),
			},
			maxiter: 7823,
			basis: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2abababab-ab^2ab^2ab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2ab^2ababab-ab^2ab^2ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bababab^2ab^2abab-ab^2ababab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babababab^2ab^2ab-ababab^2ab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ababababa-b^2ab^2ab^2ab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababababab^2a-bab^2ab^2ab^2ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2abababa-b^2ab^2ab^2ababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ab^2ababa-b^2ab^2ababab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababab^2ab^2aba-b^2ababab^2ab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababababab^2ab^2a-babab^2ab^2ab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ab^2ab^2ababa-ab^2ab^2abababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ab^2ababab^2a-abab^2ab^2ababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ababab^2ab^2a-ababab^2ab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ababababab^2-ab^2ab^2ab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ababab^2ab^2ab^2a-abababab^2ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ababababab^2ab^2-abab^2ab^2ab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2ab^2aba-ab^2ababababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2ab^2ab^2ab^2a-ababababab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2abab-b^2ab^2ababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ababab^2-bab^2ab^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ababab^2ab^2-babab^2ab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababab^2ab^2ab^2-bababab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2ab^2ab^2ab-b^2ababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ab^2ab^2ab^2-babababab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababab^2abababab-babababab^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababababab^2ababab-bababab^2ababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ab^2abab^2ababa-ab^2ab^2abab^2ababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2abab^2ababab^2a-abab^2ab^2abab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ababab^2abab^2a-abab^2abab^2ab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abab^2ab^2ab^2ababa-ab^2ab^2abababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abab^2ababab^2ab^2a-ababab^2ab^2abab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ababab^2abab^2ab^2a-ababab^2abab^2ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2abab^2aba-ab^2abab^2abababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2abab^2ababab-ab^2ab^2abab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2abab^2ab^2ab^2aba-ab^2abababab^2abab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2ab^2ab^2abab^2a-abab^2abababab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2ab^2abab^2abab-ab^2abab^2ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2abab^2ab^2ab^2a-abababab^2abab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2abab^2ab^2abab-ab^2ababab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2abab^2abababa-ab^2ababababababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bababab^2abab^2ab^2ab-ababab^2abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bababab^2abab^2ababa-abababababababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2abab^2abab-b^2ab^2abab^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abab^2ababab^2-bab^2ab^2abab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abab^2abababa-b^2ab^2ab^2abab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ababab^2abab^2-bab^2abab^2ab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abababab^2aba-b^2abab^2ab^2ab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2ab^2ab^2abab-b^2ab^2abababab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2ababab^2ab^2-babab^2ab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2abababab^2a-bab^2ab^2ab^2abab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababab^2abab^2ab^2-babab^2abab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababab^2abab^2a-bab^2abab^2ab^2ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababab^2ababab-bababab^2abababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababababababa-babab^2abab^2ababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2ab^2abab^2ab-b^2abab^2abababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2abab^2ababa-b^2ab^2abab^2ababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2ab^2ab^2ab-b^2abababab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2ab^2ababa-b^2ab^2ababab^2abab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abababab^2ab^2a-babab^2ab^2ab^2abab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ab^2ab^2abab^2-bab^2abababab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ab^2abab^2aba-b^2abab^2ababab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2abab^2ab^2ab^2-bababab^2abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2abab^2ab^2aba-b^2ababab^2abab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2abab^2ababab-b^2abababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ababababab^2-b^2ababababab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababab^2abab^2ab^2a-babab^2abab^2ab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababab^2abab^2abab-bababababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababababababab^2a-bababab^2abab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2abab^2abababab^2-ab^2ab^2ab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2abababab^2abab^2-ab^2abab^2ab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ababababababab-abab^2abab^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abab^2abababab^2ab^2-abab^2ab^2ab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abababab^2abab^2ab^2-abab^2abab^2ab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abababab^2abab^2aba-ababababababab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abababababababab^2-ababab^2abab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bababababababab^2ab^2-abababab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2abababab^2-b^2ab^2ababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2abababab^2ab^2-b^2ab^2abababab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ababab^2ababab^2-bab^2ab^2abab^2ab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababab^2ababab^2ab^2-babab^2ab^2abab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2abab^2ab^2abab-b^2ab^2ababab^2ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ab^2abab^2ab^2ab-b^2ababab^2ababab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2ab^2abab^2aba-b^2ababababababababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2ababababab^2-b^2ababababab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababab^2ab^2ab^2ababab-bababab^2ab^2ab^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2abab^2ab^2ab^2ab-bab^2ab^2ab^2abab^2ab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababab^2ababab^2ababab-bababab^2ababab^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2ab^2ab^2ab^2ab^2a-ab^2ab^2abab^2abab^2ab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2abab^2ab^2ababa-ab^2ab^2ababab^2abababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bababab^2abab^2abab^2ababa-ababab^2ababab^2abababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2ab^2ab^2ab^2ab-b^2ab^2ab^2abab^2abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ababab^2abababab^2a-bab^2ab^2ab^2abab^2ab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2abab^2abab^2ababab-babababab^2ababab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ababab^2abababab^2-b^2abababab^2ababab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ababab^2ababababa-bababab^2abab^2abab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ab^2abab^2abab^2ab^2ab^2-ab^2ab^2ab^2ab^2ab^2ab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ababab^2ababab^2abababab-ababab^2abab^2abab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abababab^2ababab^2ab^2aba-ababab^2abababababab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2abab^2abab^2ababab^2a-ab^2ababab^2ababab^2ababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bababab^2abababababab^2ab^2-abababab^2ababab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2ababab^2abab^2abab-babab^2abab^2ababab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababab^2abab^2abab^2abab-bababab^2ababab^2ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababab^2ababab^2abababa-babab^2abab^2abab^2ababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2abab^2abab^2ab^2ab^2-b^2ab^2ab^2abab^2abab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2ababab^2abababab^2-b^2ab^2abababababab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ababab^2ababab^2ababab-abab^2abab^2abab^2ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bababab^2abababababab^2abab-ab^2ababab^2abababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2abababababab^2ababa-b^2ab^2ababab^2abababababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ababab^2abababababab^2a-ababab^2abababababab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ababab^2abababababab^2-babab^2abababababab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababab^2ababab^2abab^2abab-babab^2abab^2ababab^2ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ab^2abab^2ab^2abab^2ab^2aba-abab^2ab^2abab^2abab^2abab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2abab^2ab^2abab^2ab^2ab^2a-ab^2ab^2abab^2abab^2abab^2ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2ababab^2abababab^2abab^2-ab^2ababab^2abab^2abab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2abab^2ab^2abab^2ab^2ab-bab^2ab^2abab^2abab^2abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abab^2ab^2abab^2ab^2ab^2ab-bab^2ab^2ab^2abab^2ab^2abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abab^2abab^2abab^2ab^2aba-bab^2ab^2abab^2ab^2abab^2ab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2abab^2abab^2abab^2ab^2a-b^2ab^2ab^2abab^2ab^2abab^2ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2abab^2ababab^2abab^2-b^2abab^2ababab^2abab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ab^2ab^2abab^2ab^2abab^2ab^2a-ab^2ab^2abab^2abab^2abab^2ab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ab^2abab^2abab^2abab^2ab^2ab-ab^2ab^2abab^2ab^2abab^2ab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ababab^2ababababababab^2a-abab^2abab^2abababababab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abab^2ababab^2abababababab^2a-ababab^2abababababab^2abab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2abab^2abab^2abab^2ab^2ab^2-ab^2ab^2ab^2abab^2ab^2abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bababab^2ababababababab^2ababa-ababababababababababababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ababab^2ababababababab^2-bab^2abab^2abababababab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2ababab^2abababababab^2-babab^2abababababab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2abababababab^2ababa-b^2ab^2ababab^2ababababababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2abab^2abab^2abab^2ababab-bababab^2abab^2abab^2abab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2abababababab^2abab^2aba-b^2abab^2ababab^2abababababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ab^2abab^2abab^2abab^2ababab-abab^2ab^2ababab^2ababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ababababababababababababab^2-ababab^2ababababababab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bababab^2abab^2abab^2abab^2ab^2ab^2-abab^2ababab^2abababababab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2abab^2abab^2ababab^2abab^2-b^2abab^2ababab^2abab^2abab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2abab^2abab^2abab^2abab^2ababab-abab^2abab^2ababababababab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2abab^2ababab^2ababab^2ababab-ab^2abab^2ababab^2ababab^2ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bababab^2abab^2abab^2abab^2abab^2ab^2-ababab^2ababababababab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abab^2abab^2abab^2abab^2ab^2ab-bab^2ab^2abab^2abab^2abab^2abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abab^2abab^2abab^2abab^2ababa-bab^2abab^2abab^2ababab^2ababababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2abab^2abab^2abab^2abab^2ab^2-b^2ab^2abab^2abab^2abab^2abab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2abababababab^2abab^2aba-b^2abababababababababababababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2abab^2abab^2abab^2abab^2ab^2a-b^2ababab^2ababababababab^2abab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ab^2abab^2abab^2abab^2abab^2abab-ab^2abab^2abab^2ababab^2abababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ab^2abab^2abab^2ababab^2abababa-ab^2ab^2abab^2abab^2ababab^2abababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2abab^2abab^2ababab^2abababab^2a-abab^2ab^2abab^2abab^2ababab^2ababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ababab^2ababababababab^2abab^2aba-ababab^2abab^2abab^2abab^2abab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abababab^2ababab^2abab^2abab^2ab^2a-abababab^2ababab^2abab^2abab^2ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2abab^2abab^2ababab^2abab^2aba-ab^2abab^2ababababababab^2abab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2abab^2abab^2ababab^2abababab-ab^2ab^2abab^2abab^2ababab^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2abab^2abab^2ababab^2abababababa-ab^2ab^2abab^2abab^2abab^2abab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2abab^2abab^2abab^2abab^2ab^2ab^2-ab^2ababab^2ababababababab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babababab^2ababab^2abab^2abab^2ab^2ab-abababab^2ababab^2abab^2abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abab^2abab^2ababab^2abababab^2-bab^2ab^2abab^2abab^2ababab^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abab^2abab^2ababab^2ababababa-b^2ab^2ab^2abab^2abab^2ababab^2ababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababab^2ababababababab^2abab^2ab-babab^2abab^2abab^2abab^2abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2abab^2abab^2ababab^2abababa-b^2ab^2abab^2abab^2ababab^2abababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2abab^2ababab^2ababababab-b^2ab^2abab^2abab^2abab^2abab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababab^2ababab^2abab^2abab^2ab^2ab^2-babababab^2ababab^2abab^2abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababab^2ababab^2abab^2abab^2ab^2aba-b^2abababab^2ababab^2abab^2abab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ab^2abab^2abab^2ababab^2ababab^2a-abab^2ab^2abab^2abab^2abab^2abab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2abab^2ababababababab^2abab^2ab^2-ab^2ab^2abab^2abab^2ababab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2ababababababab^2abab^2ab^2-b^2ab^2abab^2ababababababab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2abab^2abab^2abab^2abab^2ababab-b^2abababababababababababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababababababababababababababab^2-bababab^2abab^2abab^2abab^2abab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abab^2ababab^2ababab^2ababab^2abab^2a-abab^2abab^2abab^2ababab^2ababab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2ababab^2ababab^2ababab^2abab^2-bab^2abab^2abab^2ababab^2ababab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2ababab^2ababababababab^2abab-babab^2abab^2ababab^2ababababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2abab^2ababab^2ababab^2ababa-b^2abab^2ababab^2ababab^2ababab^2abab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ababababababababababababababababa-ababab^2ababab^2ababab^2ababab^2ababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2abab^2abab^2abab^2abab^2abab^2ababab-ab^2ab^2abab^2abab^2abab^2abab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2abab^2abab^2abab^2abab^2abab^2abab-ab^2abab^2abab^2abab^2abab^2abab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bababab^2abab^2abab^2abab^2abab^2abab^2ab-abab^2abab^2abab^2abab^2abab^2abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bababab^2ababab^2ababab^2ababab^2ababab-ababababababababababababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2abab^2abab^2abab^2abab^2ababa-b^2ab^2abab^2abab^2abab^2abab^2abab^2abab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2abab^2abab^2abab^2abab^2abab^2aba-b^2abab^2abab^2abab^2abab^2abab^2abab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2abab^2abab^2abab^2abab^2abab^2abab^2a-abab^2abab^2abab^2abab^2abab^2abab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abab^2abab^2abab^2abab^2abab^2abab^2ab^2a-ababab^2abab^2abab^2abab^2abab^2abab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abab^2abab^2abab^2abab^2abab^2abab^2-bab^2abab^2abab^2abab^2abab^2abab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abab^2abab^2abab^2abab^2abab^2ababa-bab^2abab^2ababab^2ababab^2ababab^2ababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2abab^2abab^2abab^2abab^2abab^2ab^2-babab^2abab^2abab^2abab^2abab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2ababab^2ababab^2ababab^2ababab^2-bab^2ab^2abab^2abab^2ababab^2ababab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababab^2ababab^2ababab^2ababab^2abab^2-bab^2abab^2abab^2ababab^2ababab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2abab^2abab^2ababab^2ababab^2abab-b^2abab^2ababab^2ababab^2ababab^2ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2abab^2ababab^2ababab^2abab^2ab-b^2ababab^2ababab^2ababab^2ababab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2ababab^2ababab^2ababab^2ababa-bab^2ab^2abab^2abab^2abab^2abab^2abab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2abab^2abab^2abab^2abab^2abab^2ab^2a-bababab^2ababab^2ababab^2ababab^2abab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ababab^2ababab^2ababab^2abab^2aba-babab^2abab^2abab^2abab^2abab^2abab^2ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ababab^2ababab^2ababab^2ababab^2abab^2ab-abab^2abab^2abab^2abab^2abab^2abab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2abab^2abab^2abab^2abab^2abab^2ab^2aba-ababab^2ababab^2ababab^2ababab^2abab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2abab^2abab^2abab^2abab^2abab^2abab-bab^2abab^2ababab^2ababab^2ababab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2abab^2abab^2abab^2abab^2abab^2ab-bab^2abab^2abab^2abab^2abab^2abab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababababababababababababababababababab-babababababababababababababababababababa"),
			},
			complete: true,
			long:     true,
		},
	}

	for i, test := range tests {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			t.Parallel()
			if testing.Short() && test.long {
				t.Skip()
			}
			ideal := make([]*Polynomial, len(test.ideal))
			copy(ideal, test.ideal)

			basis, complete := Buchberger(ideal, test.maxiter)
			if len(basis) != len(test.basis) {
				t.Fatalf("%d %v", len(basis), basis)
			}
			for i := range basis {
				if basis[i].Cmp(test.basis[i]) != 0 {
					t.Errorf("%d %v", i, basis[i])
				}
			}
			if complete != test.complete {
				t.Errorf("got %v want %v", complete, test.complete)
			}
		})
	}
}

func TestAddObstructions(t *testing.T) {
	tests := []struct {
		b   []obstruction
		g   []*Polynomial
		obs []obstruction
	}{
		// OBS(2), Example 5.12, Mora.
		{
			b: []obstruction{
				{i: 0, j: 0, iLeft: Monomial{1, 2}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2, 1}},
			},
			g: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 2, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 1, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2}},
				),
			},
			obs: []obstruction{
				{i: 0, j: 1, iLeft: Monomial{2}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{1}},
				{i: 0, j: 1, iLeft: Monomial{}, iRight: Monomial{2}, jLeft: Monomial{1}, jRight: Monomial{}},
			},
		},
		// OBS(4), Example 5.12, Mora.
		{
			b: []obstruction{
				{i: 2, j: 2, iLeft: Monomial{2}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2}},
				{i: 1, j: 2, iLeft: Monomial{2}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{1, 2}},
				{i: 1, j: 2, iLeft: Monomial{}, iRight: Monomial{2}, jLeft: Monomial{2, 1}, jRight: Monomial{}},
			},
			g: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 2, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 1, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1, 2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1, 2}},
				),
			},
			obs: []obstruction{
				{i: 2, j: 2, iLeft: Monomial{2}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2}},
				{i: 0, j: 3, iLeft: Monomial{}, iRight: Monomial{}, jLeft: Monomial{1}, jRight: Monomial{}},
				{i: 1, j: 3, iLeft: Monomial{}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2}},
				{i: 2, j: 3, iLeft: Monomial{}, iRight: Monomial{1}, jLeft: Monomial{2}, jRight: Monomial{}},
			},
		},
		// Example 4.15, Clemens Hofstadler.
		{
			b: []obstruction{
				{i: 0, j: 1, iLeft: Monomial{}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2, 1}},
				{i: 0, j: 1, iLeft: Monomial{}, iRight: Monomial{}, jLeft: Monomial{2}, jRight: Monomial{1}},
				{i: 2, j: 3, iLeft: Monomial{1, 2}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2}},
				{i: 2, j: 3, iLeft: Monomial{}, iRight: Monomial{1}, jLeft: Monomial{}, jRight: Monomial{}},
			},
			g: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 2, 1}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 2, 1}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 1}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(2, 1), Monomial: Monomial{1}},
				),
			},
			obs: []obstruction{
				{i: 2, j: 3, iLeft: Monomial{1, 2}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2}},
				{i: 2, j: 3, iLeft: Monomial{}, iRight: Monomial{1}, jLeft: Monomial{}, jRight: Monomial{}},
				{i: 0, j: 4, iLeft: Monomial{}, iRight: Monomial{}, jLeft: Monomial{2, 2}, jRight: Monomial{}},
				{i: 2, j: 4, iLeft: Monomial{}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2}},
			},
		},
	}

	for i, test := range tests {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			t.Parallel()
			b := make([]obstruction, len(test.b))
			copy(b, test.b)
			buf := &Monomial{}
			obs := addObstructions(b, test.g, buf)

			if len(obs) != len(test.obs) {
				t.Fatalf("%v", obs)
			}
			for i := range obs {
				if !obsEq(obs[i], test.obs[i]) {
					t.Errorf("%d %v", i, obs[i])
				}
			}
		})
	}
}

func TestRemove4d(t *testing.T) {
	tests := []struct {
		b       []obstruction
		sPObs   []obstruction
		g       []*Polynomial
		removed []bool
	}{
		// OBS(2), Example 5.12, Mora.
		{
			b: []obstruction{
				{i: 0, j: 0, iLeft: Monomial{1, 2}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2, 1}},
			},
			sPObs: []obstruction{
				{i: 0, j: 1, iLeft: Monomial{}, iRight: Monomial{2}, jLeft: Monomial{1}, jRight: Monomial{}},
				{i: 0, j: 1, iLeft: Monomial{2}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{1}},
			},
			g: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 2, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 1, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2}},
				),
			},
			removed: []bool{true},
		},
		// OBS(4), Example 5.12, Mora.
		{
			b: []obstruction{
				{i: 2, j: 2, iLeft: Monomial{2}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2}},
				{i: 1, j: 2, iLeft: Monomial{2}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{1, 2}},
				{i: 1, j: 2, iLeft: Monomial{}, iRight: Monomial{2}, jLeft: Monomial{2, 1}, jRight: Monomial{}},
			},
			sPObs: []obstruction{
				{i: 0, j: 3, iLeft: Monomial{}, iRight: Monomial{}, jLeft: Monomial{1}, jRight: Monomial{}},
				{i: 1, j: 3, iLeft: Monomial{}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2}},
				{i: 2, j: 3, iLeft: Monomial{}, iRight: Monomial{1}, jLeft: Monomial{2}, jRight: Monomial{}},
			},
			g: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 2, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 1, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1, 2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1, 2}},
				),
			},
			removed: []bool{false, true, true},
		},
		// Example 3.13, Kreuzer.
		{
			b: []obstruction{
				{i: 0, j: 1, iLeft: Monomial{}, iRight: Monomial{}, jLeft: Monomial{1}, jRight: Monomial{2, 1}},
			},
			sPObs: []obstruction{
				{i: 0, j: 2, iLeft: Monomial{}, iRight: Monomial{}, jLeft: Monomial{1, 1, 1, 2}, jRight: Monomial{}},
				{i: 1, j: 2, iLeft: Monomial{}, iRight: Monomial{}, jLeft: Monomial{1}, jRight: Monomial{}},
				{i: 1, j: 2, iLeft: Monomial{}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{1}},
			},
			g: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 1, 1, 2, 1}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 1}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1}},
				),
			},
			removed: []bool{true},
		},
		{
			b: []obstruction{
				{i: 0, j: 1, iLeft: Monomial{1}, iRight: Monomial{2, 1}, jLeft: Monomial{}, jRight: Monomial{}},
			},
			sPObs: []obstruction{
				{i: 0, j: 2, iLeft: Monomial{}, iRight: Monomial{}, jLeft: Monomial{1}, jRight: Monomial{}},
				{i: 0, j: 2, iLeft: Monomial{}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{1}},
				{i: 1, j: 2, iLeft: Monomial{}, iRight: Monomial{}, jLeft: Monomial{1, 1, 1, 2}, jRight: Monomial{}},
			},
			g: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 1}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 1, 1, 2, 1}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1}},
				),
			},
			removed: []bool{true},
		},
	}

	for i, test := range tests {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			t.Parallel()
			buf := &Monomial{}
			remove4d(test.b, test.sPObs, test.g, buf)
			for i := range test.b {
				if test.b[i].removed != test.removed[i] {
					t.Errorf("%d", i)
				}
			}
		})
	}
}

func TestRemove4c(t *testing.T) {
	tests := []struct {
		sPObs   []obstruction
		b       []obstruction
		ltgs    Monomial
		removed []obstruction
	}{
		// Example 3.10, Kreuzer.
		{
			sPObs: []obstruction{
				{i: 1, j: 3, iLeft: Monomial{}, iRight: Monomial{1}, jLeft: Monomial{}, jRight: Monomial{}},
				{i: 1, j: 3, iLeft: Monomial{1, 2}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2}},
				{i: 2, j: 3, iLeft: Monomial{}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2}},
				{i: 2, j: 3, iLeft: Monomial{1, 2}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2, 1, 2}},
				{i: 3, j: 3, iLeft: Monomial{}, iRight: Monomial{2, 1}, jLeft: Monomial{1, 2}, jRight: Monomial{}},
			},
			b: []obstruction{
				{i: 1, j: 2, iLeft: Monomial{1, 2}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{}},
				{i: 1, j: 2, iLeft: Monomial{}, iRight: Monomial{1, 2}, jLeft: Monomial{}, jRight: Monomial{}},
				{i: 2, j: 2, iLeft: Monomial{}, iRight: Monomial{1, 2}, jLeft: Monomial{1, 2}, jRight: Monomial{}},
			},
			ltgs: Monomial{1, 2, 1},
			removed: []obstruction{
				{i: 1, j: 3, iLeft: Monomial{}, iRight: Monomial{1}, jLeft: Monomial{}, jRight: Monomial{}},
				{i: 1, j: 3, iLeft: Monomial{1, 2}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2}},
				{i: 2, j: 3, iLeft: Monomial{}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2}},
				{i: 3, j: 3, iLeft: Monomial{}, iRight: Monomial{2, 1}, jLeft: Monomial{1, 2}, jRight: Monomial{}},
			},
		},
		{
			sPObs: []obstruction{
				{i: 2, j: 3, iLeft: Monomial{}, iRight: Monomial{1}, jLeft: Monomial{}, jRight: Monomial{}},
			},
			b: []obstruction{
				{i: 1, j: 2, iLeft: Monomial{}, iRight: Monomial{1, 2}, jLeft: Monomial{}, jRight: Monomial{}},
			},
			ltgs:    Monomial{1, 2, 1},
			removed: []obstruction{},
		},
	}

	for i, test := range tests {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			t.Parallel()
			sPObs := make([]obstruction, len(test.sPObs))
			copy(sPObs, test.sPObs)

			remove4c(sPObs, test.b, test.ltgs)
			sPObs = slices.DeleteFunc(sPObs, func(o obstruction) bool { return o.removed })

			if len(sPObs) != len(test.removed) {
				t.Fatalf("%v", sPObs)
			}
			for i := range sPObs {
				if !obsEq(sPObs[i], test.removed[i]) {
					t.Errorf("%d %v", i, sPObs[i])
				}
			}
		})
	}
}

func TestRemove4b(t *testing.T) {
	tests := []struct {
		obs     []obstruction
		order   Order
		removed []obstruction
	}{
		// OBS(2) Example 5.10, Mora.
		{
			obs: []obstruction{
				{i: 2, j: 2, iLeft: Monomial{2, 1}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{1, 2}},
				{i: 1, j: 2, iLeft: Monomial{}, iRight: Monomial{2}, jLeft: Monomial{1}, jRight: Monomial{}},
				{i: 1, j: 2, iLeft: Monomial{2}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{1}},
			},
			order: Deglex,
			removed: []obstruction{
				{i: 1, j: 2, iLeft: Monomial{}, iRight: Monomial{2}, jLeft: Monomial{1}, jRight: Monomial{}},
				{i: 1, j: 2, iLeft: Monomial{2}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{1}},
			},
		},
		// OBS(4) Example 5.10, Mora.
		{
			obs: []obstruction{
				{i: 1, j: 4, iLeft: Monomial{}, iRight: Monomial{}, jLeft: Monomial{1}, jRight: Monomial{}},
				{i: 1, j: 4, iLeft: Monomial{2}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2, 1}},
				{i: 2, j: 4, iLeft: Monomial{}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2}},
				{i: 2, j: 4, iLeft: Monomial{}, iRight: Monomial{1}, jLeft: Monomial{2, 1}, jRight: Monomial{}},
				{i: 3, j: 4, iLeft: Monomial{}, iRight: Monomial{1}, jLeft: Monomial{2}, jRight: Monomial{}},
			},
			order: Deglex,
			removed: []obstruction{
				{i: 1, j: 4, iLeft: Monomial{}, iRight: Monomial{}, jLeft: Monomial{1}, jRight: Monomial{}},
				{i: 2, j: 4, iLeft: Monomial{}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2}},
				{i: 3, j: 4, iLeft: Monomial{}, iRight: Monomial{1}, jLeft: Monomial{2}, jRight: Monomial{}},
			},
		},
		// OBS(5) Example 5.10, Mora.
		{
			obs: []obstruction{
				{i: 1, j: 5, iLeft: Monomial{1}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{1}},
				{i: 1, j: 5, iLeft: Monomial{}, iRight: Monomial{1, 2}, jLeft: Monomial{1, 2}, jRight: Monomial{}},
				{i: 2, j: 5, iLeft: Monomial{1, 1}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{1, 2}},
				{i: 3, j: 5, iLeft: Monomial{1, 1}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2}},
				{i: 4, j: 5, iLeft: Monomial{}, iRight: Monomial{1, 2}, jLeft: Monomial{2}, jRight: Monomial{}},
				{i: 4, j: 5, iLeft: Monomial{1, 1}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{1}},
			},
			order: Deglex,
			removed: []obstruction{
				{i: 1, j: 5, iLeft: Monomial{1}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{1}},
				{i: 3, j: 5, iLeft: Monomial{1, 1}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2}},
				{i: 4, j: 5, iLeft: Monomial{}, iRight: Monomial{1, 2}, jLeft: Monomial{2}, jRight: Monomial{}},
			},
		},
		// Example 3.10, Kreuzer.
		{
			obs: []obstruction{
				{i: 1, j: 2, iLeft: Monomial{1, 2}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{}},
				{i: 1, j: 2, iLeft: Monomial{}, iRight: Monomial{1, 2}, jLeft: Monomial{}, jRight: Monomial{}},
				{i: 2, j: 2, iLeft: Monomial{}, iRight: Monomial{1, 2}, jLeft: Monomial{1, 2}, jRight: Monomial{}},
			},
			order: Deglex,
			removed: []obstruction{
				{i: 1, j: 2, iLeft: Monomial{}, iRight: Monomial{1, 2}, jLeft: Monomial{}, jRight: Monomial{}},
			},
		},
		{
			obs: []obstruction{
				{i: 1, j: 1, iLeft: Monomial{1, 2, 3}, iRight: Monomial{}, jLeft: Monomial{1, 2}, jRight: Monomial{2, 1}},
				{i: 1, j: 1, iLeft: Monomial{1, 2}, iRight: Monomial{}, jLeft: Monomial{1, 2}, jRight: Monomial{2, 1}},
			},
			order: Deglex,
			removed: []obstruction{
				{i: 1, j: 1, iLeft: Monomial{1, 2}, iRight: Monomial{}, jLeft: Monomial{1, 2}, jRight: Monomial{2, 1}},
			},
		},
	}

	for i, test := range tests {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			t.Parallel()
			obs := make([]obstruction, len(test.obs))
			copy(obs, test.obs)

			remove4b(obs, test.order)
			obs = slices.DeleteFunc(obs, func(o obstruction) bool { return o.removed })

			if len(obs) != len(test.removed) {
				t.Fatalf("%v", obs)
			}
			for i := range obs {
				if !obsEq(obs[i], test.removed[i]) {
					t.Errorf("%d %v", i, obs[i])
				}
			}
		})
	}
}

func TestShrink(t *testing.T) {
	tests := []struct {
		o      obstruction
		shrunk obstruction
	}{
		// Example 4.2.3, Xiu Xingqiang.
		{
			o:      obstruction{i: 1, j: 2, iLeft: Monomial{1, 2, 1}, iRight: Monomial{}, jLeft: Monomial{1}, jRight: Monomial{1, 1, 2, 1, 2}},
			shrunk: obstruction{i: 1, j: 2, iLeft: Monomial{2, 1}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{1, 1, 2, 1, 2}},
		},
		// Example 3.5, Kreuzer.
		{
			o:      obstruction{i: 1, j: 2, iLeft: Monomial{1, 2, 1, 1}, iRight: Monomial{}, jLeft: Monomial{1, 2}, jRight: Monomial{2}},
			shrunk: obstruction{i: 1, j: 2, iLeft: Monomial{1, 1}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2}},
		},
		{
			o:      obstruction{i: 2, j: 3, iLeft: Monomial{2}, iRight: Monomial{2, 1, 2}, jLeft: Monomial{}, jRight: Monomial{1, 1, 1, 2}},
			shrunk: obstruction{i: 2, j: 3, iLeft: Monomial{2}, iRight: Monomial{2}, jLeft: Monomial{}, jRight: Monomial{1, 1}},
		},
		{
			o:      obstruction{i: 3, j: 4, iLeft: Monomial{2, 2, 1, 2}, iRight: Monomial{1, 1, 2}, jLeft: Monomial{2, 2, 2}, jRight: Monomial{2, 1, 2}},
			shrunk: obstruction{i: 3, j: 4, iLeft: Monomial{1, 2}, iRight: Monomial{1}, jLeft: Monomial{2}, jRight: Monomial{2}},
		},
	}

	for i, test := range tests {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			t.Parallel()
			o := obstruction{i: test.o.i, j: test.o.j, iLeft: make(Monomial, len(test.o.iLeft)), iRight: make(Monomial, len(test.o.iRight)), jLeft: make(Monomial, len(test.o.jLeft)), jRight: make(Monomial, len(test.o.jRight))}
			copy(o.iLeft, test.o.iLeft)
			copy(o.iRight, test.o.iRight)
			copy(o.jLeft, test.o.jLeft)
			copy(o.jRight, test.o.jRight)
			shrunk := shrink(o)
			if !obsEq(shrunk, test.shrunk) {
				t.Errorf("%v", shrunk)
			}
		})
	}
}

func TestHasOverlap(t *testing.T) {
	tests := []struct {
		obstruction obstruction
		im          Monomial
		jm          Monomial
		overlap     bool
	}{
		// Example 4.2.3, Xiu Xingqiang.
		{
			obstruction: obstruction{i: 1, j: 3, iLeft: Monomial{1, 2, 1}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{1, 2}},
			im:          Monomial{1, 2, 1, 2},
			jm:          Monomial{1, 2, 1, 1, 2},
			overlap:     true,
		},
		// Example 4.2.3, Xiu Xingqiang.
		{
			obstruction: obstruction{i: 1, j: 2, iLeft: Monomial{1, 2, 1}, iRight: Monomial{}, jLeft: Monomial{1}, jRight: Monomial{1, 1, 2, 1, 2}},
			im:          Monomial{1, 2, 1, 2},
			jm:          Monomial{2},
			overlap:     false,
		},
		// Example 4.2.5, Xiu Xingqiang.
		{
			obstruction: obstruction{i: 2, j: 3, iLeft: Monomial{1, 2, 1, 1}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2, 2}},
			im:          Monomial{2, 2, 2},
			jm:          Monomial{1, 2, 1, 1, 2},
			overlap:     true,
		},
		// Example 4.2.5, Xiu Xingqiang.
		{
			obstruction: obstruction{i: 1, j: 2, iLeft: Monomial{1, 2}, iRight: Monomial{2}, jLeft: Monomial{1, 2, 1, 1}, jRight: Monomial{}},
			im:          Monomial{1, 1, 2, 2},
			jm:          Monomial{2, 2, 2},
			overlap:     true,
		},
		// Case d, example 4.1.8, Xiu Xingqiang.
		{
			obstruction: obstruction{i: 1, j: 2, iLeft: Monomial{}, iRight: Monomial{5, 6, 5, 6, 5, 6, 2}, jLeft: Monomial{4, 6, 5, 6, 5, 6, 5}, jRight: Monomial{}},
			im:          Monomial{4, 6, 5, 6, 5, 6, 5},
			jm:          Monomial{5, 6, 5, 6, 5, 6, 2},
			overlap:     false,
		},
		// Case e, example 4.1.8, Xiu Xingqiang.
		{
			obstruction: obstruction{i: 1, j: 2, iLeft: Monomial{}, iRight: Monomial{6, 5, 6, 5, 6, 5, 6, 2}, jLeft: Monomial{4, 6, 5, 6, 5, 6, 5, 6}, jRight: Monomial{}},
			im:          Monomial{4, 6, 5, 6, 5, 6, 5},
			jm:          Monomial{5, 6, 5, 6, 5, 6, 2},
			overlap:     false,
		},
		// Case f, example 4.1.8, Xiu Xingqiang.
		{
			obstruction: obstruction{i: 1, j: 2, iLeft: Monomial{}, iRight: Monomial{6, 5, 6, 5, 6, 5, 6, 5, 6, 2}, jLeft: Monomial{4, 6, 5, 6, 5, 6, 5, 6, 5, 6}, jRight: Monomial{}},
			im:          Monomial{4, 6, 5, 6, 5, 6, 5},
			jm:          Monomial{5, 6, 5, 6, 5, 6, 2},
			overlap:     false,
		},
		// Case a, Example 2.9, Kreuzer.
		{
			obstruction: obstruction{i: 1, j: 1, iLeft: Monomial{}, iRight: Monomial{1, 1, 1}, jLeft: Monomial{1, 1, 1}, jRight: Monomial{}},
			im:          Monomial{1, 1},
			jm:          Monomial{1, 1},
			overlap:     false,
		},
		// Case a, Example 2.9, Kreuzer.
		{
			obstruction: obstruction{i: 1, j: 1, iLeft: Monomial{}, iRight: Monomial{2, 1, 1}, jLeft: Monomial{1, 1, 2}, jRight: Monomial{}},
			im:          Monomial{1, 1},
			jm:          Monomial{1, 1},
			overlap:     false,
		},
		// Case b, Example 2.9, Kreuzer.
		{
			obstruction: obstruction{i: 1, j: 2, iLeft: Monomial{}, iRight: Monomial{1, 1, 2}, jLeft: Monomial{1, 1, 1}, jRight: Monomial{}},
			im:          Monomial{1, 1},
			jm:          Monomial{1, 2},
			overlap:     false,
		},
		// Case b, Example 2.9, Kreuzer.
		{
			obstruction: obstruction{i: 1, j: 2, iLeft: Monomial{1, 2, 2}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2, 1, 1}},
			im:          Monomial{1, 1},
			jm:          Monomial{1, 2},
			overlap:     false,
		},
		// Case c, Example 2.9, Kreuzer.
		{
			obstruction: obstruction{i: 2, j: 2, iLeft: Monomial{}, iRight: Monomial{1, 1, 2}, jLeft: Monomial{1, 2, 1}, jRight: Monomial{}},
			im:          Monomial{1, 2},
			jm:          Monomial{1, 2},
			overlap:     false,
		},
		// Example 3.10, Kreuzer.
		{
			obstruction: obstruction{i: 1, j: 3, iLeft: Monomial{1, 2, 1, 2}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2, 1, 2}},
			im:          Monomial{1, 2},
			jm:          Monomial{1, 2, 1},
			overlap:     false,
		},
	}

	for i, test := range tests {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			t.Parallel()
			// Check if test.obstruction is a syzygy.
			wiw := append(append(append(Monomial{}, test.obstruction.iLeft...), test.im...), test.obstruction.iRight...)
			wjw := append(append(append(Monomial{}, test.obstruction.jLeft...), test.jm...), test.obstruction.jRight...)
			if !monomialEq(wiw, wjw) {
				t.Errorf("%v %v", wiw, wjw)
			}

			// Check overlap.
			overlap := hasOverlap(test.obstruction, test.im, test.jm)
			if overlap != test.overlap {
				t.Errorf("%v", overlap)
			}
		})
	}
}

func TestObsExample4_1_15_XiuXingQiang(t *testing.T) {
	g := []*Polynomial{
		NewPolynomial(
			Deglex,
			PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{4, 6, 5, 6, 5, 6, 5}},
			PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{3, 6, 5, 6, 5}},
			PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{4}},
			PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{3}},
		),
		NewPolynomial(
			Deglex,
			PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{5, 6, 5, 6, 5, 6, 2}},
			PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{5, 6, 5, 6, 2}},
			PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2}},
			PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1}},
		),
		NewPolynomial(
			Deglex,
			PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{4, 6, 1}},
			PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{3, 6, 2}},
		),
		NewPolynomial(
			Deglex,
			PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{4, 6, 5, 6, 1}},
			PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{3, 6, 5, 6, 2}},
		),
		NewPolynomial(
			Deglex,
			PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{4, 6, 5, 6, 5, 6, 1}},
			PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{3, 6, 5, 6, 5, 6, 2}},
		),
	}
	var obs []obstruction
	for l := 3; l <= len(g); l++ {
		obs = overlapObstruction(obs, g[:l])
	}
	if len(obs) != 0 {
		t.Errorf("%v", obs)
	}
}

func TestIjObs(t *testing.T) {
	tests := []struct {
		g           []*Polynomial
		obstruction []obstruction
	}{
		// Example 4.1.8, Xiu Xingqiang.
		{
			g: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{4, 6, 5, 6, 5, 6, 5}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{3, 6, 5, 6, 5}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{4}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{3}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{5, 6, 5, 6, 5, 6, 2}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{5, 6, 5, 6, 2}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1}},
				),
			},
			obstruction: []obstruction{
				obstruction{i: 0, j: 1, iLeft: Monomial{}, iRight: Monomial{6, 2}, jLeft: Monomial{4, 6}, jRight: Monomial{}},
				obstruction{i: 0, j: 1, iLeft: Monomial{}, iRight: Monomial{6, 5, 6, 2}, jLeft: Monomial{4, 6, 5, 6}, jRight: Monomial{}},
				obstruction{i: 0, j: 1, iLeft: Monomial{}, iRight: Monomial{6, 5, 6, 5, 6, 2}, jLeft: Monomial{4, 6, 5, 6, 5, 6}, jRight: Monomial{}},
			},
		},
		// Example 5.10, Mora.
		{
			g: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 2, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 1, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1, 2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1, 2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 1, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2}},
				),
			},
			obstruction: []obstruction{
				obstruction{i: 0, j: 1, iLeft: Monomial{2}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{1}},
				obstruction{i: 0, j: 1, iLeft: Monomial{}, iRight: Monomial{2}, jLeft: Monomial{1}, jRight: Monomial{}},
				obstruction{i: 1, j: 2, iLeft: Monomial{2}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{1, 2}},
				obstruction{i: 1, j: 2, iLeft: Monomial{}, iRight: Monomial{2}, jLeft: Monomial{2, 1}, jRight: Monomial{}},
				obstruction{i: 0, j: 3, iLeft: Monomial{2}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2, 1}},
				obstruction{i: 0, j: 3, iLeft: Monomial{}, iRight: Monomial{}, jLeft: Monomial{1}, jRight: Monomial{}},
				obstruction{i: 1, j: 3, iLeft: Monomial{}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2}},
				obstruction{i: 1, j: 3, iLeft: Monomial{}, iRight: Monomial{1}, jLeft: Monomial{2, 1}, jRight: Monomial{}},
				obstruction{i: 2, j: 3, iLeft: Monomial{}, iRight: Monomial{1}, jLeft: Monomial{2}, jRight: Monomial{}},
				obstruction{i: 0, j: 4, iLeft: Monomial{1}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{1}},
				obstruction{i: 0, j: 4, iLeft: Monomial{}, iRight: Monomial{1, 2}, jLeft: Monomial{1, 2}, jRight: Monomial{}},
				obstruction{i: 1, j: 4, iLeft: Monomial{1, 1}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{1, 2}},
				obstruction{i: 2, j: 4, iLeft: Monomial{1, 1}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2}},
				obstruction{i: 3, j: 4, iLeft: Monomial{1, 1}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{1}},
				obstruction{i: 3, j: 4, iLeft: Monomial{}, iRight: Monomial{1, 2}, jLeft: Monomial{2}, jRight: Monomial{}},
			},
		},
		// Example 2.6, Green.
		{
			g: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 2, 3, 4, 1, 2, 3, 4}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1, 2, 5, 6, 7, 4}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{3, 4, 1, 2, 3, 4, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{5, 6, 7, 4, 1}},
				),
			},
			obstruction: []obstruction{
				obstruction{i: 0, j: 1, iLeft: Monomial{3, 4}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2, 3, 4}},
				obstruction{i: 0, j: 1, iLeft: Monomial{3, 4, 1, 2, 3, 4}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2, 3, 4, 1, 2, 3, 4}},
				obstruction{i: 0, j: 1, iLeft: Monomial{}, iRight: Monomial{1}, jLeft: Monomial{1, 2}, jRight: Monomial{}},
				obstruction{i: 0, j: 1, iLeft: Monomial{}, iRight: Monomial{1, 2, 3, 4, 1}, jLeft: Monomial{1, 2, 3, 4, 1, 2}, jRight: Monomial{}},
			},
		},
		{
			g: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 1, 1, 2, 2, 1, 1, 2, 2, 1, 2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 1, 2, 2, 1}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 1, 1, 2, 2, 1, 1, 2, 2, 1}},
				),
			},
			obstruction: []obstruction{
				obstruction{i: 0, j: 1, iLeft: Monomial{1, 1, 2, 2}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{1, 1, 2, 2, 1, 1, 2, 2, 1, 2}},
				obstruction{i: 0, j: 1, iLeft: Monomial{}, iRight: Monomial{}, jLeft: Monomial{1}, jRight: Monomial{1, 2, 2, 1, 2}},
				obstruction{i: 0, j: 1, iLeft: Monomial{}, iRight: Monomial{}, jLeft: Monomial{1, 1, 1, 2, 2}, jRight: Monomial{2}},
				obstruction{i: 0, j: 2, iLeft: Monomial{2, 1, 1, 2, 2, 1, 1, 2, 2}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{1, 1, 2, 2, 1, 1, 2, 2, 1, 2}},
				obstruction{i: 0, j: 2, iLeft: Monomial{}, iRight: Monomial{1, 1, 2, 2, 1, 1, 2, 2, 1}, jLeft: Monomial{1, 1, 1, 2, 2, 1, 1, 2, 2, 1}, jRight: Monomial{}},
				obstruction{i: 1, j: 2, iLeft: Monomial{2, 1, 1, 2, 2}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{}},
				obstruction{i: 1, j: 2, iLeft: Monomial{2, 1, 1, 2, 2, 1, 1, 2, 2}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{1, 2, 2, 1}},
				obstruction{i: 1, j: 2, iLeft: Monomial{}, iRight: Monomial{1, 2, 2, 1, 1, 2, 2, 1}, jLeft: Monomial{1, 1, 2}, jRight: Monomial{}},
				obstruction{i: 1, j: 2, iLeft: Monomial{2}, iRight: Monomial{1, 2, 2, 1}, jLeft: Monomial{}, jRight: Monomial{}},
			},
		},
	}

	for i, test := range tests {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			t.Parallel()
			var obs []obstruction
			for j := range test.g {
				for i := range j {
					obs = leftObstruction(obs, i, j, test.g)
					obs = rightObstruction(obs, i, j, test.g)
					obs = centerObstruction(obs, i, j, test.g)
				}
			}

			if len(obs) != len(test.obstruction) {
				t.Errorf("%v", obs)
			}
			for k, o := range obs {
				if !obsEq(o, test.obstruction[k]) {
					t.Errorf("%d %v", k, o)
				}
				im := test.g[o.i].LeadingTerm().Monomial
				jm := test.g[o.j].LeadingTerm().Monomial
				if !hasOverlap(o, im, jm) {
					t.Errorf("%d %v", k, o)
				}
				// Check that o is a syzygy.
				gi := test.g[o.i].LeadingTerm().Monomial
				gj := test.g[o.j].LeadingTerm().Monomial
				iwgw := append(append(append(Monomial{}, o.iLeft...), gi...), o.iRight...)
				jwgw := append(append(append(Monomial{}, o.jLeft...), gj...), o.jRight...)
				if !monomialEq(iwgw, jwgw) {
					t.Errorf("%v %v", iwgw, jwgw)
				}
			}
		})
	}
}

func TestIiObs(t *testing.T) {
	tests := []struct {
		g           []*Polynomial
		obstruction []obstruction
	}{
		// Example 5.10, Mora.
		{
			g: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 2, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 1, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1, 2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1, 2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 1, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2}},
				),
			},
			obstruction: []obstruction{
				obstruction{i: 0, j: 0, iLeft: Monomial{}, iRight: Monomial{2, 1}, jLeft: Monomial{1, 2}, jRight: Monomial{}},
				obstruction{i: 1, j: 1, iLeft: Monomial{}, iRight: Monomial{1, 2}, jLeft: Monomial{2, 1}, jRight: Monomial{}},
				obstruction{i: 2, j: 2, iLeft: Monomial{}, iRight: Monomial{2}, jLeft: Monomial{2}, jRight: Monomial{}},
			},
		},
		// Example 2.6, Green.
		{
			g: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 2, 3, 4, 1, 2, 3, 4}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1, 2, 5, 6, 7, 4}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{3, 4, 1, 2, 3, 4, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{5, 6, 7, 4, 1}},
				),
			},
			obstruction: []obstruction{
				obstruction{i: 0, j: 0, iLeft: Monomial{}, iRight: Monomial{1, 2, 3, 4}, jLeft: Monomial{1, 2, 3, 4}, jRight: Monomial{}},
				obstruction{i: 1, j: 1, iLeft: Monomial{}, iRight: Monomial{2, 3, 4, 1}, jLeft: Monomial{3, 4, 1, 2}, jRight: Monomial{}},
			},
		},
	}

	for i, test := range tests {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			t.Parallel()
			var obs []obstruction
			for i := range test.g {
				obs = rightObstruction(obs, i, i, test.g)
			}

			if len(obs) != len(test.obstruction) {
				t.Errorf("%v", obs)
			}
			for k, o := range obs {
				if !obsEq(o, test.obstruction[k]) {
					t.Errorf("%d %v", k, o)
				}
				im := test.g[o.i].LeadingTerm().Monomial
				jm := test.g[o.j].LeadingTerm().Monomial
				if !hasOverlap(o, im, jm) {
					t.Errorf("%d %v", k, o)
				}
				// Check that o is a syzygy.
				gi := test.g[o.i].LeadingTerm().Monomial
				gj := test.g[o.j].LeadingTerm().Monomial
				iwgw := append(append(append(Monomial{}, o.iLeft...), gi...), o.iRight...)
				jwgw := append(append(append(Monomial{}, o.jLeft...), gj...), o.jRight...)
				if !monomialEq(iwgw, jwgw) {
					t.Errorf("%v %v", iwgw, jwgw)
				}
			}
		})
	}
}

func TestInterreduce(t *testing.T) {
	tests := []struct {
		g       []*Polynomial
		reduced []*Polynomial
	}{
		{
			g: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-2, 1)},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2, 1}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 2, 2, 1}},
					PolynomialTerm{Coefficient: big.NewRat(3, 1), Monomial: Monomial{2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-4, 1)},
				),
			},
			reduced: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(2, 1), Monomial: Monomial{2, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-2, 1)},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1)},
				),
			},
		},
	}

	for i, test := range tests {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			t.Parallel()
			testg := make([]*Polynomial, 0, len(test.g))
			for _, gi := range test.g {
				testg = append(testg, NewPolynomial(Deglex).Set(gi))
			}
			reduced := interreduce(testg)

			if len(reduced) != len(test.reduced) {
				t.Errorf("%v", reduced)
			}
			for j := range reduced {
				if reduced[j].Cmp(test.reduced[j]) != 0 {
					t.Errorf("%d %v %v", j, reduced[j], test.reduced[j])
				}
			}
		})
	}
}

func TestDivide(t *testing.T) {
	tests := []struct {
		f         *Polynomial
		g         []*Polynomial
		quotient  [][]Quotient
		remainder *Polynomial
	}{
		// Example 3.2.2, Xiu Xingqiang.
		{
			f: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 3, 3, 2, 3}},
			),
			g: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{3, 2}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{3}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{3, 3}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{3, 1}},
				),
			},
			remainder: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 3, 1, 3}},
			),
		},
		// Example 3.2.5, Xiu Xingqiang.
		{
			f: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 3, 3, 2, 3}},
			),
			g: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{3, 3}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{3, 1}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{3, 2}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{3}},
				),
			},
			remainder: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1, 3, 1, 2, 3}},
			),
		},
		// Example 1.4, Mora.
		{
			f: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(2, 1), Monomial: Monomial{3, 2, 1}},
				PolynomialTerm{Coefficient: big.NewRat(-3, 1), Monomial: Monomial{1, 3, 3}},
			),
			g: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1, 2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{3, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1, 3}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{3, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2, 3}},
				),
			},
			remainder: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(-3, 1), Monomial: Monomial{1, 3, 3}},
				PolynomialTerm{Coefficient: big.NewRat(2, 1), Monomial: Monomial{1, 2, 3}},
			),
		},
		// Section 6.4 Simplifying Polynomial Expressions, NCAlgebra.
		{
			f: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 2, 1, 1}},
				PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1, 1, 2, 2}},
				PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 2, 1}},
			),
			g: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 2, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 2, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1, 2, 2}},
				),
			},
			remainder: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2}},
			),
		},
		// Section 6.5 Poynomials and Rules, NCAlgebra.
		{
			f: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 2, 1, 1}},
				PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 1, 2, 2, 2}},
				PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 2, 1}},
			),
			g: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1, 2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 1, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 1, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1}},
				),
			},
			remainder: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2}},
				PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 2}},
			),
		},
	}

	for i, test := range tests {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			t.Parallel()
			var remainder *Polynomial
			quotient := [][]Quotient{}
			testf := NewPolynomial(Deglex).Set(test.f)
			remainder, quotient = Divide(quotient, testf, test.g)

			if remainder.Cmp(test.remainder) != 0 {
				t.Errorf("%v", remainder)
			}

			// Check if quotient * g + remainder == f.
			f := NewPolynomial(test.f.order)
			for i := range quotient {
				for j := range quotient[i] {
					c := NewPolynomial(Deglex, PolynomialTerm{Coefficient: quotient[i][j].Coefficient})
					left := NewPolynomial(Deglex, PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: quotient[i][j].Left})
					right := NewPolynomial(Deglex, PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: quotient[i][j].Right})
					cwgw := mul(c, left, test.g[i], right)
					f.Add(f, cwgw)
				}
			}
			f.Add(f, remainder)
			if f.Cmp(test.f) != 0 {
				t.Errorf("%v", f)
			}
		})
	}
}

func TestPolynomialCmp(t *testing.T) {
	tests := []struct {
		x *Polynomial
		y *Polynomial
		c int
	}{
		{
			x: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 1, 2}},
				PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2}},
			),
			y: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 1}},
				PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1, 2}},
			),
			c: 1,
		},
		{
			x: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(9, 1), Monomial: Monomial{2}},
				PolynomialTerm{Coefficient: big.NewRat(7, 1), Monomial: Monomial{}},
			),
			y: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(3, 1), Monomial: Monomial{1, 2}},
				PolynomialTerm{Coefficient: big.NewRat(2, 1), Monomial: Monomial{2, 1}},
			),
			c: -1,
		},
	}

	for i, test := range tests {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			t.Parallel()
			if c := test.x.Cmp(test.y); c != test.c {
				t.Errorf("%d", c)
			}
		})
	}
}

func TestLeadingTerm(t *testing.T) {
	tests := []struct {
		x           *Polynomial
		leadingTerm PolynomialTerm
	}{
		{
			x: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(1, 2), Monomial: Monomial{1}},
				PolynomialTerm{Coefficient: big.NewRat(1, 3), Monomial: Monomial{2}},
				PolynomialTerm{Coefficient: big.NewRat(1, -4), Monomial: Monomial{1, 2}},
				PolynomialTerm{Coefficient: big.NewRat(1, -5), Monomial: Monomial{1, 1}},
			),
			leadingTerm: PolynomialTerm{Coefficient: big.NewRat(-1, 4), Monomial: Monomial{1, 2}},
		},
	}

	for i, test := range tests {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			t.Parallel()
			lt := test.x.LeadingTerm()
			if !termEq(lt, test.leadingTerm) {
				t.Errorf("%v", lt)
			}
		})
	}
}

func TestPolynomialAdd(t *testing.T) {
	type testcase struct {
		x *Polynomial
		y *Polynomial
		z *Polynomial
	}
	tests := []testcase{
		{
			x: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 1, 1}},
				PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 1, 1}},
				PolynomialTerm{Coefficient: big.NewRat(3, 1), Monomial: Monomial{1, 2}},
				PolynomialTerm{Coefficient: big.NewRat(4, 1), Monomial: Monomial{2, 1}},
				PolynomialTerm{Coefficient: big.NewRat(0, 1), Monomial: Monomial{2, 2}},
			),
			y: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(-3, 1), Monomial: Monomial{1, 2}},
				PolynomialTerm{Coefficient: big.NewRat(-7, 1), Monomial: Monomial{2, 1}},
				PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2, 2}},
			),
			z: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(2, 1), Monomial: Monomial{1, 1, 1}},
				PolynomialTerm{Coefficient: big.NewRat(-3, 1), Monomial: Monomial{2, 1}},
				PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2, 2}},
			),
		},
	}

	for i, test := range tests {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			t.Parallel()
			z := NewPolynomial(test.z.order)
			z.Add(test.x, test.y)
			if z.Cmp(test.z) != 0 {
				t.Errorf("%v", z)
			}

			// z = x.
			x := NewPolynomial(Deglex).Set(test.x)
			z = x
			z.Add(x, test.y)
			if z.Cmp(test.z) != 0 {
				t.Errorf("%v", z)
			}

			// z = y.
			y := NewPolynomial(Deglex).Set(test.y)
			z = y
			z.Add(test.x, y)
			if z.Cmp(test.z) != 0 {
				t.Errorf("%v", z)
			}
		})
	}
}

func TestPolynomialAddZEqXY(t *testing.T) {
	tests := []struct {
		x *Polynomial
		z *Polynomial
	}{
		{
			x: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(2, 1), Monomial: Monomial{1, 1, 1}},
				PolynomialTerm{Coefficient: big.NewRat(-3, 1), Monomial: Monomial{1, 2}},
			),
			z: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(4, 1), Monomial: Monomial{1, 1, 1}},
				PolynomialTerm{Coefficient: big.NewRat(-6, 1), Monomial: Monomial{1, 2}},
			),
		},
	}

	for i, test := range tests {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			t.Parallel()
			x := NewPolynomial(Deglex).Set(test.x)
			z := x
			z.Add(x, x)
			if z.Cmp(test.z) != 0 {
				t.Errorf("%v", z)
			}
		})
	}
}

func TestPolynomialMul(t *testing.T) {
	tests := []struct {
		x *Polynomial
		y *Polynomial
		z *Polynomial
	}{
		{
			x: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(2, 1), Monomial: Monomial{1}},
				PolynomialTerm{Coefficient: big.NewRat(3, 1), Monomial: Monomial{2}},
				PolynomialTerm{Coefficient: big.NewRat(-4, 1), Monomial: Monomial{1, 2}},
			),
			y: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(5, 2), Monomial: Monomial{1}},
				PolynomialTerm{Coefficient: big.NewRat(-6, 1), Monomial: Monomial{2}},
			),
			z: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(5, 1), Monomial: Monomial{1, 1}},
				PolynomialTerm{Coefficient: big.NewRat(-12, 1), Monomial: Monomial{1, 2}},
				PolynomialTerm{Coefficient: big.NewRat(15, 2), Monomial: Monomial{2, 1}},
				PolynomialTerm{Coefficient: big.NewRat(-18, 1), Monomial: Monomial{2, 2}},
				PolynomialTerm{Coefficient: big.NewRat(-10, 1), Monomial: Monomial{1, 2, 1}},
				PolynomialTerm{Coefficient: big.NewRat(24, 1), Monomial: Monomial{1, 2, 2}},
			),
		},
		{
			x: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(3, 1), Monomial: Monomial{1}},
				PolynomialTerm{Coefficient: big.NewRat(-4, 1), Monomial: Monomial{2}},
			),
			y: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(-2, 1)},
			),
			z: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(-6, 1), Monomial: Monomial{1}},
				PolynomialTerm{Coefficient: big.NewRat(8, 1), Monomial: Monomial{2}},
			),
		},
	}

	for i, test := range tests {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			t.Parallel()
			z := NewPolynomial(test.z.order).Mul(test.x, test.y)
			if z.Cmp(test.z) != 0 {
				t.Errorf("%v", z)
			}
		})
	}
}

func TestPow(t *testing.T) {
	tests := []struct {
		x *Polynomial
		y int
		z *Polynomial
	}{
		{
			x: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1}},
				PolynomialTerm{Coefficient: big.NewRat(-2, 1), Monomial: Monomial{2}},
			),
			y: 2,
			z: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(4, 1), Monomial: Monomial{2, 2}},
				PolynomialTerm{Coefficient: big.NewRat(-2, 1), Monomial: Monomial{2, 1}},
				PolynomialTerm{Coefficient: big.NewRat(-2, 1), Monomial: Monomial{1, 2}},
				PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 1}},
			),
		},
		{
			x: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1}},
				PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2}},
			),
			y: 3,
			z: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2, 2, 2}},
				PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 2, 1}},
				PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 1, 2}},
				PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2, 1, 1}},
				PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 2, 2}},
				PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1, 2, 1}},
				PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1, 1, 2}},
				PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 1, 1}},
			),
		},
	}

	for i, test := range tests {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			t.Parallel()
			z := NewPolynomial(Deglex)
			z.Pow(test.x, test.y)
			if z.Cmp(test.z) != 0 {
				t.Errorf("%v", z)
			}
		})
	}
}

func TestMulScalar(t *testing.T) {
	tests := []struct {
		x      *Polynomial
		scalar *big.Rat
		z      *Polynomial
	}{
		{
			x: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(-2, 3), Monomial: Monomial{2, 1, 2}},
				PolynomialTerm{Coefficient: big.NewRat(5, 1), Monomial: Monomial{2}},
			),
			scalar: big.NewRat(-6, 1),
			z: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(4, 1), Monomial: Monomial{2, 1, 2}},
				PolynomialTerm{Coefficient: big.NewRat(-30, 1), Monomial: Monomial{2}},
			),
		},
	}
	for i, test := range tests {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			t.Parallel()
			z := NewPolynomial(test.z.order)
			z.mulScalar(test.scalar, test.x)
			if z.Cmp(test.z) != 0 {
				t.Errorf("%v", z)
			}

			z = NewPolynomial(Deglex).Set(test.x)
			z.mulScalar(test.scalar, z)
			if z.Cmp(test.z) != 0 {
				t.Errorf("%v", z)
			}
		})
	}
}

func TestOrder(t *testing.T) {
	tests := []struct {
		words  []Monomial
		order  Order
		sorted []Monomial
	}{
		{
			words:  []Monomial{{2}, {1, 1, 2}, {1, 2}},
			order:  lexicographic,
			sorted: []Monomial{{1, 1, 2}, {1, 2}, {2}},
		},
		{
			words:  []Monomial{{1, 3, 4}, {1, 2, 4}, {1, 3, 6}, {3, 4, 6}, {2, 5, 6}, {1, 5, 6}, {2, 3, 4}, {3, 4, 5}, {1, 4, 5}, {1, 4, 6}, {2, 3, 6}, {2, 4, 5}, {1, 3, 5}, {1, 2, 5}, {2, 3, 5}, {1, 2, 3}, {4, 5, 6}, {3, 5, 6}, {2, 4, 6}, {1, 2, 6}},
			order:  lexicographic,
			sorted: []Monomial{{1, 2, 3}, {1, 2, 4}, {1, 2, 5}, {1, 2, 6}, {1, 3, 4}, {1, 3, 5}, {1, 3, 6}, {1, 4, 5}, {1, 4, 6}, {1, 5, 6}, {2, 3, 4}, {2, 3, 5}, {2, 3, 6}, {2, 4, 5}, {2, 4, 6}, {2, 5, 6}, {3, 4, 5}, {3, 4, 6}, {3, 5, 6}, {4, 5, 6}},
		},
		{
			words:  []Monomial{{2, 3}, {1, 3}, {1, 1}, {3, 3}, {1, 2}, {2, 2}},
			order:  lexicographic,
			sorted: []Monomial{{1, 1}, {1, 2}, {1, 3}, {2, 2}, {2, 3}, {3, 3}},
		},
		{
			words:  []Monomial{{'b', 'e', 'n', 'e', 'f', 'i', 't'}, {'b', 'e'}, {'b', 'e', 'n', 't'}, {'b', 'a', 'r', 'n', 'a', 'c', 'l', 'e'}, {'b', 'e', 'e', 'n'}},
			order:  lexicographic,
			sorted: []Monomial{{'b', 'a', 'r', 'n', 'a', 'c', 'l', 'e'}, {'b', 'e'}, {'b', 'e', 'e', 'n'}, {'b', 'e', 'n', 'e', 'f', 'i', 't'}, {'b', 'e', 'n', 't'}},
		},
		{
			words:  []Monomial{{1, 1}, {1, 1, 2}, {2, 1, 1}, {2, 2}, {1, 2, 2}, {2, 2, 2}, {}, {1, 2}, {1, 1, 1}, {2, 1, 2}, {2}, {1}, {2, 2, 1}, {1, 2, 1}, {2, 1}},
			order:  Deglex,
			sorted: []Monomial{{}, {1}, {2}, {1, 1}, {1, 2}, {2, 1}, {2, 2}, {1, 1, 1}, {1, 1, 2}, {1, 2, 1}, {1, 2, 2}, {2, 1, 1}, {2, 1, 2}, {2, 2, 1}, {2, 2, 2}},
		},
	}
	for i, test := range tests {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			t.Parallel()
			slices.SortFunc(test.words, test.order)
			if !slices.EqualFunc(test.words, test.sorted, monomialEq) {
				t.Errorf("%v", test.words)
			}
		})
	}
}

func TestMain(m *testing.M) {
	flag.Parse()
	log.SetFlags(log.Lmicroseconds | log.Llongfile | log.LstdFlags)

	m.Run()
}

func mul(x *Polynomial, y ...*Polynomial) *Polynomial {
	z := x
	for i := range y {
		z = NewPolynomial(z.order).Mul(z, y[i])
	}
	return z
}

func termEq(a, b PolynomialTerm) bool {
	if eq := monomialEq(a.Monomial, b.Monomial); !eq {
		return false
	}
	if a.Coefficient.Cmp(b.Coefficient) != 0 {
		return false
	}
	return true
}

func parseMust(variables map[string]Symbol, order Order, input string) *Polynomial {
	p, err := Parse(variables, order, input)
	if err != nil {
		panic(fmt.Sprintf("%+v", err))
	}
	return p
}
