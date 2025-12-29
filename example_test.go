package nag_test

import (
	"fmt"
	"math"

	"github.com/fumin/nag"
)

func Example() {
	// This example shows how to simplify the expression
	//
	//   bbaa - aabb + aba
	//
	// assuming
	//
	//   aba - b = 0
	//
	// where a, b are noncommutative variables.
	// Using Gröbner bases we will eventually see that bbaa-aabb+aba -> b.
	expr := "bbaa - aabb + aba"
	rules := []string{"aba - b"}

	// Compute the Gröbner basis using the Buchberger algorithm.
	variables := map[string]nag.Symbol{"a": 1, "b": 2}
	ideal := make([]*nag.Polynomial[*nag.Rat], len(rules))
	ideal[0], _ = nag.Parse(variables, nag.ElimOrder(), rules[0])
	basis, _ := nag.Buchberger(ideal, 50)
	// Print the Gröbner basis and notice that we have found an additional
	// useful relation: bba = abb
	fmt.Printf("Gröbner basis:\n")
	for _, b := range basis {
		fmt.Printf("  %v = 0\n", b)
	}
	fmt.Printf("\n")

	// Use the Gröbner basis to simplify the target expression.
	exprP, _ := nag.Parse(variables, nag.ElimOrder(), expr)
	_, simplified := nag.Divide(nil, exprP, basis)
	fmt.Printf("Simplified result: %v\n", simplified)

	// Output:
	// Gröbner basis:
	//   aba-b = 0
	//   b^2a-ab^2 = 0
	//
	// Simplified result: b
}

func Example_equation_solving() {
	// This example shows solving a set of equations using Gröbner bases.
	//
	// Given the set of equations below, we want to obtain an expression for
	// the variable G in terms of other variables D, X, A, B,... etc.
	//
	// The context behind this example is about the boundary value problem of differential equations.
	// In fact, the variable G represents the Green's function operator.
	// For more details, please see:
	// Rosenkranz, M., Buchberger, B., & Engl, H. W. (2003). Solving linear boundary value problems via non-commutative Gröbner bases. Applicable Analysis, 82(7), 655-675.
	equations := []string{
		"D^2GD^2 - D^2",
		"GD^2G - G",
		"GD^2 - 1 + (1-X)L + XR",
		"D^2G - 1",
		"DX - XD - 1",
		"DA - 1",
		"AD - 1 + L",
		"DB + 1",
		"BD - R + 1",
		"RX - R",
		"LX",
	}

	variables := map[string]nag.Symbol{"D": 1, "L": 2, "X": 3, "A": 4, "B": 5, "R": 6, "G": 7}
	ideal := make([]*nag.Polynomial[*nag.Rat], len(equations))
	for i, eq := range equations {
		ideal[i], _ = nag.Parse(variables, nag.ElimOrder(), eq)
	}
	basis, _ := nag.Buchberger(ideal, 50)
	solution := basis[len(basis)-1]
	fmt.Printf("Solution: %v = 0\n", solution)

	// Output:
	// Solution: G-XBX+XB-XAX+AX = 0
}

func Example_minimal_polynomial() {
	// This example shows how to find the minimal polynomial for α = √2+√3+√5.
	// The minimal polynomial, P(x), for an algebraic number α is one whose
	// coefficients are integers and which evaluates to zero when α is
	// passed in, i.e. P(α) = 0.
	ideal := []string{
		"x^2 - 2",
		"y^2 - 3",
		"z^2 - 5",
		"α - x - y - z",
		// Equations below express the fact that all variables commute.
		"xy - yx",
		"xz - zx",
		"xα - αx",
		"yz - zy",
		"yα - αy",
		"zα - αz",
	}

	// Compute the Gröbner basis.
	variables := map[string]nag.Symbol{"x": 4, "y": 3, "z": 2, "α": 1}
	idealP := make([]*nag.Polynomial[*nag.Rat], len(ideal))
	for i, p := range ideal {
		idealP[i], _ = nag.Parse(variables, nag.ElimOrder(), p)
	}
	basis, _ := nag.Buchberger(idealP, 50)
	fmt.Printf("Gröbner basis:\n")
	for _, b := range basis {
		fmt.Printf("  %v = 0\n", b)
	}
	fmt.Printf("\n")

	// The minimal polynomial for α is one that contains only α without other variables.
	// In this case, it is the first basis:
	//
	// α^8-40α^6+352α^4-960α^2+576
	//
	// Verify that the above polynomial is indeed zero when α = √2+√3+√5.
	minPoly := basis[0]
	α := math.Sqrt(2) + math.Sqrt(3) + math.Sqrt(5)
	var res float64
	for c, w := range minPoly.Terms() {
		cf, _ := c.Float64()
		pow := float64(len(w))
		res += cf * math.Pow(α, pow)
	}
	fmt.Printf("minPoly(α): %f\n", math.Abs(res))

	// Output:
	// Gröbner basis:
	//   α^8-40α^6+352α^4-960α^2+576 = 0
	//   z-5/576α^7+97/288α^5-95/36α^3+53/12α = 0
	//   y+1/96α^7-37/96α^5+61/24α^3-15/4α = 0
	//   x-1/576α^7+7/144α^5+7/72α^3-5/3α = 0
	//
	// minPoly(α): 0.000000
}

func ExampleElimOrder() {
	order := nag.ElimOrder()
	fmt.Println(order(nag.Monomial{2, 2, 2}, nag.Monomial{2, 2, 1, 1}))
	fmt.Println(order(nag.Monomial{2, 2}, nag.Monomial{2, 2, 1, 1}))
	fmt.Println(order(nag.Monomial{1, 1, 2, 2}, nag.Monomial{2, 2, 1, 1}))

	// Output:
	// 1
	// -1
	// -1
}

func ExampleDeglex() {
	fmt.Println(nag.Deglex(nag.Monomial{2, 2, 2}, nag.Monomial{2, 2, 1, 1}))
	fmt.Println(nag.Deglex(nag.Monomial{2, 2, 2}, nag.Monomial{2, 2, 1}))
	fmt.Println(nag.Deglex(nag.Monomial{2, 2, 2}, nag.Monomial{2, 3, 1}))

	// Output:
	// -1
	// 1
	// -1
}

func ExamplePolynomial_Terms() {
	p := nag.NewPolynomial(
		nag.NewRat(0, 1), nag.Deglex,
		nag.PolynomialTerm[*nag.Rat]{Coefficient: nag.NewRat(3, 2), Monomial: nag.Monomial{1, 1, 1}},
		nag.PolynomialTerm[*nag.Rat]{Coefficient: nag.NewRat(-1, 1), Monomial: nag.Monomial{2, 2}},
		nag.PolynomialTerm[*nag.Rat]{Coefficient: nag.NewRat(5, 1), Monomial: nag.Monomial{1, 2}},
	)
	for coeffi, monomial := range p.Terms() {
		fmt.Printf("coefficient: %s, monomial: %v\n", coeffi.RatString(), monomial)
	}

	// Output:
	// coefficient: 3/2, monomial: [1 1 1]
	// coefficient: -1, monomial: [2 2]
	// coefficient: 5, monomial: [1 2]
}

func ExamplePolynomial_LeadingTerm() {
	terms := []nag.PolynomialTerm[*nag.Rat]{
		{Coefficient: nag.NewRat(1, 2), Monomial: nag.Monomial{1, 1, 1}},
		{Coefficient: nag.NewRat(1, 3), Monomial: nag.Monomial{2, 2}},
	}

	p0 := nag.NewPolynomial(nag.NewRat(0, 1), nag.Deglex, terms...)
	fmt.Println(p0.LeadingTerm())

	p1 := nag.NewPolynomial(nag.NewRat(0, 1), nag.ElimOrder(), terms...)
	fmt.Println(p1.LeadingTerm())

	// Output:
	// {1/2 [1 1 1]}
	// {1/3 [2 2]}
}

func ExampleDivide() {
	variables := map[string]nag.Symbol{"x": 3, "y": 2, "z": 1}
	f, _ := nag.Parse(variables, nag.Deglex, "zx^2yx")
	g := make([]*nag.Polynomial[*nag.Rat], 2)
	g[0], _ = nag.Parse(variables, nag.Deglex, "xy + x")
	g[1], _ = nag.Parse(variables, nag.Deglex, "x^2 + xz")

	// Create a copy of f since nag.Divide modifies f upon return.
	fCopy := nag.NewPolynomial(nag.NewRat(0, 1), nag.Deglex).Set(f)
	var remainder *nag.Polynomial[*nag.Rat]
	quotient := make([][]nag.Quotient[*nag.Rat], 0)
	// Perfom the division.
	quotient, remainder = nag.Divide(quotient, fCopy, g)
	fmt.Println("remainder:", remainder)

	// Check that f = quotient*g + remainder.
	ff := nag.NewPolynomial(nag.NewRat(0, 1), nag.Deglex)
	ff.SymbolStringer = f.SymbolStringer
	cw := nag.NewPolynomial(nag.NewRat(0, 1), nag.Deglex)
	cwg := nag.NewPolynomial(nag.NewRat(0, 1), nag.Deglex)
	cwgw := nag.NewPolynomial(nag.NewRat(0, 1), nag.Deglex)
	for i := range quotient {
		for j := range quotient[i] {
			cij := nag.NewPolynomial(nag.NewRat(0, 1), nag.Deglex, nag.PolynomialTerm[*nag.Rat]{Coefficient: quotient[i][j].Coefficient})
			wij := nag.NewPolynomial(nag.NewRat(0, 1), nag.Deglex, nag.PolynomialTerm[*nag.Rat]{Monomial: quotient[i][j].Left})
			wPij := nag.NewPolynomial(nag.NewRat(0, 1), nag.Deglex, nag.PolynomialTerm[*nag.Rat]{Monomial: quotient[i][j].Right})
			cwgw.Mul(cwg.Mul(cw.Mul(cij, wij), g[i]), wPij)
			ff.Add(ff, cwgw)
		}
	}
	ff.Add(ff, remainder)
	fmt.Println("g*quotient + remainder:", ff, "==", f)

	// Output:
	// remainder: zxzx
	// g*quotient + remainder: zx^2yx == zx^2yx
}

func ExampleBuchberger() {
	ideal := []string{
		"aba - b",
		"bab - b",
	}

	// Run the Buchberger algorithm.
	variables := map[string]nag.Symbol{"a": 1, "b": 2}
	idealP := make([]*nag.Polynomial[*nag.Rat], len(ideal))
	idealP[0], _ = nag.Parse(variables, nag.Deglex, ideal[0])
	idealP[1], _ = nag.Parse(variables, nag.Deglex, ideal[1])
	basis, complete := nag.Buchberger(idealP, 10)

	// Print the computed Gröbner basis.
	fmt.Println("Gröbner basis:")
	fmt.Println("")
	for _, b := range basis {
		fmt.Println(" ", b)
	}
	fmt.Println("")
	fmt.Println("Basis is complete:", complete)

	// Output:
	// Gröbner basis:
	//
	//   ba-ab
	//   b^2-ab
	//   a^2b-b
	//
	// Basis is complete: true
}

func ExampleBuchbergerHomogeneous() {
	ideal := []string{
		"x^2 - 2y^2",
		"xy - 3z^2",
	}
	variables := map[string]nag.Symbol{"x": 1, "y": 2, "z": 3}
	idealP := make([]*nag.Polynomial[*nag.Rat], len(ideal))
	for i := range ideal {
		idealP[i], _ = nag.Parse(variables, nag.Deglex, ideal[i])
	}

	// Run the homogeneous Buchberger algorithm and truncate at degree 3.
	var maxDeg int = 3
	basis3, complete := nag.BuchbergerHomogeneous(idealP, maxDeg)
	fmt.Printf("Gröbner basis truncated at degree %d:\n", maxDeg)
	for _, b := range basis3 {
		fmt.Println(" ", b)
	}
	fmt.Println("Basis is complete:", complete)
	fmt.Println("")

	// Run Buchberger again and truncate at a higher degree.
	// We will get the full basis this time round, since maxDeg is higher than the basis' maximum degree.
	maxDeg = 5
	basis, complete := nag.BuchbergerHomogeneous(idealP, maxDeg)
	fmt.Printf("Gröbner basis truncated at degree %d:\n", maxDeg)
	// Skip printing bases that have been printed above.
	fmt.Printf("  ...\n")
	for _, b := range basis[len(basis3):] {
		fmt.Println(" ", b)
	}
	fmt.Println("Basis is complete:", complete)

	// Output:
	// Gröbner basis truncated at degree 3:
	//   y^2-1/2x^2
	//   z^2-1/3xy
	//   yx^2-x^2y
	//   zxy-xyz
	// Basis is complete: false
	//
	// Gröbner basis truncated at degree 5:
	//   ...
	//   zx^3-2xyzy
	// Basis is complete: true
}

func ExampleParse() {
	pStr := "-x^2y^3 + 5/3(y-x)x"

	variables := map[string]nag.Symbol{"x": 1, "y": 2}
	p, err := nag.Parse(variables, nag.Deglex, pStr)
	if err != nil {
		fmt.Println("error:", err)
		return
	}

	for coefficient, monomial := range p.Terms() {
		fmt.Printf("coefficient: %s, monomial: %v\n", coefficient.RatString(), monomial)
	}

	// Output:
	// coefficient: -1, monomial: [1 1 2 2 2]
	// coefficient: 5/3, monomial: [2 1]
	// coefficient: -5/3, monomial: [1 1]
}
