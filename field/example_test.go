package field_test

import (
	"cmp"
	"fmt"
	"math/big"
	"slices"

	"github.com/fumin/nag"
	"github.com/fumin/nag/field"
)

func Example() {
	// This example checks the freshman's dream identity for finite fields:
	//
	//   (x + y)^p = x^p + y^p
	//
	// where p is the characteristic of the field.
	// Create the Galois field GF(11^3), and pick any two elements x and y.
	p, n := 11, 3
	k := field.NewPrimeExtDeg(big.NewInt(int64(p)), n)
	x := k.NewZero().SetCoeffs(big.NewInt(10), big.NewInt(2), big.NewInt(1))
	y := k.NewZero().SetCoeffs(big.NewInt(1), big.NewInt(2), big.NewInt(6))

	// Compute (x + y)^p.
	xPlusY := k.NewZero().Add(x, y)
	xyp := k.NewOne()
	for range p {
		xyp.Mul(xyp, xPlusY)
	}

	// Compute x^p + y^p.
	xp, yp := k.NewOne(), k.NewOne()
	for range p {
		xp.Mul(xp, x)
		yp.Mul(yp, y)
	}
	xpPlusyp := k.NewZero().Add(xp, yp)

	fmt.Println("(x + y)^p == x^p + y^p")
	fmt.Println(xyp, "==", xpPlusyp)

	// Output:
	// (x + y)^p == x^p + y^p
	// [6 4 5] == [6 4 5]
}

func ExampleFactor() {
	variables := map[string]nag.Symbol{"x": 0}
	order, _ := new(big.Int).SetString("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10)
	k := field.NewPrimeExtDeg(order, 1)
	p, _ := field.Parse(variables, k, "x^2+2x-8")

	// Factor p into irreducible polynomials.
	factors := field.Factor(p)
	slices.SortFunc(factors, func(a, b field.IrrFactor[*field.PrimeExt]) int {
		return cmp.Compare(len(a.P.String()), len(b.P.String()))
	})

	fmt.Println(p)
	fmt.Println("==", factors)

	// Output:
	// x^2+2x+21888242871839275222246405745257275088548364400416034343698204186575808495609
	// == [{x+4 1} {x+21888242871839275222246405745257275088548364400416034343698204186575808495615 1}]
}
