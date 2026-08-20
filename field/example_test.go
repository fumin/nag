package field_test

import (
	"cmp"
	"fmt"
	"math/big"
	"slices"

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
	irr := field.NewIrreduciblePoly(big.NewInt(int64(p)), n)
	x, y := irr.Ext(big.NewInt(153)), irr.Ext(big.NewInt(749))

	// Compute (x + y)^p.
	xPlusY := irr.Ext(big.NewInt(0)).Add(x, y)
	xyp := irr.Ext(big.NewInt(1))
	for range p {
		xyp.Mul(xyp, xPlusY)
	}

	// Compute x^p + y^p.
	xp, yp := irr.Ext(big.NewInt(1)), irr.Ext(big.NewInt(1))
	for range p {
		xp.Mul(xp, x)
		yp.Mul(yp, y)
	}
	xpPlusyp := irr.Ext(big.NewInt(0)).Add(xp, yp)

	fmt.Println("(x + y)^p == x^p + y^p")
	fmt.Println(xyp, "==", xpPlusyp)

	// Output:
	// (x + y)^p == x^p + y^p
	// 655 == 655
}

func ExampleFactor() {
	order, _ := new(big.Int).SetString("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10)
	k := field.NewIrreduciblePoly(order, 1).Ext(big.NewInt(0))
	p, _ := field.Parse(k, "x^2+2x-8")

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
