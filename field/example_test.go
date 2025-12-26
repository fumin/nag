package field_test

import (
	"fmt"

	"github.com/fumin/nag/field"
)

func Example() {
	// This example checks the freshman's dream identity for finite fields:
	//
	//   (x + y)^p = x^p + y^p
	//
	// where p is the characteristic of the field.
	// Create the Galois field GF(101^2), and pick any two elements x and y.
	p, n := 101, 2
	irr := field.NewIrreduciblePoly(p, n)
	x, y := irr.Ext(6158), irr.Ext(8033)

	// Compute (x + y)^p.
	xPlusY := irr.Ext(0).Add(x, y)
	xyp := irr.Ext(1)
	for range p {
		xyp.Mul(xyp, xPlusY)
	}

	// Compute x^p + y^p.
	xp, yp := irr.Ext(1), irr.Ext(1)
	for range p {
		xp.Mul(xp, x)
		yp.Mul(yp, y)
	}
	xpPlusyp := irr.Ext(0).Add(xp, yp)

	fmt.Println("(x + y)^p == x^p + y^p")
	fmt.Println(xyp, "==", xpPlusyp)

	// Output:
	// (x + y)^p == x^p + y^p
	// 6414 == 6414
}
