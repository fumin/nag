package field_test

import (
	"fmt"

	"github.com/fumin/nag/field"
)

func Example() {
	order := 1759
	a := field.NewPrime(order, 550)
	b := field.NewPrime(order, 0)
	b.Inv(a)
	fmt.Printf("1/%v = %v\n", a, b)

	c := field.NewPrime(order, 355)
	d := field.NewPrime(order, 0)
	d.Mul(a, c)
	fmt.Printf("%v * %v = %v\n", a, c, d)

	// Output:
	// 1/550 = 355
	// 550 * 355 = 1
}
