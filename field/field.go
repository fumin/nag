// Package field implements [finite field] arithmetic.
//
// [finite field]: https://en.wikipedia.org/wiki/Finite_field
package field

import (
	"math/big"

	"github.com/fumin/nag"
)

// A Finite is an element of a finite field with order q = p^n where p is a prime number.
type Finite[K nag.Field[K]] interface {
	nag.Field[K]
	// Characteristic returns p in the finite field GF(p^n).
	Characteristic() int
	// PrimePower returns n in the finite field GF(p^n).
	PrimePower() int
}

// A Prime is an element of a finite field whose order is a prime number.
type Prime struct {
	// Order is the number of elements in the finite field which must be a prime number.
	Order *big.Int
	// I is the integer representation of the element.
	I *big.Int
}

// NewPrime returns the i'th element in the finite field GF(p) where p is a prime number.
func NewPrime(p, i int) *Prime {
	return &Prime{Order: big.NewInt(int64(p)), I: big.NewInt(int64(i))}
}

// Equal reports whether x and y are equal.
func (x *Prime) Equal(y *Prime) bool {
	return (x.Order.Cmp(y.Order) == 0) && (x.I.Cmp(y.I) == 0)
}

// Add sets z to the sum x+y and returns z.
func (z *Prime) Add(x, y *Prime) *Prime {
	z.Order.Set(x.Order)
	z.I.Add(x.I, y.I)
	z.I.Mod(z.I, z.Order)
	return z
}

// Sub sets z to the difference x-y and returns z.
func (z *Prime) Sub(x, y *Prime) *Prime {
	z.Order.Set(x.Order)
	z.I.Sub(x.I, y.I)
	z.I.Mod(z.I, z.Order)
	return z
}

// Mul sets z to the product x*y and returns z.
func (z *Prime) Mul(x, y *Prime) *Prime {
	z.Order.Set(x.Order)
	z.I.Mul(x.I, y.I)
	z.I.Mod(z.I, z.Order)
	return z
}

// Div sets z to the quotient x/y and returns z.
func (z *Prime) Div(x, y *Prime) *Prime {
	z.Inv(y)
	z.I.Mul(x.I, z.I)
	z.I.Mod(z.I, z.Order)
	return z
}

// Inv sets z to 1/x and returns z.
func (z *Prime) Inv(x *Prime) *Prime {
	if x.I.Sign() == 0 {
		panic("division by zero")
	}
	z.Order.Set(x.Order)
	z.I.ModInverse(x.I, z.Order)
	return z
}

// String returns the integer representation of x as a string.
func (x *Prime) String() string {
	return x.I.String()
}

// Characteristic returns p in the finite field GF(p^n) where x belongs.
func (x *Prime) Characteristic() int {
	return int(x.Order.Int64())
}

// Characteristic returns n in the finite field GF(p^n) where x belongs.
func (x *Prime) PrimePower() int { return 1 }

// A Primes is a finite field whose order is a prime number.
type Primes struct {
	// Order is the number of elements in the field which must be a prime number.
	Order *big.Int
}

// NewPrime returns the the finite field GF(p) where p is a prime number.
func NewPrimes(p int) *Primes {
	return &Primes{Order: big.NewInt(int64(p))}
}

// NewZero returns the additive identity 0.
func (x Primes) NewZero() *Prime {
	return &Prime{Order: big.NewInt(0).Set(x.Order), I: big.NewInt(0)}
}

// NewOne returns the multiplicative identity 1.
func (x Primes) NewOne() *Prime {
	return &Prime{Order: big.NewInt(0).Set(x.Order), I: big.NewInt(1)}
}
