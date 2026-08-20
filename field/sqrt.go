package field

import (
	"math/big"

	"github.com/fumin/nag"
)

// Sqrt returns the square root of x.
func Sqrt[K Finite[K]](x K) (K, bool) {
	k := x.NewZero()
	if x.Equal(k) {
		return x, true
	}
	q := new(big.Int).Exp(k.Characteristic(), big.NewInt(int64(k.Degree())), nil)

	// Characteristic 2.
	two := big.NewInt(2)
	if k.Characteristic().Cmp(two) == 0 {
		return nag.Pow(k.Set(x), q.Div(q, two)), true
	}

	// Euler's criterion: a is a nonzero square iff a^((q-1)/2) == 1.
	one := big.NewInt(1)
	qm1 := q.Sub(q, one)
	half := new(big.Int).Div(qm1, two)
	tmp := new(big.Int)
	kone := k.NewOne()
	if euler := nag.Pow(k.Set(x), tmp.Set(half)); !euler.Equal(kone) {
		return x, false
	}
	z := findNonResidue(k, half, kone, one)

	// Factor q-1 = Q * 2^S with Q odd.
	Q := qm1
	S := half.SetInt64(0)
	for tmp.Mod(Q, two).Sign() == 0 {
		Q.Div(Q, two)
		S.Add(S, one)
	}

	c := nag.Pow(k.NewZero().Set(z), tmp.Set(Q))
	t := nag.Pow(k.NewZero().Set(x), tmp.Set(Q))
	Q.Add(Q, one)
	Q.Div(Q, two)
	R := nag.Pow(k.NewZero().Set(x), Q)
	M := S

	i, j, bound := new(big.Int), new(big.Int), new(big.Int)
	b := k.NewZero()
	for {
		if t.Equal(kone) {
			return R, true
		}

		i.SetInt64(0)
		tmp := k.Set(t)
		for !tmp.Equal(kone) {
			tmp.Mul(tmp, tmp)
			i.Add(i, one)
			if i.Cmp(M) == 0 {
				panic("euler criterion violated")
			}
		}

		b.Set(c)
		bound.Sub(M, i)
		bound.Sub(bound, one)
		for j.SetInt64(0); j.Cmp(bound) < 0; j.Add(j, one) {
			b.Mul(b, b)
		}

		M.Set(i)
		c.Mul(b, b)
		t.Mul(t, c)
		R.Mul(R, b)
	}
}

func findNonResidue[K Finite[K]](k K, half *big.Int, kone K, one *big.Int) K {
	order := new(big.Int).Exp(k.Characteristic(), big.NewInt(int64(k.Degree())), nil)
	tmp := new(big.Int)
	z := k.NewZero()
	for i := big.NewInt(2); i.Cmp(order) < 0; i.Add(i, one) {
		setIth(z, i)
		if isNonResidue(z, tmp.Set(half), kone) {
			return setIth(z, i)
		}
	}
	panic("unable to find non-residue")
}

func isNonResidue[K Finite[K]](a K, half *big.Int, kone K) bool {
	euler := nag.Pow(a, half)
	return !euler.Equal(kone)
}
