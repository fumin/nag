package nag

import (
	"fmt"
	"testing"
)

func TestGaussElim(t *testing.T) {
	tests := []struct {
		m    [][]*Rat
		want [][]*Rat
	}{
		// https://en.wikipedia.org/wiki/Gaussian_elimination#Example_of_the_algorithm
		{
			m: [][]*Rat{
				{NewRat(2, 1), NewRat(1, 1), NewRat(-1, 1), NewRat(8, 1)},
				{NewRat(-3, 1), NewRat(-1, 1), NewRat(2, 1), NewRat(-11, 1)},
				{NewRat(-2, 1), NewRat(1, 1), NewRat(2, 1), NewRat(-3, 1)},
			},
			want: [][]*Rat{
				{NewRat(1, 1), NewRat(0, 1), NewRat(0, 1), NewRat(2, 1)},
				{NewRat(0, 1), NewRat(1, 1), NewRat(0, 1), NewRat(3, 1)},
				{NewRat(0, 1), NewRat(0, 1), NewRat(1, 1), NewRat(-1, 1)},
			},
		},
		{
			// https://en.wikipedia.org/wiki/Gaussian_elimination#Finding_the_inverse_of_a_matrix
			m: [][]*Rat{
				{NewRat(2, 1), NewRat(-1, 1), NewRat(0, 1), NewRat(1, 1), NewRat(0, 1), NewRat(0, 1)},
				{NewRat(-1, 1), NewRat(2, 1), NewRat(-1, 1), NewRat(0, 1), NewRat(1, 1), NewRat(0, 1)},
				{NewRat(0, 1), NewRat(-1, 1), NewRat(2, 1), NewRat(0, 1), NewRat(0, 1), NewRat(1, 1)},
			},
			want: [][]*Rat{
				{NewRat(1, 1), NewRat(0, 1), NewRat(0, 1), NewRat(3, 4), NewRat(1, 2), NewRat(1, 4)},
				{NewRat(0, 1), NewRat(1, 1), NewRat(0, 1), NewRat(1, 2), NewRat(1, 1), NewRat(1, 2)},
				{NewRat(0, 1), NewRat(0, 1), NewRat(1, 1), NewRat(1, 4), NewRat(1, 2), NewRat(3, 4)},
			},
		},
		{
			// https://www.cs.bu.edu/fac/snyder/cs132-book/L03RowReductions.html#example
			m: [][]*Rat{
				{NewRat(0, 1), NewRat(3, 1), NewRat(-6, 1), NewRat(6, 1), NewRat(4, 1), NewRat(-5, 1)},
				{NewRat(3, 1), NewRat(-7, 1), NewRat(8, 1), NewRat(-5, 1), NewRat(8, 1), NewRat(9, 1)},
				{NewRat(3, 1), NewRat(-9, 1), NewRat(12, 1), NewRat(-9, 1), NewRat(6, 1), NewRat(15, 1)},
			},
			want: [][]*Rat{
				{NewRat(1, 1), NewRat(0, 1), NewRat(-2, 1), NewRat(3, 1), NewRat(0, 1), NewRat(-24, 1)},
				{NewRat(0, 1), NewRat(1, 1), NewRat(-2, 1), NewRat(2, 1), NewRat(0, 1), NewRat(-7, 1)},
				{NewRat(0, 1), NewRat(0, 1), NewRat(0, 1), NewRat(0, 1), NewRat(1, 1), NewRat(4, 1)},
			},
		},
		{
			// System with no solutions.
			m: [][]*Rat{
				{NewRat(1, 1), NewRat(1, 1), NewRat(1, 1), NewRat(2, 1)},
				{NewRat(1, 1), NewRat(2, 1), NewRat(3, 1), NewRat(7, 1)},
				{NewRat(2, 1), NewRat(3, 1), NewRat(4, 1), NewRat(13, 1)},
			},
			want: [][]*Rat{
				{NewRat(1, 1), NewRat(0, 1), NewRat(-1, 1), NewRat(0, 1)},
				{NewRat(0, 1), NewRat(1, 1), NewRat(2, 1), NewRat(0, 1)},
				{NewRat(0, 1), NewRat(0, 1), NewRat(0, 1), NewRat(1, 1)},
			},
		},
	}
	for i, test := range tests {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			got := GaussElim(cloneM(test.m))
			if !mEq(got, test.want) {
				t.Errorf("GaussElim(%v) = %v want %v", test.m, got, test.want)
			}
		})
	}
}

func cloneM[K Field[K]](m [][]K) [][]K {
	k := m[0][0]
	c := make([][]K, len(m))
	for i := range m {
		c[i] = make([]K, len(m[i]))
		for j := range m[i] {
			c[i][j] = k.NewZero().Set(m[i][j])
		}
	}
	return c
}

func mEq[T Magma[T]](a, b [][]T) bool {
	for i := range a {
		for j := range a[i] {
			if !a[i][j].Equal(b[i][j]) {
				return false
			}
		}
	}
	return true
}
