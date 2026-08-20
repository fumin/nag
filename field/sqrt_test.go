package field

import (
	"fmt"
	"math/big"
	"slices"
	"testing"
)

func TestSqrt(t *testing.T) {
	tests := []struct {
		p     int64
		n     int
		sqrts map[int64][]int64
	}{
		{
			p: 7, n: 1,
			sqrts: map[int64][]int64{
				1: {1, 6},
				2: {4, 3},
				4: {2, 5},
			},
		},
		{
			p: 11, n: 1,
			sqrts: map[int64][]int64{
				1: {1, 10},
				3: {5, 6},
				4: {9, 2},
				5: {4, 7},
				9: {3, 8},
			},
		},
		{
			p: 3, n: 3,
			sqrts: map[int64][]int64{
				1:  {1, 2},
				6:  {22, 17},
				7:  {20, 10},
				8:  {25, 14},
				9:  {6, 3},
				11: {13, 26},
				12: {16, 23},
				13: {7, 5},
				15: {9, 18},
				16: {8, 4},
				20: {15, 21},
				22: {12, 24},
				25: {11, 19},
			},
		},
		{
			p: 41, n: 1,
			sqrts: map[int64][]int64{
				1:  {1, 40},
				2:  {17, 24},
				4:  {2, 39},
				5:  {28, 13},
				8:  {34, 7},
				9:  {38, 3},
				10: {16, 25},
				16: {37, 4},
				18: {10, 31},
				20: {26, 15},
				21: {12, 29},
				23: {33, 8},
				25: {36, 5},
				31: {20, 21},
				32: {14, 27},
				33: {19, 22},
				36: {6, 35},
				37: {18, 23},
				39: {11, 30},
				40: {32, 9},
			},
		},
	}
	for i, test := range tests {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			k := NewPrimeExtDeg(big.NewInt(test.p), test.n)
			order := new(big.Int).Exp(k.Characteristic(), big.NewInt(int64(k.Degree())), nil)
			for i := big.NewInt(1); i.Cmp(order) < 0; i.Add(i, big.NewInt(1)) {
				want := make([]*PrimeExt, 0)
				for _, ith := range test.sqrts[i.Int64()] {
					want = append(want, setIth(k.NewZero(), big.NewInt(ith)))
				}

				ki := setIth(k.NewZero(), i)
				var got []*PrimeExt
				s, ok := Sqrt(ki)
				if ok {
					got = append(got, s)
					ms := k.NewZero()
					got = append(got, ms.Sub(ms, s))
				}

				if !slices.EqualFunc(got, want, func(a, b *PrimeExt) bool { return a.Equal(b) }) {
					t.Errorf("GF(%d, %d) Sqrt(%d) = %v want %v", test.p, test.n, i, got, want)
				}
				for _, s := range got {
					if s2 := k.Mul(s, s); !s2.Equal(ki) {
						t.Errorf("GF(%d, %d) %v^2 = %v want %v", test.p, test.n, s, s2, ki)
					}
				}
			}
		})
	}
}

func TestSqrtLoop(t *testing.T) {
	tests := []struct {
		p       int64
		n       int
		numSqrt int
	}{
		{p: 41, n: 2, numSqrt: 841},
	}
	for i, test := range tests {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			k := NewPrimeExtDeg(big.NewInt(test.p), test.n)
			order := new(big.Int).Exp(k.Characteristic(), big.NewInt(int64(k.Degree())), nil)
			ith := func(x *PrimeExt) *big.Int {
				xi := big.NewInt(0)
				for d, c := range x.Coeffs() {
					ccd := new(big.Int).Exp(k.Characteristic(), big.NewInt(int64(d)), nil)
					ccd.Mul(ccd, c)
					xi.Add(xi, ccd)
				}
				return xi
			}

			sqrts := make(map[int64][]int64)
			for i := big.NewInt(0); i.Cmp(order) < 0; i.Add(i, big.NewInt(1)) {
				ki := setIth(k.NewZero(), i)
				ki2 := k.NewZero().Mul(ki, ki)
				i2 := ith(ki2).Int64()
				sqrts[i2] = append(sqrts[i2], i.Int64())
			}

			numSqrt := 1 // for i = 0
			for i := big.NewInt(1); i.Cmp(order) < 0; i.Add(i, big.NewInt(1)) {
				ki := setIth(k.NewZero(), i)
				var got []int64
				s, ok := Sqrt(ki)
				if ok {
					got = append(got, ith(s).Int64())
					ms := k.NewZero()
					got = append(got, ith(ms.Sub(ms, s)).Int64())
					numSqrt++
				}
				slices.Sort(got)

				want := sqrts[i.Int64()]
				if !slices.Equal(got, want) {
					t.Errorf("GF(%d, %d) Sqrt(%v) = %v want %v", test.p, test.n, i, got, want)
				}
			}
			if numSqrt != test.numSqrt {
				t.Errorf("GF(%d, %d) %d != %d", test.p, test.n, numSqrt, test.numSqrt)
			}
		})
	}
}
