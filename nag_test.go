package nag

import (
	"flag"
	"fmt"
	"log"
	"math/big"
	"slices"
	"testing"
)

func TestAddObstructions(t *testing.T) {
	tests := []struct {
		b   []obstruction
		g   []*Polynomial
		obs []obstruction
	}{
		// OBS(2), Example 5.12.
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
		// OBS(4), Example 5.12.
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
		// OBS(2), Example 5.12.
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
		// OBS(4), Example 5.12.
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
		// Example 3.13, M. Kreuzer, arXiv:1302.3805
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
		// Example 3.10, M. Kreuzer, arXiv:1302.3805
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
		obs      []obstruction
		ordering Ordering
		removed  []obstruction
	}{
		// OBS(2) Example 5.10, Mora.
		{
			obs: []obstruction{
				{i: 2, j: 2, iLeft: Monomial{2, 1}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{1, 2}},
				{i: 1, j: 2, iLeft: Monomial{}, iRight: Monomial{2}, jLeft: Monomial{1}, jRight: Monomial{}},
				{i: 1, j: 2, iLeft: Monomial{2}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{1}},
			},
			ordering: Deglex,
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
			ordering: Deglex,
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
			ordering: Deglex,
			removed: []obstruction{
				{i: 1, j: 5, iLeft: Monomial{1}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{1}},
				{i: 3, j: 5, iLeft: Monomial{1, 1}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2}},
				{i: 4, j: 5, iLeft: Monomial{}, iRight: Monomial{1, 2}, jLeft: Monomial{2}, jRight: Monomial{}},
			},
		},
		// Example 3.10, M. Kreuzer, arXiv:1302.3805
		{
			obs: []obstruction{
				{i: 1, j: 2, iLeft: Monomial{1, 2}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{}},
				{i: 1, j: 2, iLeft: Monomial{}, iRight: Monomial{1, 2}, jLeft: Monomial{}, jRight: Monomial{}},
				{i: 2, j: 2, iLeft: Monomial{}, iRight: Monomial{1, 2}, jLeft: Monomial{1, 2}, jRight: Monomial{}},
			},
			ordering: Deglex,
			removed: []obstruction{
				{i: 1, j: 2, iLeft: Monomial{}, iRight: Monomial{1, 2}, jLeft: Monomial{}, jRight: Monomial{}},
			},
		},
		{
			obs: []obstruction{
				{i: 1, j: 1, iLeft: Monomial{1, 2, 3}, iRight: Monomial{}, jLeft: Monomial{1, 2}, jRight: Monomial{2, 1}},
				{i: 1, j: 1, iLeft: Monomial{1, 2}, iRight: Monomial{}, jLeft: Monomial{1, 2}, jRight: Monomial{2, 1}},
			},
			ordering: Deglex,
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

			remove4b(obs, test.ordering)
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
		// Example 4.2.3.
		{
			o:      obstruction{i: 1, j: 2, iLeft: Monomial{1, 2, 1}, iRight: Monomial{}, jLeft: Monomial{1}, jRight: Monomial{1, 1, 2, 1, 2}},
			shrunk: obstruction{i: 1, j: 2, iLeft: Monomial{2, 1}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{1, 1, 2, 1, 2}},
		},
		// Example 3.5, M. Kreuzer, arXiv:1302.3805
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
		// Example 4.2.3.
		{
			obstruction: obstruction{i: 1, j: 3, iLeft: Monomial{1, 2, 1}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{1, 2}},
			im:          Monomial{1, 2, 1, 2},
			jm:          Monomial{1, 2, 1, 1, 2},
			overlap:     true,
		},
		// Example 4.2.3.
		{
			obstruction: obstruction{i: 1, j: 2, iLeft: Monomial{1, 2, 1}, iRight: Monomial{}, jLeft: Monomial{1}, jRight: Monomial{1, 1, 2, 1, 2}},
			im:          Monomial{1, 2, 1, 2},
			jm:          Monomial{2},
			overlap:     false,
		},
		// Example 4.2.5.
		{
			obstruction: obstruction{i: 2, j: 3, iLeft: Monomial{1, 2, 1, 1}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2, 2}},
			im:          Monomial{2, 2, 2},
			jm:          Monomial{1, 2, 1, 1, 2},
			overlap:     true,
		},
		// Example 4.2.5.
		{
			obstruction: obstruction{i: 1, j: 2, iLeft: Monomial{1, 2}, iRight: Monomial{2}, jLeft: Monomial{1, 2, 1, 1}, jRight: Monomial{}},
			im:          Monomial{1, 1, 2, 2},
			jm:          Monomial{2, 2, 2},
			overlap:     true,
		},
		// Case d, example 4.1.8.
		{
			obstruction: obstruction{i: 1, j: 2, iLeft: Monomial{}, iRight: Monomial{5, 6, 5, 6, 5, 6, 2}, jLeft: Monomial{4, 6, 5, 6, 5, 6, 5}, jRight: Monomial{}},
			im:          Monomial{4, 6, 5, 6, 5, 6, 5},
			jm:          Monomial{5, 6, 5, 6, 5, 6, 2},
			overlap:     false,
		},
		// Case e, example 4.1.8.
		{
			obstruction: obstruction{i: 1, j: 2, iLeft: Monomial{}, iRight: Monomial{6, 5, 6, 5, 6, 5, 6, 2}, jLeft: Monomial{4, 6, 5, 6, 5, 6, 5, 6}, jRight: Monomial{}},
			im:          Monomial{4, 6, 5, 6, 5, 6, 5},
			jm:          Monomial{5, 6, 5, 6, 5, 6, 2},
			overlap:     false,
		},
		// Case f, example 4.1.8.
		{
			obstruction: obstruction{i: 1, j: 2, iLeft: Monomial{}, iRight: Monomial{6, 5, 6, 5, 6, 5, 6, 5, 6, 2}, jLeft: Monomial{4, 6, 5, 6, 5, 6, 5, 6, 5, 6}, jRight: Monomial{}},
			im:          Monomial{4, 6, 5, 6, 5, 6, 5},
			jm:          Monomial{5, 6, 5, 6, 5, 6, 2},
			overlap:     false,
		},
		// Case a, Example 2.9, M. Kreuzer, arXiv:1302.3805
		{
			obstruction: obstruction{i: 1, j: 1, iLeft: Monomial{}, iRight: Monomial{1, 1, 1}, jLeft: Monomial{1, 1, 1}, jRight: Monomial{}},
			im:          Monomial{1, 1},
			jm:          Monomial{1, 1},
			overlap:     false,
		},
		// Case a, Example 2.9, M. Kreuzer, arXiv:1302.3805
		{
			obstruction: obstruction{i: 1, j: 1, iLeft: Monomial{}, iRight: Monomial{2, 1, 1}, jLeft: Monomial{1, 1, 2}, jRight: Monomial{}},
			im:          Monomial{1, 1},
			jm:          Monomial{1, 1},
			overlap:     false,
		},
		// Case b, Example 2.9, M. Kreuzer, arXiv:1302.3805
		{
			obstruction: obstruction{i: 1, j: 2, iLeft: Monomial{}, iRight: Monomial{1, 1, 2}, jLeft: Monomial{1, 1, 1}, jRight: Monomial{}},
			im:          Monomial{1, 1},
			jm:          Monomial{1, 2},
			overlap:     false,
		},
		// Case b, Example 2.9, M. Kreuzer, arXiv:1302.3805
		{
			obstruction: obstruction{i: 1, j: 2, iLeft: Monomial{1, 2, 2}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2, 1, 1}},
			im:          Monomial{1, 1},
			jm:          Monomial{1, 2},
			overlap:     false,
		},
		// Case c, Example 2.9, M. Kreuzer, arXiv:1302.3805
		{
			obstruction: obstruction{i: 2, j: 2, iLeft: Monomial{}, iRight: Monomial{1, 1, 2}, jLeft: Monomial{1, 2, 1}, jRight: Monomial{}},
			im:          Monomial{1, 2},
			jm:          Monomial{1, 2},
			overlap:     false,
		},
		// Example 3.10, M. Kreuzer, arXiv:1302.3805
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
		// Example 4.1.8.
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
		// Example 3.2.2.
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
		// Example 3.2.5.
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
		// NCAlgebra 6.4
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
		// NCAlgebra 6.5
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

			f := NewPolynomial(test.f.ordering)
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
			z := NewPolynomial(test.z.ordering)
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
			z := NewPolynomial(test.z.ordering).Mul(test.x, test.y)
			if z.Cmp(test.z) != 0 {
				t.Errorf("%v", z)
			}
		})
	}
}

func TestOrdering(t *testing.T) {
	tests := []struct {
		words    []Monomial
		ordering Ordering
		sorted   []Monomial
	}{
		{
			words:    []Monomial{{2}, {1, 1, 2}, {1, 2}},
			ordering: lexicographic,
			sorted:   []Monomial{{1, 1, 2}, {1, 2}, {2}},
		},
		{
			words:    []Monomial{{1, 3, 4}, {1, 2, 4}, {1, 3, 6}, {3, 4, 6}, {2, 5, 6}, {1, 5, 6}, {2, 3, 4}, {3, 4, 5}, {1, 4, 5}, {1, 4, 6}, {2, 3, 6}, {2, 4, 5}, {1, 3, 5}, {1, 2, 5}, {2, 3, 5}, {1, 2, 3}, {4, 5, 6}, {3, 5, 6}, {2, 4, 6}, {1, 2, 6}},
			ordering: lexicographic,
			sorted:   []Monomial{{1, 2, 3}, {1, 2, 4}, {1, 2, 5}, {1, 2, 6}, {1, 3, 4}, {1, 3, 5}, {1, 3, 6}, {1, 4, 5}, {1, 4, 6}, {1, 5, 6}, {2, 3, 4}, {2, 3, 5}, {2, 3, 6}, {2, 4, 5}, {2, 4, 6}, {2, 5, 6}, {3, 4, 5}, {3, 4, 6}, {3, 5, 6}, {4, 5, 6}},
		},
		{
			words:    []Monomial{{2, 3}, {1, 3}, {1, 1}, {3, 3}, {1, 2}, {2, 2}},
			ordering: lexicographic,
			sorted:   []Monomial{{1, 1}, {1, 2}, {1, 3}, {2, 2}, {2, 3}, {3, 3}},
		},
		{
			words:    []Monomial{{'b', 'e', 'n', 'e', 'f', 'i', 't'}, {'b', 'e'}, {'b', 'e', 'n', 't'}, {'b', 'a', 'r', 'n', 'a', 'c', 'l', 'e'}, {'b', 'e', 'e', 'n'}},
			ordering: lexicographic,
			sorted:   []Monomial{{'b', 'a', 'r', 'n', 'a', 'c', 'l', 'e'}, {'b', 'e'}, {'b', 'e', 'e', 'n'}, {'b', 'e', 'n', 'e', 'f', 'i', 't'}, {'b', 'e', 'n', 't'}},
		},
		{
			words:    []Monomial{{1, 1}, {1, 1, 2}, {2, 1, 1}, {2, 2}, {1, 2, 2}, {2, 2, 2}, {}, {1, 2}, {1, 1, 1}, {2, 1, 2}, {2}, {1}, {2, 2, 1}, {1, 2, 1}, {2, 1}},
			ordering: Deglex,
			sorted:   []Monomial{{}, {1}, {2}, {1, 1}, {1, 2}, {2, 1}, {2, 2}, {1, 1, 1}, {1, 1, 2}, {1, 2, 1}, {1, 2, 2}, {2, 1, 1}, {2, 1, 2}, {2, 2, 1}, {2, 2, 2}},
		},
	}
	for i, test := range tests {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			t.Parallel()
			slices.SortFunc(test.words, test.ordering)
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
		z = NewPolynomial(z.ordering).Mul(z, y[i])
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
