package nag

// GaussElim performs [Gaussian elimination] on the matrix m, reducing it to
// reduced row echelon form. The matrix m is modified in place and returned.
//
// [Gaussian elimination]: https://en.wikipedia.org/wiki/Gaussian_elimination
func GaussElim[K Field[K]](m [][]K) [][]K {
	rows, cols := len(m), len(m[0])
	k := m[0][0]
	zero := k.NewZero()
	tmp0, tmp1 := k.NewZero(), k.NewZero()

	lead := 0
	for r := range rows {
		if lead >= cols {
			break
		}

		i := r
		for m[i][lead].Equal(zero) {
			i++
			if i == rows {
				i = r
				lead++
				if lead == cols {
					return m
				}
			}
		}
		m[i], m[r] = m[r], m[i]

		pivot := tmp0.Set(m[r][lead])
		for j := range cols {
			m[r][j].Div(m[r][j], pivot)
		}

		for k := range rows {
			if k == r {
				continue
			}
			factor := tmp0.Set(m[k][lead])
			if factor.Equal(zero) {
				continue
			}
			for j := range cols {
				m[k][j].Sub(m[k][j], tmp1.Mul(factor, m[r][j]))
			}
		}

		lead++
	}
	return m
}
