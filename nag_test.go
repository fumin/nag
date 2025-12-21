// Tests in this file come from the following references:
//
// Clemens Hofstadler. Certifying operator identities and ideal membership of noncommutative polynomials. Master’s thesis. Johannes Kepler University Linz, Austria, 2020.
// Green, Edward L. "Noncommutative Gröbner bases, and projective resolutions." Computational Methods for Representations of Groups and Algebras: Euroconference in Essen (Germany), April 1–5, 1977. Basel: Birkhäuser Basel, 1999. 29-60.
// Kreuzer, Martin, and Xingqiang Xiu. "Non-commutative Gebauer-Moeller criteria." arXiv preprint arXiv:1302.3805 (2013).
// Mora, Teo. "An introduction to commutative and noncommutative Gröbner bases." Theoretical Computer Science 134.1 (1994): 131-173.
// Xiu, Xingqiang. "Non-commutative Gröbner bases and applications." PhD diss., Universität Passau, 2012.
// Bergman software by Jörgen Backelin, https://servus.math.su.se/bergman
// NCAlgebra Mathematica library, version 6.0.3, https://mathweb.ucsd.edu/~ncalg
package nag

import (
	"flag"
	"fmt"
	"log"
	"math/big"
	"slices"
	"testing"
)

func TestBuchbergerHomogeneous(t *testing.T) {
	tests := []struct {
		ideal    []*Polynomial
		maxDeg   int
		basis    []*Polynomial
		complete bool
	}{
		// Section 1.2.2 Commutative algebras, Bergman manual.
		{
			ideal: []*Polynomial{
				parseMust(map[string]Symbol{"x": 2, "y": 1}, Deglex, "x^2-y^2"),
				parseMust(map[string]Symbol{"x": 2, "y": 1}, Deglex, "xy"),
				parseMust(map[string]Symbol{"x": 2, "y": 1}, Deglex, "xy-yx"),
			},
			maxDeg: 4,
			basis: []*Polynomial{
				parseMust(map[string]Symbol{"x": 2, "y": 1}, Deglex, "yx"),
				parseMust(map[string]Symbol{"x": 2, "y": 1}, Deglex, "xy"),
				parseMust(map[string]Symbol{"x": 2, "y": 1}, Deglex, "x^2-y^2"),
				parseMust(map[string]Symbol{"x": 2, "y": 1}, Deglex, "y^3"),
			},
			complete: true,
		},
		// Section 1.2.3 Non-commutative algebras, Bergman manual.
		{
			ideal: []*Polynomial{
				parseMust(map[string]Symbol{"x": 1, "y": 2}, Deglex, "x^2-y^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2}, Deglex, "xy"),
			},
			maxDeg: 4,
			basis: []*Polynomial{
				parseMust(map[string]Symbol{"x": 1, "y": 2}, Deglex, "xy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2}, Deglex, "y^2-x^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2}, Deglex, "x^3"),
				parseMust(map[string]Symbol{"x": 1, "y": 2}, Deglex, "yx^2"),
			},
			complete: true,
		},
		// Section 2.4.3 Using eliminating ordering, Bergman manual.
		{
			ideal: []*Polynomial{
				parseMust(map[string]Symbol{"r": 4, "l": 3, "q": 2, "h": 1}, ElimOrder(), "r^2-q^2+lr-rh"),
				parseMust(map[string]Symbol{"r": 4, "l": 3, "q": 2, "h": 1}, ElimOrder(), "rq-qr+lq-qh"),
				parseMust(map[string]Symbol{"r": 4, "l": 3, "q": 2, "h": 1}, ElimOrder(), "lr-rl"),
				parseMust(map[string]Symbol{"r": 4, "l": 3, "q": 2, "h": 1}, ElimOrder(), "lq+ql-2qh"),
				parseMust(map[string]Symbol{"r": 4, "l": 3, "q": 2, "h": 1}, ElimOrder(), "rr-rl-2lr+2l^2-rh+lh"),
				parseMust(map[string]Symbol{"r": 4, "l": 3, "q": 2, "h": 1}, ElimOrder(), "rh-hr"),
				parseMust(map[string]Symbol{"r": 4, "l": 3, "q": 2, "h": 1}, ElimOrder(), "qh-hq"),
				parseMust(map[string]Symbol{"r": 4, "l": 3, "q": 2, "h": 1}, ElimOrder(), "lh-hl"),
			},
			maxDeg: 10,
			basis: []*Polynomial{
				parseMust(map[string]Symbol{"r": 4, "l": 3, "q": 2, "h": 1}, ElimOrder(), "qh-hq"),
				parseMust(map[string]Symbol{"r": 4, "l": 3, "q": 2, "h": 1}, ElimOrder(), "q^7-13h^2q^5+39h^4q^3-27h^6q"),
				parseMust(map[string]Symbol{"r": 4, "l": 3, "q": 2, "h": 1}, ElimOrder(), "lh-hl"),
				parseMust(map[string]Symbol{"r": 4, "l": 3, "q": 2, "h": 1}, ElimOrder(), "lq+ql-2hq"),
				parseMust(map[string]Symbol{"r": 4, "l": 3, "q": 2, "h": 1}, ElimOrder(), "q^3l-3h^2ql-hq^3+3h^3q"),
				parseMust(map[string]Symbol{"r": 4, "l": 3, "q": 2, "h": 1}, ElimOrder(), "h^2ql^2-2h^3ql+1/48q^5-5/24h^2q^3+19/16h^4q"),
				parseMust(map[string]Symbol{"r": 4, "l": 3, "q": 2, "h": 1}, ElimOrder(), "h^2q^2l^2-2h^3q^2l+1/48q^6-5/24h^2q^4+19/16h^4q^2"),
				parseMust(map[string]Symbol{"r": 4, "l": 3, "q": 2, "h": 1}, ElimOrder(), "hl^3+1/3q^2l^2-5/3hq^2l-1/4h^3l+1/12q^4+1/2h^2q^2"),
				parseMust(map[string]Symbol{"r": 4, "l": 3, "q": 2, "h": 1}, ElimOrder(), "ql^3-3hql^2+11/4h^2ql-3/4h^3q"),
				parseMust(map[string]Symbol{"r": 4, "l": 3, "q": 2, "h": 1}, ElimOrder(), "l^4-2/3q^2l^2-1/4h^2l^2-1/6hq^2l+1/12q^4"),
				parseMust(map[string]Symbol{"r": 4, "l": 3, "q": 2, "h": 1}, ElimOrder(), "rh-hr"),
				parseMust(map[string]Symbol{"r": 4, "l": 3, "q": 2, "h": 1}, ElimOrder(), "rq-qr-ql+hq"),
				parseMust(map[string]Symbol{"r": 4, "l": 3, "q": 2, "h": 1}, ElimOrder(), "hqr-ql^2+5/2hql-1/4q^3-9/4h^2q"),
				parseMust(map[string]Symbol{"r": 4, "l": 3, "q": 2, "h": 1}, ElimOrder(), "q^2r+3l^3-5/2q^2l-3/4h^2l-3/4hq^2"),
				parseMust(map[string]Symbol{"r": 4, "l": 3, "q": 2, "h": 1}, ElimOrder(), "lr-1/2l^2-1/4hl-1/4q^2"),
				parseMust(map[string]Symbol{"r": 4, "l": 3, "q": 2, "h": 1}, ElimOrder(), "rl-1/2l^2-1/4hl-1/4q^2"),
				parseMust(map[string]Symbol{"r": 4, "l": 3, "q": 2, "h": 1}, ElimOrder(), "r^2-hr+1/2l^2+1/4hl-3/4q^2"),
			},
			complete: true,
		},
		// Section 2.5 Homogenisation, Bergman manual.
		{
			ideal: []*Polynomial{
				parseMust(map[string]Symbol{"x": 5, "y": 4, "z": 3, "t": 2, "f": 1}, ElimOrder(), "x^2-2f^2"),
				parseMust(map[string]Symbol{"x": 5, "y": 4, "z": 3, "t": 2, "f": 1}, ElimOrder(), "y^2-3f^2"),
				parseMust(map[string]Symbol{"x": 5, "y": 4, "z": 3, "t": 2, "f": 1}, ElimOrder(), "z^2-5f^2"),
				parseMust(map[string]Symbol{"x": 5, "y": 4, "z": 3, "t": 2, "f": 1}, ElimOrder(), "t-x-y-z"),
				// Commutative relations.
				parseMust(map[string]Symbol{"x": 5, "y": 4, "z": 3, "t": 2, "f": 1}, ElimOrder(), "yx-xy"),
				parseMust(map[string]Symbol{"x": 5, "y": 4, "z": 3, "t": 2, "f": 1}, ElimOrder(), "zx-xz"),
				parseMust(map[string]Symbol{"x": 5, "y": 4, "z": 3, "t": 2, "f": 1}, ElimOrder(), "zy-yz"),
				parseMust(map[string]Symbol{"x": 5, "y": 4, "z": 3, "t": 2, "f": 1}, ElimOrder(), "tx-xt"),
				parseMust(map[string]Symbol{"x": 5, "y": 4, "z": 3, "t": 2, "f": 1}, ElimOrder(), "ty-yt"),
				parseMust(map[string]Symbol{"x": 5, "y": 4, "z": 3, "t": 2, "f": 1}, ElimOrder(), "tz-zt"),
				parseMust(map[string]Symbol{"x": 5, "y": 4, "z": 3, "t": 2, "f": 1}, ElimOrder(), "fx-xf"),
				parseMust(map[string]Symbol{"x": 5, "y": 4, "z": 3, "t": 2, "f": 1}, ElimOrder(), "fy-yf"),
				parseMust(map[string]Symbol{"x": 5, "y": 4, "z": 3, "t": 2, "f": 1}, ElimOrder(), "fz-zf"),
				parseMust(map[string]Symbol{"x": 5, "y": 4, "z": 3, "t": 2, "f": 1}, ElimOrder(), "ft-tf"),
			},
			maxDeg: 11,
			basis: []*Polynomial{
				parseMust(map[string]Symbol{"x": 5, "y": 4, "z": 3, "t": 2, "f": 1}, ElimOrder(), "tf-ft"),
				parseMust(map[string]Symbol{"x": 5, "y": 4, "z": 3, "t": 2, "f": 1}, ElimOrder(), "t^8-40f^2t^6+352f^4t^4-960f^6t^2+576f^8"),
				parseMust(map[string]Symbol{"x": 5, "y": 4, "z": 3, "t": 2, "f": 1}, ElimOrder(), "zf-fz"),
				parseMust(map[string]Symbol{"x": 5, "y": 4, "z": 3, "t": 2, "f": 1}, ElimOrder(), "f^6z-5/576t^7+194/576f^2t^5-1520/576f^4t^3+2544/576f^6t"),
				parseMust(map[string]Symbol{"x": 5, "y": 4, "z": 3, "t": 2, "f": 1}, ElimOrder(), "zt-tz"),
				parseMust(map[string]Symbol{"x": 5, "y": 4, "z": 3, "t": 2, "f": 1}, ElimOrder(), "f^4tz-1/96t^6+40/96f^2t^4-376/96f^4t^2+480/96f^6"),
				parseMust(map[string]Symbol{"x": 5, "y": 4, "z": 3, "t": 2, "f": 1}, ElimOrder(), "f^2t^2z-96/80f^4z+1/80t^5-60/80f^2t^3-24/80f^4t"),
				parseMust(map[string]Symbol{"x": 5, "y": 4, "z": 3, "t": 2, "f": 1}, ElimOrder(), "t^3z-1/4t^4-20/4f^2t^2+24/4f^4"),
				parseMust(map[string]Symbol{"x": 5, "y": 4, "z": 3, "t": 2, "f": 1}, ElimOrder(), "z^2-5f^2"),
				parseMust(map[string]Symbol{"x": 5, "y": 4, "z": 3, "t": 2, "f": 1}, ElimOrder(), "yf-fy"),
				parseMust(map[string]Symbol{"x": 5, "y": 4, "z": 3, "t": 2, "f": 1}, ElimOrder(), "f^2y-3/2t^2z+6/2f^2z+1/2t^3+4/2f^2t"),
				parseMust(map[string]Symbol{"x": 5, "y": 4, "z": 3, "t": 2, "f": 1}, ElimOrder(), "yt-ty"),
				parseMust(map[string]Symbol{"x": 5, "y": 4, "z": 3, "t": 2, "f": 1}, ElimOrder(), "f^2ty+3f^2tz+1/8t^4-11/2f^2t^2+9f^4"),
				parseMust(map[string]Symbol{"x": 5, "y": 4, "z": 3, "t": 2, "f": 1}, ElimOrder(), "t^2y-7t^2z+12f^2z+2t^3+12f^2t"),
				parseMust(map[string]Symbol{"x": 5, "y": 4, "z": 3, "t": 2, "f": 1}, ElimOrder(), "zy-ty-tz+1/2t^2+3f^2"),
				parseMust(map[string]Symbol{"x": 5, "y": 4, "z": 3, "t": 2, "f": 1}, ElimOrder(), "yz-ty-tz+1/2t^2+3f^2"),
				parseMust(map[string]Symbol{"x": 5, "y": 4, "z": 3, "t": 2, "f": 1}, ElimOrder(), "y^2-3f^2"),
				parseMust(map[string]Symbol{"x": 5, "y": 4, "z": 3, "t": 2, "f": 1}, ElimOrder(), "x+y+z-t"),
			},
			complete: true,
		},
		// braid3-9, Example 4.2, Kreuzer, expected basis from bergman.
		{
			ideal: []*Polynomial{
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxy - zyz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyx - zxy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zxz - yzx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "x^3 + y^3 + z^3 + xyz"),
			},
			maxDeg: 9,
			basis: []*Polynomial{
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zxy-xyx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zxz-yzx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zyz-yxy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "z^3+y^3+xyz+x^3"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zy^3+zx^3-y^3z-xyz^2+xyxz-x^3z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyxz-xyx^2y"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^3xy+xyxyx+xyx^3+x^4y"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^3xz+yxyxy+xy^2zx+x^4z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yzx^2z-xyxzx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^2yx-yzx^2y"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^2yz+zx^4+y^3zx+xyxy^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zyxyx-yxyxy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zyxyz+zyx^3-zx^3y+y^3zy+yxyz^2+xyz^2y-xyxzy+x^3zy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zy^2xy-yxy^2z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zy^2zx-xyx^2y"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "z^2yxy+y^4z+xy^2xy+x^3yz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyx^2z-xyx^2zx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyxyz-xyx^3y"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxy^2zx-xyx^2y^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^4zx+yxyxy^2+xy^2xyx+x^3yzx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^3zx-y^3z^2x+xy^4x+xyxz^2x+xyxyzx+xyx^2yz+xyx^4-x^3z^2x"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "z^2yx^3-z^2x^3y-zyx^3z+zx^3yz+y^3z^2y-y^4z^2+yxy^4+yxyxyx+2yxyx^3+yx^4y-xyzyxy-xy^5-xy^2xyz-xyxyzy+xyxyxy-xyx^2zy+x^3z^2y-x^3yz^2-x^3yxy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyx^2yxz-xyx^3yx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yx^2yx^2y-xyx^3zx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyx^2yx-xyx^2yxy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyxyx^2+yxyx^2yz+2yxyx^4+yx^4yx-xy^2xyxy+xyx^2yz^2-x^3yxyx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyxyxy-xyx^2y^2z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyxy^3+yxyx^2yz+yxyx^4+xyx^2yz^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyxy^2z-xyx^4y"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxy^2xyx-xyx^2y^3"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxy^3zx-xyx^2y^2z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^2xyx^2y+xyx^4z+x^2yx^2zx+x^4yxz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^2xyxy^2+yx^3yzx-xyxyx^2y+xyx^2y^3-xyx^4y-x^3y^2zx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^3x^2yx+xy^2zx^2y+xyx^2y^2z+x^5yx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^3x^2yz+y^3x^4-xyx^3yz-xyx^3y^2+x^2y^2xyx+x^4yzx+x^5yz+x^7"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^2zx^2y^2-xyx^3yx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yzx^5-xy^2xyx^2+xyxzx^2z-x^3yzx^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yzx^3yx-xyxzx^2y"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yzx^2y^3+yzx^3yz+xy^2xyx^2-xyxzx^2z+x^2yx^2y^2+x^3yzx^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^5y+yzx^2y^2x+y^3zx^2y+xyxy^2xy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^4yx+xyx^2yzy+xyx^2yx^2+xyx^5"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^3yzx-yxyx^2yz-yxyx^4+xyxyxyx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^2y^2xy+zx^4yz+xyxy^3z-xyx^4z-x^2yx^2zx-x^4yxz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^2y^4-zx^4z^2+yzx^2y^2z+y^3zx^2y-xyxy^2z^2+xyxyx^2y+xyx^4y+x^3y^2zx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^2y^2zx+zx^5z+xyxy^2xz+xyx^2yxy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zyx^4y-zx^3yxy+y^3zyxy+yxyxy^2x-yxyx^2yz-2yxyx^4-yx^4yx-xy^5z+xy^2xyxy-xyxzyxy-xyxy^2xy-xyx^2yz^2-xyx^3yz+x^3zyxy+x^3yxyx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zyx^4z-zx^3yxz+y^3zyxz+xyz^2yxz-xyxzyxz+xyx^2y^3+xyx^3y^2+x^3zyxz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zyxy^2xy+zyx^3yz-zx^3y^2z+y^3zy^2z+yxyzyxy+xyz^2y^2z-xyxzy^2z+x^3zy^2z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zyxy^4-zyx^3z^2+zx^3yz^2+yxy^4z-yxyx^2yz+2yxyx^3z-2yxyx^4+yx^4yz-yx^4yx+xyzyx^3-xyzx^3y+xy^4zy+xy^2xyz^2+xy^2xyxy+xyxyz^2y-xyx^2yz^2+xyx^3zy+xyx^3zx+xyx^3yz+xyx^4y-x^2yx^2y^2+x^2yx^3y-x^3yxyz+x^3yxyx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyx^4zx+xyx^4yx+x^2yx^2zx^2+x^4yxzx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyx^3yxz-xyx^4yx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyx^3yzx-xyx^4zy+x^4y^2zx-x^4yxzy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyx^2yxy^2-xyx^4yx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyx^2y^2xz-xyx^3zx^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyx^2y^2zx-xyx^3yxy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyx^2zx^2y-xyx^2yxyz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyxyx^3z+xyx^2zx^3+xyx^5z+x^4yx^2z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyxzx^2zy-xyx^3yx^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yx^2yx^3y-xyx^2y^2z^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yx^2yxyx^2+yx^2yx^4+yx^5yx+xyx^2y^2zy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yx^2y^2xyx+yx^4yzx+xyx^2y^2z^2+xyx^2y^2zy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyx^5+1/2yx^3y^2zx-1/2yx^4yxy+1/2yx^4yx^2-1/2xy^2xyxyx-1/2xyxyx^2yz+1/2xyx^2zx^2z+1/2xyx^2yz^2x-1/2xyx^2y^2z^2-1/2xyx^3yz^2-1/2xyx^4zy-1/2xyx^4yz+1/2x^2yxyx^2y-1/2x^2yx^2y^3+1/2x^2yx^4y-1/2x^3y^3zx+1/2x^3yxyxy-1/2x^3yxyx^2+x^4y^2zx-1/2x^4yxzy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyx^4y+yx^3y^2zx-xyxyx^2yz+xyx^2y^3z-xyx^4yz-x^3y^3zx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyx^4z+yx^2yx^2zx+yx^4yxz-xyxyx^3y-xyx^5y-x^4yx^2y"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyx^3yx-xyx^2yxyz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyx^2y^3+yxyx^3yz-1/2yx^3y^2zx+1/2yx^4yxy-1/2yx^4yx^2+1/2xy^2xyxyx+1/2xyxyx^2yz-1/2xyx^2zx^2z-1/2xyx^2yz^2x+1/2xyx^2y^2z^2+1/2xyx^3yz^2+xyx^3yxy+1/2xyx^4zy+1/2xyx^4yz-1/2x^2yxyx^2y+1/2x^2yx^2y^3-1/2x^2yx^4y+1/2x^3y^3zx-1/2x^3yxyxy+1/2x^3yxyx^2-x^4y^2zx+1/2x^4yxzy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyx^2yzx-xyx^2zx^2z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyx^2yzy-2yx^3y^2zx+yx^4yxy+2xyxyx^2yz+xyx^2yz^2y-xyx^2y^3z+xyx^4zy+2xyx^4yz-x^2yxyx^2y+x^2yx^2y^3-x^2yx^4y+2x^3y^3zx-x^3yxyxy-2x^4y^2zx+x^4yxzy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyx^2yz^2-yx^2yx^2zx-yx^4yxz+xyxyx^3y-xyx^2y^4-xyx^2yxyz-xyx^2yx^3+2xyx^5y+x^4yx^2y"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyxy^2xy-xyx^3y^2z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyxy^2xz+xyx^4yx+x^2yx^2zx^2+x^4yxzx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^2x^3yzx-yxyx^3yz+1/2yx^3y^2zx-1/2yx^4yxy+1/2yx^4yx^2-1/2xy^2xyxyx-3/2xyxyx^2yz-xyxyx^2y^2+1/2xyx^2zx^2z+1/2xyx^2yz^2x-1/2xyx^2y^2z^2-1/2xyx^3yz^2-xyx^3yxy-1/2xyx^4zy-3/2xyx^4yz-xyx^4y^2+1/2x^2yxyx^2y-1/2x^2yx^2y^3+1/2x^2yx^4y-3/2x^3y^3zx+1/2x^3yxyxy-1/2x^3yxyx^2+x^4y^2zx-1/2x^4yxzy-x^4yxy^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^2xyx^4+yx^2yx^2zx-yx^3yxyx+yx^4yxz+xyxyx^2y^2-xyxyx^3y+xyx^2yxyz+xyx^2yx^3-xyx^4z^2+xyx^4y^2-2xyx^5y-x^2yx^2yzx+x^3y^2xyx-x^4yxz^2-x^4yx^2y"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^2xyx^3y+xyxyx^2yz+xyx^4yz+x^4yxyz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^2zx^3yz+yx^3yzx^2+xyx^2y^3x-xyx^2yx^2z+xyx^3yxy+xyx^4yx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yzx^3yzy+yzx^4z^2-xyxyx^2yz+xyx^3y^2z-xyx^3yx^2-xyx^4yx+x^2yx^2y^3-2x^2yx^4z+2x^3yzx^2y-x^3y^3zx-2x^3yx^2zx-2x^5yxz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^4yzx-zx^5zy-xyxy^2xzy+x^2yx^2y^2z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^3yx^2y-yzx^2y^2xz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^3yxyx-yx^3y^2zx+yx^4yxy+xyxyx^2yz+xyx^2yz^2y+xyx^4zy+xyx^4yz-x^2yxyx^2y+x^2yx^2y^2z+x^2yx^2y^3-x^2yx^4y+x^3y^3zx-x^3yxyxy-2x^4y^2zx+x^4yxzy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^3yxzx+zx^6z+y^3zx^3z+xyxy^2x^2z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^3y^2zx-xyx^2y^4-xyx^2yxyz-xyx^2yx^3+xyx^5y+x^2yx^3zx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^2y^3zx+zx^5z^2+xyxy^2xz^2+xyx^2yxyz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zyx^3yzx-xyx^2yxyz-xyx^2yx^3+xyx^3y^3"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zyx^2yxzx+zyx^5z-zx^3yx^2z+y^3zyx^2z+yxyz^2x^2z+xyz^2yx^2z-xyxzyx^2z+x^3zyx^2z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zy^2x^4y-zx^3y^2xy-yxyx^3yz-yx^4y^2z+xyzyxy^2z-xyx^2y^2z^2+xyx^2y^3x+xyx^2y^2x^2-x^2yx^4y+x^3yxy^2z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zy^2x^4z-zx^3y^2xz+y^3zy^2xz+xyz^2y^2xz-xyxzy^2xz+xyx^2yzy^2+xyx^2y^2zy+x^3zy^2xz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "x^4y^3zx-x^4yxyxy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "x^2yx^2zx^2z-x^2yx^2y^2z^2-x^2yx^2y^2zy+x^4yxyzx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyx^4yxy+xyx^5yx+x^2yx^2yxyz+x^4yx^2yx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyx^4yxz-xyx^5yx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyx^4yzx+xyx^5yx+x^2yx^2y^2z^2+x^2yx^2y^2zy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyx^3yx^2y-x^2yx^3zx^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyx^3yxyx-xyx^4zy^2+x^4y^2xyx-x^4yxzy^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyx^3yxy^2-xyx^5yx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyx^3y^2xz+xyx^4yx^2+x^2yx^2zx^3+x^4yxzx^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyx^3y^2zx+xyx^5yx+x^2yx^2yxyz+x^4yx^2yx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyx^3y^2zy+xyx^5yx-x^2yx^4zy+x^3yzx^2y^2-x^3yx^3yx-x^5yxzy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyx^3zx^2y-xyx^3yxyz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyx^2yx^4-x^2yx^2yxyz+x^2yx^4zy-x^4yx^2yx-x^5y^2zx+x^5yxzy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyx^2yx^2zy-xyx^3yx^2z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyx^2yx^2z^2-xyx^4yx^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyx^2yxyxy-xyx^3yxyz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyx^2yxyzx-xyx^3yx^2z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyx^2yxyz^2+xyx^5yx+x^2yx^2yxyz+x^4yx^2yx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyx^2y^2xy^2-xyx^3zx^2z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyx^2y^3zx-xyx^3yxyz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyx^2y^2z^2x-xyx^3zx^2z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "+xyx^2y^2z^2y-xyx^2y^2zy^2-xyx^2y^3z^2+xyx^4zy^2-x^4y^2xyx+x^4yxzy^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyx^2yzx^2y-xyx^3yx^2z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyx^2yzyxz-xyx^3zx^3"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyx^2yzy^2z-xyx^2y^2x^2y"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyx^2zx^4+xyx^3yxyz-x^2yx^2yxyz-x^4yx^2yx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyx^2zx^3y-xyx^5yx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyx^2zx^3z+xyx^4z^2x+x^2yx^2yzx^2+x^4yxz^2x"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyx^2zx^2zy-xyx^4yx^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyx^2zx^2z^2-xyx^2yxyzy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyxyzx^2y^2-xyx^3yx^2z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyxzx^2y^2x+xyxyx^2y^2z+xyx^4z^2y+xyx^4y^2z+xyx^5yx+x^2yx^2yxyx+x^4yxz^2y+x^4yxy^2z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yx^5yxy-xyxyx^3y^2+xyx^2y^2zy^2+xyx^2y^5+2xyx^2y^2xyz+2xyx^2y^2x^3+xyx^3zx^2z-xyx^5y^2+xyx^5yx+x^2yx^2yxyz-x^4yx^2y^2+x^4yx^2yx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yx^4yx^2y-xyxyx^2y^2z+xyx^2y^4z-2xyx^4y^2z-xyx^4yx^2+xyx^5yx+x^2yx^2y^2zy+x^2yx^2yxyz-x^2yx^4zy-x^3yxyx^2y+x^4yx^2yx+x^5y^2zx-x^5yxzy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yx^4yx^2z+xyx^2y^5+xyx^2y^2xyz+xyx^2y^2x^3-xyx^2yxyzy+xyx^3y^4+xyx^3yx^3-xyx^5zx-xyx^5yx-xyx^6y+x^2yx^2yxyz+x^4yx^2yx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yx^4yxyx-xyx^2y^2zy^2-xyx^2y^3z^2+xyx^3zx^2z+xyx^3yz^2x+xyx^3yxyz+xyx^4zyx+x^2yx^2y^3x-x^2yx^4yx-x^3yxyxyx-x^3yx^2yxy-2x^4y^2zx^2+x^4yxzyx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yx^4y^2zx-xyxyx^3y^2+xyx^2y^2xyz+xyx^2y^2x^3+xyx^3zx^2z-xyx^5y^2+xyx^5yx+x^2yx^2yxyz-x^4yx^2y^2+x^4yx^2yx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yx^3yxyx^2-yx^4yxzx-xyxyx^2y^2x+xyx^2yxyzy-xyx^2yx^2zx+xyx^3yxyz-xyx^3yx^2z+xyx^4z^2x-xyx^4zy^2-xyx^4yz^2-xyx^4y^2x+xyx^6z+x^2yxyx^3y+x^2yx^2yzx^2-x^2yx^2y^4-2x^2yx^2yxyz-x^2yx^2yx^3-x^2yx^3yxy+x^2yx^4zy+2x^2yx^5y-x^3y^2xyx^2+x^4yxz^2x-x^4yxzy^2-x^4yxyz^2-x^4yx^2yx-x^5y^2zx+x^5yxzy+x^5yx^2y"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yx^3yxyxy-yx^4yxzy+xyxyx^2y^2z+xyxyx^3yz+xyxyx^3y^2-xyx^2yxyzy-xyx^3zx^2z+xyx^4z^2y+xyx^4y^2z-xyx^4y^3-1/2xyx^4yx^2+2xyx^5y^2+1/2x^2y^2xyxyx+1/2x^2yxyx^2yz-1/2x^2yx^2yz^2x-x^2yx^2y^2z^2-1/2x^2yx^2y^2zy+x^2yx^2yxyx+1/2x^2yx^3yz^2+x^2yx^3yxy+1/2x^2yx^4zy+1/2x^2yx^4yz-x^3y^2xyxy-1/2x^3yxyx^2y+1/2x^3yx^2y^3-1/2x^3yx^4y+x^4yxz^2y+1/2x^4yxyzx+x^4yxy^2z+1/2x^4yxyx^2+x^4yx^2y^2-x^5y^2zx+1/2x^5yxzy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yx^3y^2xyx-yx^4yxy^2-xyx^2y^2zy^2-xyx^2y^3z^2-xyx^2y^4z-xyx^3yz^2y-xyx^4yzy+xyx^5yx+2x^2yxyx^2yz+x^2yxyx^2y^2+x^2yx^2yz^2y-x^2yx^2y^3z-x^2yx^2y^4+x^2yx^2yxyz+x^2yx^4zy+2x^2yx^4yz+x^2yx^4y^2+x^3yxyxy^2-x^3yxyx^2y+x^3yx^2y^3-x^3yx^4y+x^4y^2xyx+x^4yxyxy+x^4yxyx^2+x^4yx^2yx+x^4yx^4-2x^5y^2zx+x^5yxzy+x^7yx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yx^3y^3zx-yx^4yxyz-xyx^2yx^2zx-xyx^4yz^2+xyx^6y+x^2yxyx^2yz+x^2yxyx^3y-x^2yx^2y^3z-x^2yx^2y^4-x^2yx^2yx^3+x^2yx^4yz+2x^2yx^5y+x^3yxyxy^2+x^4y^2xyx+x^4yxyxy+x^4yx^2yx+x^4yx^3y+x^5yx^2y+x^6yzx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yx^3y^2zx^2+xyx^3yxyz-x^2yx^2yxyz-x^3y^3zx^2+x^4yxyzx-x^4yx^2yx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yx^3yzx^3-xyxyx^3yz+xyx^2y^3x^2-xyx^2yxyzy-xyx^2yx^2zx+xyx^4zy^2+3/2xyx^4yx^2-xyx^5yz-xyx^7-1/2x^2y^2xyxyx-1/2x^2yxyx^2yz+1/2x^2yx^2yz^2x+1/2x^2yx^2y^2zy-1/2x^2yx^3yz^2-1/2x^2yx^4zy-1/2x^2yx^4yz+1/2x^3yxyx^2y-1/2x^3yx^2y^3+1/2x^3yx^4y-x^4y^2xyx+x^4yxzy^2-1/2x^4yxyzx-1/2x^4yxyx^2-x^4yx^2yz-x^4yx^4+x^5y^2zx-1/2x^5yxzy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yx^2yx^4y+xyx^3yx^3+xyx^4zy^2+xyx^6y-x^4y^2xyx+x^4yxzy^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yx^2yx^4z+yx^3yx^2zx+yx^5yxz-xyx^3yx^3-xyx^4zy^2-xyx^6y+x^4y^2xyx-x^4yxzy^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yx^2yx^3zx-xyx^3yxyz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yx^2yx^2zx^2+xyx^3yxyz-xyx^5yx+xyx^6z-x^2yx^2yxyz-x^4yx^2yx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyx^3y^2z+yx^3yzx^2y-xyx^3yx^3-xyx^4zy^2+xyx^5yx-xyx^6y+x^2yx^2yxyz-x^2yx^4yx-x^3y^2zx^2y+x^4y^2xyx-x^4yxzy^2+x^4yx^2yx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyx^3yzy-yx^4yxz^2+xyxyx^3yz-xyx^2y^4z-xyx^2yx^3z-xyx^4zyx-xyx^4y^2z-xyx^4yx^2+2xyx^5yz+3xyx^5yx+2x^2yx^2y^2zy+3x^2yx^2yxyz-2x^2yx^4zy+x^4y^2zx^2-x^4yxzyx+x^4yx^2yz+3x^4yx^2yx+2x^5y^2zx-2x^5yxzy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyx^2y^2xy-yx^4yxyz+xyx^2yzyxy-xyx^2y^3z^2+xyx^5yx+2xyx^6y+x^2yxyx^2yz-x^2yx^2y^3z+x^2yx^2yxyz+x^2yx^4yz+x^4yxyxy+x^4yx^2yx+x^4yx^3y"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyxy^2x^3-xyx^2y^2zy^2-xyx^2y^3z^2-xyx^2y^4z-xyx^3yz^2y+xyx^3y^2z^2+xyx^4zy^2+xyx^4yz^2-x^4y^2xyx+x^4yxzy^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "+yxyxy^2x^2y-xyx^3y^3z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^2x^4yx^2-yxyx^2y^2z^2-yxyx^3yz^2+yx^4yxzy-xyxyx^2y^2z-xyxyx^2y^2x-xyxyx^3y^2+xyx^2y^5+xyx^2yx^2zx+xyx^3zx^2z-xyx^4z^2y+xyx^4zy^2+xyx^4yz^2-xyx^4y^2z-xyx^4y^2x-2xyx^5y^2+xyx^5yx-x^2yxyx^3y+x^2yx^2y^2z^2+x^2yx^2y^4+x^2yx^2yxyz-x^2yx^2yxyx+x^2yx^2yx^3+x^2yx^3yxy-2x^2yx^5y-x^3y^2xyx^2-x^4yxz^2y+x^4yxzy^2+x^4yxyz^2-x^4yxy^2z-x^4yx^2y^2-x^5yx^2y"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^2x^4yxy-yxyx^3yz^2-xyxyx^2y^2z+xyxyx^3yz+xyx^2yxyzy-xyx^3yxyz-xyx^4y^2z-xyx^4y^3-1/2xyx^4yx^2+xyx^5zx+xyx^5yx+1/2x^2y^2xyxyx+1/2x^2yxyx^2yz-1/2x^2yx^2yz^2x-1/2x^2yx^2y^2zy+1/2x^2yx^3yz^2+x^2yx^3yxy+1/2x^2yx^4zy+1/2x^2yx^4yz-x^3y^2xyxy-1/2x^3yxyx^2y+1/2x^3yx^2y^3-1/2x^3yx^4y+1/2x^4yxyzx-x^4yxy^2z+1/2x^4yxyx^2-x^5y^2zx+1/2x^5yxzy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^2x^3yxyx-yx^4yxz^2+2xyxyx^3yz-xyx^2y^4z-xyx^2yx^3z-xyx^4zyx-xyx^4yzy-xyx^4y^3-1/2xyx^4yx^2+2xyx^5yz+2xyx^5yx+1/2x^2y^2xyxyx+5/2x^2yxyx^2yz+x^2yx^2yz^2y-1/2x^2yx^2yz^2x+1/2x^2yx^2y^2zy-x^2yx^2y^3z+3x^2yx^2yxyz+1/2x^2yx^3yz^2+x^2yx^3yxy+1/2x^2yx^4zy+5/2x^2yx^4yz-3/2x^3yxyx^2y+3/2x^3yx^2y^3-3/2x^3yx^4y+x^4y^2zx^2-x^4yxzyx+1/2x^4yxyzx-x^4yxy^3+x^4yxyxy+3/2x^4yxyx^2+x^4yx^2yz+3x^4yx^2yx+x^4yx^4-2x^5y^2zx+1/2x^5yxzy+x^7yx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "+y^2x^3y^2zx-yxyx^3yz^2-xyxyx^2y^2z+xyx^2yxyzy-xyx^2yx^2zx-xyx^3yxyz-xyx^4yz^2-xyx^4y^2z+xyx^5zx+x^2yxyx^3y-x^2yx^2y^4-x^2yx^2yxyz-x^2yx^2yx^3+2x^2yx^5y+x^3yxyxy^2+x^4y^2xyx-x^4yxy^2z+x^5yx^2y+x^6yzx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^2x^2yx^2zx+y^2x^4yxz-yx^4yxz^2-xyxyx^2y^2z+xyxyx^3yz-xyx^2y^4z-xyx^2yx^3z-xyx^4zyx-xyx^4y^3+xyx^4yx^2+xyx^5yz-xyx^7-x^2yx^2y^2zy+x^2yx^3yxy+x^2yx^4zy-x^3yxyx^2y+x^4y^2zx^2-x^4yxzyx-x^4yxy^3-x^4yx^4-x^5y^2zx+x^5yxzy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^2xyx^3zx+xyx^5yx+x^2yx^2yxyz+x^5yx^2y"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^3x^5y-xyx^2y^2zyx-xyx^3y^2xy-xyx^4zy^2-x^2yx^3yx^2-x^3yx^4z+x^4yzx^2y+x^4y^2xyx-x^4yxzy^2-x^4yx^2zx-x^6yxz+x^8y"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^3x^4yx-xyx^4yzy-1/2xyx^4yx^2+xyx^5yx-xyx^7-1/2x^2y^2xyxyx+3/2x^2yxyx^2yz+x^2yx^2yz^2y+1/2x^2yx^2yz^2x+1/2x^2yx^2y^2zy-x^2yx^2y^3z+x^2yx^2yxyz-x^2yx^2yxyx-1/2x^2yx^3yz^2+1/2x^2yx^4zy+3/2x^2yx^4yz-1/2x^3yxyx^2y+1/2x^3yx^2y^3-1/2x^3yx^4y-1/2x^4yxyzx+x^4yxyxy-1/2x^4yxyx^2+x^4yx^2yx-x^5y^2zx+1/2x^5yxzy+x^7yx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^3x^2y^2xy+y^3x^4yz-xyx^3y^3z-xyx^3y^2xy+x^2yxyx^3y+x^5y^2xy+x^5yx^2y+x^7yz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^3x^2y^4-y^3x^4z^2+xyx^3y^2z^2-xyx^3y^4-xyx^3yxyz+xyx^6y-x^2yxyx^2yz+x^2yx^3y^2z-x^2yx^4yx-x^3yx^4z+x^4yzx^2y-x^4yxyxy-x^4yx^2zx+x^5y^4-x^6yxz-x^7z^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^3x^2y^2zx+y^3x^5z+xyx^4yx^2+xyx^5yx+x^2yx^2zx^3+x^2yx^2yxyz+x^3yx^2zx^2+x^4yxzx^2+x^4yx^2yx+x^5y^2zx+x^5yxzx+x^8z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^2zx^4yz-xyxyx^2y^2z-2xyx^4y^2z-xyx^4yx^2+xyx^5yx+x^2yx^2y^2zy+x^2yx^2yxyz+x^2yx^3zx^2-x^2yx^4zy-x^3yxyx^2y+x^4yx^2yx+x^5y^2zx-x^5yxzy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^2zx^4z^2-yx^4yxyz-xyx^2y^3z^2-xyx^2yx^2zx-xyx^4yz^2-xyx^5yx+xyx^6y+x^2yxyx^2yz+x^2yxyx^3y-x^2yx^2y^3z-x^2yx^2y^4-x^2yx^2yx^3+x^2yx^4yz+x^2yx^4yx+2x^2yx^5y+x^3y^2zx^2y+x^3yxyxy^2+x^4y^2xyx+x^4yxyxy+x^4yx^2yx+x^4yx^3y+x^5yx^2y+x^6yzx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yzx^3y^2xy-yzx^4y^3+xyxzx^2zx^2+xyx^3y^2z^2-xyx^4zy^2-xyx^5yx+x^2yxyx^2y^2+x^2yx^2y^3z-x^2yx^2y^4-2x^2yx^4z^2+x^2yx^4y^2-2x^3yzx^4+x^3yxyxy^2-2x^3yx^2yzx+4x^4y^2xyx-x^4yxzy^2-2x^5yxz^2+2x^6yzx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yzx^3yz^2x-xyxzx^2z^2x-xyxzx^2y^3-xyxzx^3yz-xyxzx^5+x^2yx^2zx^3+x^2yx^3yxy+x^4yxzx^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yzx^2y^2x^2y-xyx^3zx^2z-xyx^3yx^2z+2xyx^5yx+2x^2yx^2y^2z^2+x^2yx^2y^2zy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yzx^2y^2xzy+xyxzx^2y^3+xyxzx^3yz+xyxzx^5"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^6yz+zx^8-yzx^2y^2xzx-xyxy^2xyzx+xyxy^2x^2yz+xyxy^2x^4+xyx^2y^5-xyx^5y^2+xyx^5yx+x^2yx^2y^2z^2+x^2yx^2yxyz-x^2yx^3zx^2-x^2yx^4yx+x^4yx^2yx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^6zy+yzx^2y^2xzx+y^3zx^3zy+xyxy^2x^2zy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^5zy^2+xyxy^2xzy^2+xyx^2yzy^2x+xyx^5yx-x^2yx^2y^2zy+x^2yx^3zx^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^4y^3x+zx^8+yzx^2y^2z^2x+yzx^4yzy+xyxy^5x+xyxy^2xyzx+xyxy^2x^4-xyxyx^3y^2-xyx^2y^5-xyx^2yxyzy+xyx^5y^2-xyx^5yx-x^2yx^2y^2z^2+x^2yx^3zx^2+x^2yx^4yx+x^3y^3zx^2-x^4yxyzx-x^4yx^2y^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^4y^2zx+xyx^3zx^2z+xyx^3yx^2z-xyx^5yx-x^2yx^2y^2z^2-x^2yx^2y^2zy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^3yx^4-zx^5z^2y+zx^6yx-xyxy^2xz^2y-xyx^2yz^2yx-xyx^2y^2zy^2-xyx^2y^3z^2-xyx^2yxyzy+xyx^3zx^2z+xyx^3yz^2x+xyx^5yx+x^2yx^2yxyz-x^2yx^3yxy+x^4yx^2yx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^3yx^3y+yzx^4yz^2-xyxyx^3yz-x^4yx^2yz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^3yx^2zx-yzx^2y^2x^2z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^3yxyzx+zx^6z^2+y^3zx^3z^2+xyxy^2x^2z^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^3y^2xyx-xyx^2y^5-xyx^2yxyzy+xyx^5y^2-x^2yx^2y^2z^2+x^2yx^4yx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^3y^3zx-xyx^2y^4z-xyx^2yx^3z+xyx^5yz+xyx^5yx+x^2yx^2yxyz+x^2yx^4zy+x^4yx^2yx-x^5y^2zx+x^5yxzy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zyx^5yx-yzx^2y^2xzx+y^3zyx^2yx+xyz^2yx^2yx-xyxzyx^2yx+xyx^3y^2xy-xyx^3yx^3-xyx^4zy^2-xyx^6y+x^3zyx^2yx+x^4y^2xyx-x^4yxzy^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zyx^5yz+zyx^7-zx^5z^2y+zx^6yx-yzx^2y^2xz^2+y^3zyx^2yz+y^3zyx^4+xyz^2yx^2yz+xyz^2x^3yx+xyzyx^3zx-xy^4z^2yx+xy^5z^2x-xy^2xy^4x-xyxzyx^2yz-xyxzyx^4+xyxy^5x-xyxy^2xz^2y+xyxy^2xyzx-xyx^2yz^2yx-xyx^2y^2zy^2-2xyx^2y^3z^2-xyx^2y^3zy-xyx^2yxyzy+2xyx^2yx^2zx-xyx^3z^2yx+xyx^3zx^2z+2xyx^3yz^2x-xyx^3y^2z^2-xyx^4zy^2+4xyx^5yx+x^2yxyx^2y^2-2x^2yxyx^3y+x^2yx^2y^4+3x^2yx^2yxyz+2x^2yx^2yx^3+x^2yx^3y^2x-x^2yx^3yxy+x^2yx^3yx^2-x^2yx^4z^2-x^2yx^4zy+x^2yx^4y^2-4x^2yx^5y+x^3zyx^2yz+x^3zyx^4+x^3yzx^2y^2-x^3yx^2yzx-x^3yx^3yx+2x^4y^2xyx-x^4yxzy^2+x^4yx^2yx-x^5yxz^2-x^5yxzy-2x^5yx^2y"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zyx^5zy-zx^3yx^2zy+y^3zyx^2zy+yxyz^2x^2zy+xyz^2yx^2zy-xyxzyx^2zy-xyx^4yx^2-x^2yx^2zx^3+x^3zyx^2zy-x^4yxzx^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zyx^3yxyx-xyx^2yxyzy+xyx^3y^4-x^2yx^2y^2z^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zyx^3y^3x+zyx^7-zx^3y^4x-zx^5z^2y+zx^6yx+y^3zyx^4+y^6zyx+y^3x^3zyx+yxy^4z^2x+yxyx^3z^2x+xyzyx^3zx+xy^5z^2x-xy^2xy^4x-xyxzyx^4+xyxzx^3yx-xyxy^3zyx-xyxy^2xz^2y+xyxy^2xyzx-2xyx^2yz^2yx-xyx^2y^2zy^2-xyx^2y^3z^2-xyx^2yxyzy+2xyx^2yx^2zx-xyx^3z^2yx+xyx^3zx^2z+3xyx^3yz^2x+xyx^3y^2z^2-xyx^3yxyz+xyx^4yx^2+2xyx^5yx+x^2yxyx^2y^2-2x^2yxyx^3y+x^2yx^2zx^3-x^2yx^2yzyx-x^2yx^2y^2z^2-x^2yx^2y^2zy+x^2yx^2y^4+2x^2yx^2yxyz+2x^2yx^2yx^3-x^2yx^3yxy+x^2yx^3yx^2-x^2yx^4z^2+x^2yx^4y^2-4x^2yx^5y+x^3zyx^4-x^3zx^3yx+x^3y^3zyx-x^3yx^2yzx+x^4y^2xyx+x^4yxzx^2+x^4yx^2yx-x^5yxz^2-2x^5yx^2y+x^6zyx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zyx^3y^2zx-xyx^2yx^3z+xyx^3y^3z+xyx^5yx+x^2yx^2yxyz+x^4yx^2yx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zyx^2yx^2zx-yxyxy^2x^2z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zyx^2yxyzx+zyx^5z^2-zx^3yx^2z^2+y^3zyx^2z^2+yxyz^2x^2z^2+xyz^2yx^2z^2-xyxzyx^2z^2+x^3zyx^2z^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zy^2x^3z^2x-zx^3y^2z^2x+y^3zy^2z^2x-yxy^2xyz^2x-yxy^2x^3zx+xyz^2y^2z^2x-xyxzy^2z^2x+xyx^2y^5+xyx^2y^2x^3-xyx^2yx^2zx-xyx^3zx^2z+x^3zy^2z^2x"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zy^2x^2yxyx+zy^2x^2yx^3+zy^2x^5y-yx^3yx^3y-yx^5y^2z+xyx^2y^5+xyx^2y^2xyz+xyx^2y^2x^3"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "+zy^2x^2yxzx-yxy^2z^2x^2z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zy^2x^2y^2zx+zy^2x^5z+yxy^2zy^2xz-xyx^3yx^3-xyx^4zy^2-xyx^6y+x^4y^2xyx-x^4yxzy^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "z^2yx^2yxyx+z^2yx^2yx^3+zyx^3zx^2y-y^3z^2yx^2y+y^4z^2x^2y+y^2x^4y^2z-yxy^4x^2y+yxyx^2y^2z^2+yxyx^3zx^2+xy^5x^2y+xy^2xyzx^2y+xy^2x^2yxyx+xy^2x^2yx^3+xy^2x^5y+xyxyzyx^2y+xyx^2zyx^2y-xyx^2yx^2zx-xyx^4yz^2+xyx^4y^2z+xyx^4yx^2-3xyx^5yx+x^2yxyx^3y-2x^2yx^2y^2zy-x^2yx^2y^4-3x^2yx^2yxyz-x^2yx^2yx^3+2x^2yx^4zy+2x^2yx^5y-x^3z^2yx^2y+x^3yz^2x^2y-x^3y^2xy^2z-x^4yxyz^2-2x^4yx^2yx-2x^5y^2zx+2x^5yxzy+x^5yx^2y"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "z^2yx^2y^2zx+z^2x^3yx^2z+zyx^3zx^2z-y^3z^2yx^2z+y^4z^2x^2z-y^4zy^2xz-yxy^4x^2z+xy^5x^2z+xy^2x^2y^2zx+xy^2x^2yxzx+xy^2x^5z+xyxyzyx^2z+xyxyx^2y^2z+xyx^2zyx^2z+xyx^2yxyzy+xyx^4y^2z+xyx^5zx+2xyx^5yx-x^3z^2yx^2z+x^3yz^2x^2z-x^3yzy^2xz"),
			},
			complete: false,
		},
		// braid4-9, Example 4.2, Kreuzer, expected basis from bergman.
		{
			ideal: []*Polynomial{
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxy - zyz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyz - zxy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zxz - yzx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "x^3 + y^3 + z^3 + xyz"),
			},
			maxDeg: 9,
			basis: []*Polynomial{
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zxy-xyz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zxz-yzx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zyz-yxy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "z^3+y^3+xyz+x^3"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zy^3+zx^3-y^3z-x^3z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyz^2-xyxyz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^3xy-xy^4-xyx^3+x^4y"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yzx^2z-xyz^2x"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^4+yzx^2y+y^3zx+xyzy^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^2yz-yzx^2y"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zyx^3-zx^3y+y^3zy+yxyxy+xyxyz+x^3zy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zyxyx+y^3xz+xy^2zx+x^4z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zyxyz-yxyxy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zy^2xy-yxy^2z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zy^2zx-yxyxz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "z^2yxy+y^4z+xy^2xy+x^3yz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyx^2z+xy^4x+xyxyzx+xyx^4"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyxyz-x^2yxyz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxy^4+yxyx^3+2x^2yxyz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxy^2z^2-yxyxzy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^4z^2-y^3xzy+x^3yz^2-x^4zy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^3zx-y^3z^2x+yxyxz^2-x^3z^2x"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xy^3xzy+xyx^3z^2-x^3yxyz-2x^4yz^2+x^5zy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yx^2yxyz-x^3yxyz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyxyxy-x^3yxyz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "+yxy^2xyx-xyxy^2zx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxy^2xyz+yxy^2x^3-yxyx^3y-2x^2yxyzy+x^3yxyz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyzyxy-xyxy^2xy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^7+y^4x^3-yx^4yz+x^3y^4-3x^3yxyz+x^3yx^3"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^2zx^2yx-xyz^2x^2z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yzx^3yz-xyz^2x^2y"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yzx^2yxy-x^3yxyz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yzx^2y^3-y^4zx^2-yxyzy^2x-xyz^2x^2z+xyz^2x^2y+xy^2xyxz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yzx^2y^2z-yxyxz^2y-xyzx^3z+3x^3yxyz+2x^4yz^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^3yzx-yxyxy^3-yxyx^2yz-yxyx^4+yx^4yx+2x^2yxyzx-x^3yxyx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^2yxyz-x^3yxyz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^2y^2xy-yxyxz^2y-xyzx^3z+3x^3yxyz+2x^4yz^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^2y^4+zx^2yx^3+y^3zx^2y+x^3yxyz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^2y^2zx-yzx^2yxz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zyxy^2xy-yxyxy^2z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zyxy^2zx-yxyxyxz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "x^4yxyz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "x^4y^2z^2-x^4yxzy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyx^3y^3+xyx^4yz+xyx^6-2x^4y^4-2x^4yx^3"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xy^4x^2z+xyx^5z-x^3yxyzx-x^4yx^2z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyz^2x^2zy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yx^4yz^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yx^3yxyz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyx^3yz+xy^4x^2y+xyxyzx^2y+xyx^5y"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyx^2yxy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyx^2y^3+yxyx^5-xy^4xz^2-xy^4x^2y-xyxyzx^2y-xyxy^3zx-xyx^4z^2-xyx^5y"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyx^2yzx+xyx^2yz^2x+x^3yxyzx+x^4yx^2z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyxyx^3-yxyx^2yzy-yxyx^4y+yx^4yxy-x^3yxyxy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyxy^2xy-x^2yxy^2xy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyxy^2zx-x^2yxy^2zx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyxzyxy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyxzyxz-yx^2y^4x-yx^2yx^4+yx^5yx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyxz^2yx-xyz^2x^2z^2-xy^2xyxz^2-2x^3yxyzx-x^4yz^2x"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxy^2x^3z+xy^4x^2y+xyxyzx^2y+xyx^5y-2x^2yxy^2xy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxy^2zyxy-yxyxzy^2z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^2xyxy^3+y^2xyx^2yz+y^2xyx^4-y^2x^4yx+yx^3yxyx-xyz^2x^2yx-2x^3yxyzx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^3x^2yzy+y^3x^4y-yxyx^2y^2z+xyz^2x^3z-xy^4x^3+xyx^4yz-xyx^6+3x^4y^4+3x^4yx^3+x^5yzy+x^7y"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^3xzyxy+x^4zyxy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^3xzyxz+y^2x^4yx-x^3y^2xyx+2x^3yxyzx+x^4zyxz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^5x^3-y^4x^3y-y^2x^4yz+yx^3y^4+yx^3yx^3+yx^4yzy-x^3y^5+3x^3yxyzy-x^3yx^3y"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^4zyxy-y^3xzy^2z+x^3yzyxy-x^4zy^2z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yzx^2yx^3+2y^4zx^2y+yxy^2xy^2z-xyz^2x^2y^2-xy^2xyxzy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yzx^2yxzy-xyzx^3z^2-2x^4y^4-2x^4yx^3"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^3yxyz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^3y^2zx-yxyxy^3z-yxyx^2yz^2-yxyx^4z+yx^4yxz+2x^2yxy^2zx-x^3yxyxz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^3yz^2x-yzx^2yx^2z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^2y^3zx-yzx^2yxz^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zyx^2yz^2x+yx^2y^4x+yx^2yx^4+x^3yxyzx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zyxy^2x^3+y^3xzx^2y+yxyxyxzy+xy^2zx^3y-2x^3yxyzy+x^4zx^2y"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zyxy^3zx-yxyxyxz^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zy^2x^4y-zx^3y^2xy-yx^4y^2z-2x^2yxy^2xy+x^3yxy^2z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "x^4yx^2zy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "x^4yxyxy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "x^4y^5+x^4y^2xyz+x^4y^2x^3"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xy^2x^4yx-xyx^3yxyx+x^3yxy^2zx+x^4y^2xyx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xy^3xzx^3-xy^3x^4z-xyx^3z^2y^2-x^2y^6z-x^2yx^3y^2z+x^3yxyzy^2+2x^4yz^2y^2+x^5zx^3-x^8z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyzx^3y^3-xy^4zx^3-xy^2xyxz^2y-xyxyzy^2x^2-x^2yz^2x^2zx-2x^4y^4z-2x^4yx^3z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyz^2x^2yxy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "xyz^2x^2z^2y+xy^2xyxz^2y"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yx^4yx^2z+2x^2yx^2yz^2x"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yx^4y^4+yx^4yx^3"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yx^2yxy^2xy-x^3yxy^2xy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yx^2yxy^2zx-x^3yxy^2zx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyx^2yxzy-xyz^2x^3z^2-2x^4y^4z-2x^4yx^3z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyx^2y^2zx+xyx^2y^2xyx+x^3yxy^2zx+x^4yx^2z^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyx^2y^2z^2-xyz^2x^3z^2-2x^4y^4z-2x^4yx^3z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyx^2yz^2x-x^2yx^2yz^2x"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyxy^2x^3-yxyx^2yzy^2-yxyx^4y^2+yx^4yxy^2-x^2yxy^2x^3+x^2yxyx^3y-x^3yxyxy^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyxy^3xz+yxyx^5z-xyx^2y^2xyx-x^3yxy^2zx-x^4yx^2z^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyxy^3zx-x^2yxy^3zx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyxy^3zy+yxyx^2yz^2y+yxyx^4zy+2x^2yxy^2x^3-2x^2yxyx^3y"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyxzy^2z^2+yx^5yxy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxy^2x^4y-yxyx^3yxy-2x^3yxy^2xy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxy^2xy^2xy+yxy^2x^3yz-yxyx^3y^2z-2x^2yxyzy^2z+x^3yxy^2xy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxy^2xy^2zx+yxy^2x^4z-yxyx^3yxz-2x^2yxyzyxz+x^3yxy^2zx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyzx^3z^2+xy^2xyxz^2y"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yxyzy^2xzx+xyz^2x^2z^2x+xyz^2x^2y^3+xyz^2x^3yz-2xy^2xyxz^2x-2xy^2xyx^2yx-x^2yxyzy^2x+x^3yxyzx^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^3x^4yx-y^2x^3yxyx-xy^4x^2yz-xy^4x^4+xyxyzx^2yx-xyx^5yz-xyx^7+2x^3yxyzy^2+x^4yxy^3+x^4yx^2yz+x^4yx^4"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^3x^2y^2xy+y^3x^4yz-xy^4x^3z-xyx^6z+x^4y^4z+x^4yx^3z+x^5y^2xy+x^7yz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^3x^2y^4+y^3x^2yx^3-y^3x^5y+xy^3x^4y-xyx^3y^2xy+x^5y^4+x^5yx^3-x^8y"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^3x^2y^2zx+y^3x^5z+xy^6xz+xyx^3y^2xz+x^5y^2zx+x^8z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^3x^2yz^2x-xyx^3zx^2z+x^5yz^2x"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^3xzx^2zx+yxyxyx^2yz+yxyx^2yzyx+yxyx^4yx-yx^4yxyx+xy^4xzyx+xy^2x^3z^2x-xyx^3yz^2x+xyx^4zyx-x^2y^4xz^2-x^2yx^4z^2+x^3yxyzy^2+x^3yxyxyx+x^4zx^2zx+x^5yxz^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^3xzy^2z^2+y^2x^4yxy-x^3y^2xyxy+x^4zy^2z^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^4x^3z^2-y^3x^2yz^2y-y^3x^3yz^2-y^3x^4zy-xy^6zy-xyx^3y^2zy+x^3yx^3z^2-x^5yz^2y-x^6yz^2-x^7zy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^6xzy+y^3x^2yz^2y+2y^3x^4zy+xy^6zy+xyx^3y^2zy-x^3yx^3z^2+x^5yz^2y+2x^6yz^2+x^7zy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^5zx^2y+1/2y^2xy^2xy^2z+1/2xyz^2x^2zx^2-1/2xyxyzx^2y^2-1/2xyx^2y^4z-1/2xyx^2yx^3z+1/2xyx^5yz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "y^4zx^3y+yxyzy^2x^2y-yxyxz^2y^2z+xyz^2x^3yz-xyzx^3yxy-xy^2xyx^2yz+3x^3yxy^2xy+2x^4yzyxy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yzx^3y^2xy-xyz^2x^2y^2z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yzx^3y^4+yzx^3yx^3+xy^2xyx^2yz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yzx^2yx^2zy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yzx^2y^2x^3+y^4zx^2y^2+yxy^2xy^2zy-xyz^2x^2y^3-xy^2xyxz^2y-xy^2xyxzy^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yzx^2y^2x^2z+xyxyzy^2xz+3xyx^3zx^2z+x^2yz^2x^2z^2-x^5yz^2x"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "yzx^2y^2xzy-yxyzy^2xz^2+xyz^2x^2y^3+xyz^2x^3yz+xy^3zx^2y^2-xy^2xyxz^2x+xy^2xyx^2yz-xy^2xyx^2yx-xy^2x^4yz+xyxyzy^2zy-x^2yz^2x^2yx+x^2y^2xyx^2y-x^2yxyzy^2x+x^4y^2xyx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^3yx^2yx-y^3zyx^2yx-y^3xz^2yxz-y^3xzyx^2z-yxyxyx^2yx-xy^2x^2yz^2x-xyxyzx^2yx-x^2y^5zx-x^2yx^3yzx-x^3zyx^2yx-x^4z^2yxz-x^4zyx^2z+x^5y^2zx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^3y^2xyx-yzx^2yx^2z^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^3y^2xyz+yx^4yxzy-x^3yxyxzy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^3y^4z+zx^3yx^3z-yzx^2yxz^2y+yzx^2yx^2yz+xyzy^2x^2yz+x^2yxyzx^2y"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^3y^3zx-yxyxy^3z^2+yxyx^2yx^3-yxyx^4z^2-yxyx^5y+yx^4yxz^2+xy^4xz^2y+xy^4x^2y^2+xyxyzx^2y^2+xyx^2y^4z+xyx^2yx^3z+xyx^4z^2y-xyx^5yz+xyx^5y^2+2x^2yxy^3zx-x^3yxyxz^2"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^2yx^3zx-yzx^2yx^2yz-yxy^2x^4z+yxyx^3yxz-xyz^2x^2y^2x-xy^2xyxzyx+2x^2yxyzyxz+3x^2yx^2yz^2x-x^3yxyzy^2-x^3yxy^2zx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^2yxy^2xy-x^3yxy^2xy"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zx^2yxy^2zx-x^3yxy^2zx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zyx^2y^2xyx+yx^2y^4xz+yx^2yx^4z+x^3yxy^2zx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zyx^2y^4x+zyx^2yx^4-y^3xzyx^2z-xy^2x^2yz^2x-x^4zyx^2z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zyx^2y^4z+zyx^2yx^3z-zx^3yx^2yz+y^3zyx^2yz-yxyxyxz^2y+yxyxyx^2yz+xyxy^2zx^2y+x^3zyx^2yz"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zy^2x^3yz^2-zy^2x^4zy-zx^3y^3z^2+zx^3y^2xzy-y^3zy^2xzy-y^3zx^3z^2-y^3x^3y^3-y^3x^4yz-y^3x^6+y^2x^3y^4+y^2x^3yx^3-x^3zy^2xzy-x^3zx^3z^2-x^3y^6-x^3y^3x^3-x^4y^4z-x^4yx^3z-x^6y^3-x^9"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zy^2x^3z^2x-zx^3y^2z^2x-yzx^2y^2xz^2+y^3zy^2z^2x-yxyxzx^2zx-xy^4zx^2y-4xy^6zx-xyxyzy^2z^2-3xyx^3y^2zx+x^3zy^2z^2x+x^4y^3zx"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zy^2x^2y^4+zy^2x^2yx^3-zy^2x^5y-yx^2y^5z-yx^2yx^3yz+yx^5y^2z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "zy^2x^2yz^2x-yxyxzyx^2z"),
				parseMust(map[string]Symbol{"x": 1, "y": 2, "z": 3}, Deglex, "z^2yx^2y^4+z^2yx^2yx^3-z^2x^3yx^2y+y^3z^2yx^2y-y^3xzyx^2y-y^2x^4y^2z-xy^2xyzx^2y+xy^2x^2y^4+xy^2x^2yx^3-xy^2x^5y+x^2yxyzx^2y+x^3z^2yx^2y+x^3y^2xy^2z-2x^3yxy^2xy-x^4zyx^2y"),
			},
			complete: false,
		},
	}

	for i, test := range tests {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			t.Parallel()
			basis, complete := BuchbergerHomogeneous(test.ideal, test.maxDeg)
			if len(basis) != len(test.basis) {
				t.Fatalf("%d %v", len(basis), basis)
			}
			for i := range basis {
				if basis[i].Cmp(test.basis[i]) != 0 {
					t.Errorf("%d %v", i, basis[i])
				}
			}
			if complete != test.complete {
				t.Errorf("got %v want %v", complete, test.complete)
			}
		})
	}
}

func TestBuchberger(t *testing.T) {
	tests := []struct {
		ideal    []*Polynomial
		maxiter  int
		basis    []*Polynomial
		complete bool
		long     bool
	}{
		// Example 5.12, Mora.
		{
			ideal: []*Polynomial{
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
			maxiter: 10,
			basis: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1, 2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1, 2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 1, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2}},
				),
			},
			complete: true,
		},
		// Example 4.1.15, Xiu Xingqiang.
		{
			ideal: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{4, 6, 5, 6, 5, 6, 5}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{4, 6, 5, 6, 5}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{4}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{3}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{5, 6, 5, 6, 5, 6, 2}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{5, 6, 5, 6, 2}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1}},
				),
			},
			maxiter: 4,
			basis: []*Polynomial{
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
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{4, 6, 5, 6, 5, 6, 5}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{4, 6, 5, 6, 5}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{4}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{3}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{5, 6, 5, 6, 5, 6, 2}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{5, 6, 5, 6, 2}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1}},
				),
			},
			complete: true,
		},
		// Section 6.5 Polynomials and Rules, NCAlgebra.
		{
			ideal: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 2, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2}},
				),
			},
			maxiter: 13,
			basis: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1, 2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 1, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 1, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2}},
				),
			},
			complete: true,
		},
		// Section 6.6 Polynomials and Rules, NCAlgebra.
		{
			ideal: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 1, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-2, 1), Monomial: Monomial{1, 2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 1, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-2, 1), Monomial: Monomial{2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1}},
				),
			},
			maxiter: 20,
			basis: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 1}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 1}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 2), Monomial: Monomial{1}},
				),
			},
			complete: true,
		},
		// Section 1.2.3 Non-commutative algebras, Bergman manual.
		{
			ideal: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2, 2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 2}},
				),
			},
			maxiter: 9,
			basis: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 2}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1, 1}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 1, 1}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 1, 1}},
				),
			},
			complete: true,
		},
		// Section 6.9.1 Lex Order, NCAlgebra.
		{
			ideal: []*Polynomial{
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2, 3}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{3, 4, 1}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{3, 2, 1, 1}},
				),
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{3, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{}},
				),
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 3}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{}},
				),
			},
			maxiter: 7,
			basis: []*Polynomial{
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 3}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{}},
				),
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{3, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{}},
				),
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{4}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 1}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1, 2, 3, 3}},
				),
			},
			complete: true,
		},
		// Section 2.5 Homogenisation, Bergman manual.
		{
			ideal: []*Polynomial{
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{5, 5}},
					PolynomialTerm{Coefficient: big.NewRat(-2, 1), Monomial: Monomial{}},
				),
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{4, 4}},
					PolynomialTerm{Coefficient: big.NewRat(-3, 1), Monomial: Monomial{}},
				),
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{3, 3}},
					PolynomialTerm{Coefficient: big.NewRat(-5, 1), Monomial: Monomial{}},
				),
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{5}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{4}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{3}},
				),
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{5, 4}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{4, 5}},
				),
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{5, 3}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{3, 5}},
				),
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{5, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2, 5}},
				),
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{4, 3}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{3, 4}},
				),
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{4, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2, 4}},
				),
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{3, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2, 3}},
				),
			},
			maxiter: 28,
			basis: []*Polynomial{
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 2, 2, 2, 2, 2, 2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-40, 1), Monomial: Monomial{2, 2, 2, 2, 2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(352, 1), Monomial: Monomial{2, 2, 2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-960, 1), Monomial: Monomial{2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(576, 1), Monomial: Monomial{}},
				),
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(576, 576), Monomial: Monomial{3}},
					PolynomialTerm{Coefficient: big.NewRat(-5, 576), Monomial: Monomial{2, 2, 2, 2, 2, 2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(194, 576), Monomial: Monomial{2, 2, 2, 2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1520, 576), Monomial: Monomial{2, 2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(2544, 576), Monomial: Monomial{2}},
				),
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(96, 96), Monomial: Monomial{4}},
					PolynomialTerm{Coefficient: big.NewRat(1, 96), Monomial: Monomial{2, 2, 2, 2, 2, 2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-37, 96), Monomial: Monomial{2, 2, 2, 2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(244, 96), Monomial: Monomial{2, 2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-360, 96), Monomial: Monomial{2}},
				),
				NewPolynomial(
					ElimOrder(),
					PolynomialTerm{Coefficient: big.NewRat(576, 576), Monomial: Monomial{5}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 576), Monomial: Monomial{2, 2, 2, 2, 2, 2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(28, 576), Monomial: Monomial{2, 2, 2, 2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(56, 576), Monomial: Monomial{2, 2, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-960, 576), Monomial: Monomial{2}},
				),
			},
			complete: true,
		},
		// G1, Example 4.2.26 Xiu Xingqiang, expected basis from bergman.
		{
			ideal: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "(ababab^2ab^2)^2-1"),
			},
			maxiter: 689,
			basis: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abab-b^2ab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababab^2-bab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2ab-b^2ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ab^2-babab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ababab-bababab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababab^2abab-babab^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababab^2ab-bab^2abababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abababab^2-b^2abababab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2abab^2abab-babab^2abab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2ab^2ab^2ab-ab^2abababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2ab^2aba-b^2ab^2abababababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2abab^2ababa-babab^2abababababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababababababab^2-bab^2ababababab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ababababab^2ab-b^2ababababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2abab^2abab^2abab-abab^2ababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2abababababab^2a-ab^2ab^2ab^2ab^2ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2abab^2abab^2ab^2-ababababababab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2ababababababa-ab^2abab^2abab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2abab^2ab^2a-b^2ababababababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ababababababab^2aba-abab^2abab^2abab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababababababab^2ab-bab^2abababababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2abababababab-babab^2abab^2abababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abababababababab^2-b^2abababababababab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2abab^2abab^2abab^2abab-ab^2abababababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2abab^2abab^2abab^2ab-abababababababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2abab^2abab^2aba-b^2abababababababababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ababababababababab-bab^2abababababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababababababababab^2ab-bababababababababab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2abababababababababa-ab^2abab^2abab^2abab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abababababababababab^2a-abab^2abab^2abab^2abab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababababababababab^2-bab^2abab^2abab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababababababababababab-babababababababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababababababababababab-bababababababababababab^2a"),
			},
			complete: true,
		},
		// G2, Example 4.2.26 Xiu Xingqiang, expected basis from bergman.
		{
			ideal: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "(ababab^2)^3-1"),
			},
			maxiter: 2395,
			basis: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2ababa-abab^2ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2aba-babab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2abab-bab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abab^2ab^2ab-abab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ababab^2aba-ab^2ab^2abab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2abab^2-ababab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ababab^2a-ab^2abab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abab^2a-b^2ababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2ab^2a-bab^2ababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababab^2ab-b^2ab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ababab^2-b^2abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2abab^2ab^2-ab^2ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababababab^2ab-bab^2ababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ababababab^2-b^2ababababab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababababab^2ab^2ab-bab^2ab^2ababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ababababab^2-b^2ababababab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2ab^2ab^2abab-babab^2ab^2ab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ab^2ab^2ab^2ab-bab^2ab^2ab^2ab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2ababab^2-b^2ababab^2ab^2ab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2abab^2abab^2a-bab^2abababababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2abab^2abab^2abab^2-ab^2abababababab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abab^2abab^2abab^2ab^2-abab^2abababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2abababababab^2aba-ab^2abab^2abab^2abab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abababab^2ab^2ab-b^2ababab^2ab^2ab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2ab^2ab^2ab^2aba-babab^2abab^2abab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abababababab^2ab-b^2abab^2abab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ab^2ab^2ababab^2-bab^2ab^2abababab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abababababababab^2aba-ab^2ab^2ababababababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ababababababab^2-abababababababab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2abababababababab^2a-ab^2ababababababab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ababababababab^2a-b^2abababababababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababababababab^2ab^2a-bab^2abababababababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababababababab^2ab-b^2ab^2ababababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abababababababab^2-b^2ababababababab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ababababababab^2ab^2-ab^2abababababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2abababababab^2ab-bab^2ab^2ab^2ab^2abababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2ab^2abababab^2ab^2a-bab^2ab^2abababababababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2abababababab^2-b^2ab^2ab^2ab^2abababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2ab^2ab^2ab^2ab^2abab-ababab^2ab^2ab^2abababab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2ab^2ab^2ab^2ab^2ab^2ab^2ab-abab^2abababab^2ab^2ab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2ab^2ab^2ab^2ab^2ab^2ababa-abab^2ab^2ab^2ab^2ab^2ab^2ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2ab^2ab^2ab^2ababa-b^2ababab^2ab^2ab^2abababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababababababababab^2ab-bab^2abababababababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2ab^2ab^2ab^2ab^2ab^2aba-babab^2ab^2ab^2ab^2ab^2ab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abababababababababab^2-b^2abababababababababab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ababab^2ab^2ab^2abababab^2aba-ab^2ab^2ab^2ab^2ab^2ab^2ab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2ab^2ab^2ababababa-bab^2ab^2ab^2ab^2ababababababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abababababababababab^2-b^2abababababababababab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ababababababababab^2aba-b^2ab^2ab^2ab^2ab^2ab^2ababababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ab^2ab^2ab^2ab^2ab^2abababab-ab^2ab^2ab^2ab^2abababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ab^2ab^2ab^2ab^2abababababa-abab^2ababababababababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2ab^2ab^2abababab^2ab^2-ababab^2ab^2ab^2abababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2ab^2abababababababa-ab^2ab^2ab^2ab^2ab^2ab^2abababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2ab^2abababab^2ab^2a-b^2ababab^2ab^2ab^2abababababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ababab^2ab^2ab^2abababababab^2a-ab^2ab^2ab^2ab^2ab^2abababab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ababababababababababab^2ab-ab^2abababab^2ab^2ab^2abababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2abababababababab^2-b^2ab^2ab^2ab^2ab^2ab^2ab^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababababababababababab^2aba-b^2ab^2abababab^2ab^2ab^2abababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ababababababababababab^2a-b^2abababab^2ab^2ab^2abababab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2abababab^2ab^2ab^2abababab^2a-ab^2ababababababababababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2ab^2ab^2ababababababab-abab^2ab^2ab^2ab^2ab^2ab^2ab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ababababababababab^2ab-bab^2abababab^2ab^2ab^2abababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababab^2ab^2ab^2abababab^2ab^2-bab^2ababababababababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abababab^2ab^2ab^2ababababab-b^2ab^2ab^2ababababababababab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ab^2ab^2abababababababab^2-bab^2ab^2ababababababababab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababababababababababababababab-bababababababababababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2ababababababababab^2-b^2abababab^2ab^2ab^2abababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababab^2ab^2ab^2abababababab^2-b^2ab^2ab^2ab^2ababababababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababababababababababababab-babababababababababababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ababababababababababababab^2aba-ab^2ab^2abababababababababababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2abababababababababababab^2-ababababababababababababab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2abababababababababababab^2aba-abababab^2ab^2ab^2abababababababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ababababababababababababab^2a-ab^2abababababababababababab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2ab^2ab^2ab^2ab^2ab^2ab^2ab-b^2ab^2ab^2ab^2ababababababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2abababababababababab-b^2ab^2ab^2ab^2ab^2ab^2ab^2ab^2ab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ababababababababababab^2-bab^2ab^2ab^2ab^2ab^2ab^2ab^2ab^2ab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abababababababababababab^2a-b^2ababababababababababababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ababababababababababababab-bababababababababababababab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababab^2ab^2ab^2ababababababab-b^2abababab^2ab^2ab^2abababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababababababababababab^2ab^2a-bab^2ababababababababababababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababababababababababababab^2ab-b^2ab^2abababababababababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2ab^2ab^2ab^2ab^2ab^2ab^2ab^2ab^2-b^2ab^2ab^2ababababababababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abababababababababababab^2ab-bababab^2ab^2ab^2ababababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ababababababababababababab^2-b^2abababababababababababab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ab^2ab^2ababababababababab^2-babab^2ab^2ab^2ababababababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababab^2ab^2ab^2ababababababababa-bab^2abababababababababababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababababababababababababab^2ab^2-b^2ab^2abababababababababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2abababababababababababab^2ab^2-ab^2ababababababababababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abababab^2ab^2ab^2abababababababab-ab^2abababababababababababab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2abababababababababababab-babab^2ab^2ab^2ababababababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ab^2ab^2abababababababababab-b^2ab^2ab^2ababababababababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abababababababababababababab^2aba-abab^2ab^2ab^2ababababababababababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2ab^2ab^2ababababababababababab-abababababababababababababab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2ab^2abababababababababababa-b^2abababababababababababababab^2ab"),
			},
			complete: true,
			long:     true,
		},
		// G3, Example 4.2.26 Xiu Xingqiang, expected basis from bergman.
		{
			ideal: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "(abab^2)^2-1"),
			},
			maxiter: 381,
			basis: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab-a^2ba^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2aba-b^2a^2b"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2a-ba^2b^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abab^2-a^2b^2a^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2a^2ba^2-ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ba^2b^2a^2-abab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababa^2b-ba^2baba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2b^2a^2b-b^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2ba^2b^2-bab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "aba^2b^2ab-bab^2a^2ba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2babab^2-b^2ababa^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2ba^2bab-baba^2ba^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2a^2bab^2-b^2aba^2b^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2aba^2bab^2-a^2bababa^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "baba^2ba^2ba-aba^2babab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ba^2b^2ab^2a^2-abababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ba^2bababa^2-aba^2bab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2a^2b^2-b^2a^2b^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2a^2ba-ba^2bababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababababab-bababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2b^2ab^2a^2b-b^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2a^2b^2ab-a^2babababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abababab^2-a^2b^2ab^2a^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2aba^2babab-aba^2ba^2ba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2a^2bababab-a^2ba^2ba^2ba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2a^2b^2ab^2-ababababa^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2a^2babab-aba^2ba^2ba^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2a^2baba^2-a^2bababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "baba^2bababa-ab^2a^2babab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ba^2babababa-ab^2a^2b^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ba^2ba^2ba^2ba-a^2bababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2a^2b^2ab^2a-b^2abababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2a^2bababa-b^2aba^2ba^2b"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "aba^2bababab-b^2abababa^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "aba^2ba^2ba^2b-babababa^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2babababa^2-ba^2ba^2ba^2b"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2ba^2ba^2ba^2-bab^2a^2b^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ababababa^2-ab^2a^2b^2ab^2"),
			},
			complete: true,
		},
		// G4, Example 4.2.26 Xiu Xingqiang, expected basis from bergman.
		{
			ideal: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "(aba^2b^2)^2-1"),
			},
			maxiter: 342,
			basis: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2a^2b-b^2aba^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "aba^2b^2-bab^2a^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2b^2ab-b^2a^2ba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2bab^2-ba^2b^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab-bab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "aba^2bab-baba^2ba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2-b^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2baba^2b-ba^2baba^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "aba^2ba^2ba-bab^2ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2aba^2-b^2aba^2ba^2b"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababa^2b-b^2abababa^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababa^2ba-bab^2ab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2bababab^2-ba^2babab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2a^2-b^2abababa^2b"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababa^2-b^2a^2ba^2ba^2b"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2aba^2ba^2b-b^2a^2bababa^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ababab-bab^2ab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababab^2a^2-b^2a^2bababa^2b"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababab^2-b^2abababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababababa-ba^2babababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababababab^2a-bababababa^2b"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababababab-babababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2babababa^2b-ba^2babababa^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2a-b^2ababababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2bababababa^2-b^2ab^2abababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ababa-b^2a^2babababab^2"),
			},
			complete: true,
		},
		// G5, Example 4.2.26 Xiu Xingqiang, expected basis from bergman.
		{
			ideal: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^5-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "(abab^2)^2-1"),
			},
			maxiter: 183,
			basis: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^5-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab-ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab-bababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abab^2-ab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab-babab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^3-ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^4-abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3abab-babab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^3ab-abab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4abab-babab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^4ab^2-abab^3aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^4ab-ab^2ab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^3ab^2-b^2ab^3ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^3aba-b^3ab^3ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^3ab^2a-b^2ab^3ab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^3ab^3-abab^4aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^3ab^2a-ab^2ab^3ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^3ab^3a-abab^3ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^3ab^2-b^2ab^3aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^3ab^3ab-bab^3ab^3aba"),
			},
			complete: true,
		},
		// G6, Example 4.2.26 Xiu Xingqiang, expected basis from bergman.
		{
			ideal: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^5-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "(ababab^4)^2-1"),
			},
			maxiter: 3095,
			basis: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^5-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^4abab-ab^4aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^4ab-abab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^3abab-babab^3aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2ab-bab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ababa-bab^2ab^3ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^3ab-bab^3ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2aba-bab^3ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^4ab-ab^4ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^3aba-ab^4ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4abab^4-abab^4aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4abab^3a-ab^2ab^4ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^4aba-ab^4abab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3abab^4a-abab^4ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^4ab^4-ababab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^4ab^3a-ab^2abab^4"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^4ab^2ab-ab^3aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ab^4a-abab^2ab^4"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ab^3ab-ababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ab^2aba-ababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^4ab-abab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^3aba-abab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4ab^3ab-b^4ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4ab^2aba-b^4ab^3ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4abab^3-b^2ab^4aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4abab^2a-b^3ab^4ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^4ab-b^4abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3abab^4-bab^4ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^4aba-b^4abab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^4a-bab^4ab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^3ab-bab^3abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^4ab^3-b^2abab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^4ab^2a-b^3abab^4"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^3ab^4-bab^2ab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^3ab^2ab-babab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^4a-bab^3ab^4"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^3ab-bab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4abababab-ab^3ab^3aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babababab^4-abab^3ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3abab^3ab-bab^3abab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^3ab-bab^3ab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4abab^2abab-ab^2ab^3aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ababab^4-ab^4ab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ababab^2ab-ab^3ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ababab^3-ab^2ab^2ab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3abababab-ab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ababab^4-abab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2abab^4-abab^3ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4ababab^3-b^3ababab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4ababab^2a-b^3ab^2abab^4"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ababab^4-b^4ababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^3abab^3-b^3abab^3ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abab^4-b^4abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababab^4a-b^4abab^2ab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababababab^2ab-bab^2ababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababababababa-bab^2ababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2abab^4-ab^4ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2abab^2ab-ab^3ababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2abababa-ab^3ababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4abab^2ab^4-ab^3ab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4abab^2ab^3a-ab^2ababab^4"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3abab^2ab^3-ab^3ab^2ab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3abab^2abab-ab^3abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ababab^2ab-ab^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2abab^2ab^4-ababababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^3abab^3-b^3abab^3ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2abab^4-b^4ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3abab^2ab^4-b^4ab^2abab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^3ab^2ab^2a-bab^2ab^2abab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2abab-babab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababab^2ab^4-bab^2ababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2ab^2ab^4-ab^3ababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2ab^2ab^2ab-ab^3ab^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2ab^2ababa-ab^3ab^2abab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^3abab^2aba-ab^3abab^2abab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3abab^2abab^3a-ab^2ab^3abab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^3abab^2ab^2a-ab^2abab^2abab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ab^2abab^3a-abab^2abab^2ab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2ab^4-abababab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2abab^2aba-ab^2ab^3ab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4ab^2ab^2ab^3-b^3ab^2ab^2ab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4ab^2ab^2abab-b^3ab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^3abab^2ab-b^3abab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ab^2ab^4-b^4ab^2ab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2abab^2aba-b^4ab^2ab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3abab^2abab^3-b^2ab^3abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3abab^2abab^2a-b^3ab^3abab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^3abab^2ab^2-b^2abab^2abab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^3abab^2aba-b^3abab^2abab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2abab^3-bab^2abab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ababab^2ab-bab^2ababab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2abab^3a-b^2ab^3abab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2abab^2ab-bab^2abab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2ab^3a-b^2ab^2ab^2abab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^3ab^2ab^4-ab^3abab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2ab^3ab^4-ab^3ab^2abab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2ab^3ab^2ab-ab^3ab^3abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ab^2ab^2ab^3a-ab^2ababab^2ab^4"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ab^2ab^2ab^2ab-ab^2abababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ab^2ab^2ababa-ab^2abababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^3ab^2ab^4-ab^2ab^3ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ababababab^3-ab^2ab^2ababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4ab^2ab^3ab^3-b^3ab^3ab^2ab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4ab^2ab^3ab^2a-b^2ab^2abab^2abab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^3ab^2ab^4-b^4ab^2ab^3ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ababababab^3-b^3ababababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^3ab^2ab^4a-b^3abab^2abab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2abab^2ab-bab^2abababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2abababa-bab^2abababababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abab^2abab^2-b^2abab^2abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^3ab^2ab^3ab-ab^3abababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^3ab^2ab^2aba-ab^3abababababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2abab^2abab^3-ab^3abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3abab^2abab^2ab^3-ab^2ab^2abab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3abab^2abab^2ab^2a-ab^2ab^3ab^2ab^4"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ab^2ab^3ab^4-ab^2abababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ab^2ab^3ab^3a-ab^3ab^2ab^2abab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ab^2ab^3ab^2ab-ab^3ab^2ab^3ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^3ab^2ab^3ab-ab^2ab^3ab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2abababab^3-ab^2ab^2abab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^3ab^2ab^3ab-b^4ab^2ab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ab^2abab^2ab-bab^3ab^2ab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2abababab^3-b^3abababab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3abababab^2ab^3-b^3ab^2abababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2ab^2ab-bab^2ab^2ab^2ab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^3ab^2ab^3ab^2a-ab^3ab^2ab^3ab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2ab^3ab^2ab^3a-ab^2ab^3ab^2ab^3ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2ab^2ab^2ababab-ab^3ab^2ab^2ab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2ababab^3-ab^2ab^2ab^2ab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2abababab^2ab-ab^2ab^2abababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2abababab^2ab^3a-ab^2abababababab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ab^3ab^2ab^3-b^2ab^3ab^2ab^3ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ab^2ab^2abab^2-b^3ab^2ab^2ababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ab^2ababab^3-b^3ab^2ab^2ab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2abababab^2ab-b^3abababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3abababababab^2a-b^3ab^2abababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^3ab^2ab^3ab^3-b^2ab^2ab^3ab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^3ab^2ab^3ab^2a-b^3ab^2ab^3ab^2ab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abababababab-bab^2ab^2abababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababababab^3a-bab^2abababab^2ab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2ab^2abab^2ab^3-ab^2ab^2abababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2ababab^2ab^4-ab^3ab^2ab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^3ab^2ab^2abab^3-ab^3abababab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2ab^2abababab^2a-ab^2ab^2abababab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2abababab^2ab^2a-ab^3ab^2ab^2abababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ab^2ab^2ab^2abab^2-ab^2ab^2ababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2ab^2abab^2a-ab^3ab^3ab^2ab^2ab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2ab^2ababab-ab^2ab^2ab^2ab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^3ab^2ab^2abab^2-b^3ab^2ab^2abab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ab^2ab^2ab^2ab^2-b^2ab^2ab^2ab^2ab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ab^2abababab^2-b^2ab^2abababab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ab^2ababababa-b^3ab^2abababab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ababab^2ab^4-b^4ab^2ababab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ababab^2ab^3a-b^3ab^2ab^2ab^2ab^2ab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abababab^2ab^3-b^2ab^2ab^2abababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abababab^2ab^2a-b^3ab^2ab^2abababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2ab^2ab^2ab^2ab^3a-ab^3ab^2ababab^2ab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2ab^2ab^2ab^2ababa-ab^2ab^2ab^2ab^2ababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2ab^2ababababab^2-ab^2ab^2ab^2ab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2abababababab^3-ab^3ab^2ab^2ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ab^2ab^2ab^2abab-b^2ab^2ab^2ab^2ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2abababababab^2-b^2ab^2ab^2ababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2ababab^2a-b^3ab^2ab^2ab^2ab^2abab"),
			},
			complete: true,
			long:     true,
		},
		// G7, Example 4.2.26 Xiu Xingqiang, expected basis from bergman.
		{
			ideal: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^5-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "(abab^2ab^4)^2-1"),
			},
			maxiter: 5323,
			basis: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^5-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^4abab-ab^4aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^4ab-abab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^3abab-babab^3aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2ab-bab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ababa-bab^2ab^3ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^3ab-bab^3ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2aba-bab^3ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^4ab-ab^4ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^3aba-ab^4ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4abab^4-abab^4aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4abab^3a-ab^2ab^4ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^4aba-ab^4abab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3abab^4a-abab^4ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^4ab^4-ababab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^4ab^3a-ab^2abab^4"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^4ab^2ab-ab^3aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ab^4a-abab^2ab^4"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ab^3ab-ababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ab^2aba-ababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^4ab-abab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^3aba-abab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4ab^3ab-b^4ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4ab^2aba-b^4ab^3ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4abab^3-b^2ab^4aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4abab^2a-b^3ab^4ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^4ab-b^4abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3abab^4-bab^4ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^4aba-b^4abab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^4a-bab^4ab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^3ab-bab^3abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^4ab^3-b^2abab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^4ab^2a-b^3abab^4"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^3ab^4-bab^2ab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^3ab^2ab-babab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^4a-bab^3ab^4"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^3ab-bab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4abababab-ab^3ab^3aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babababab^4-abab^3ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3abab^3ab-bab^3abab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^3ab-bab^3ab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4abab^2abab-ab^2ab^3aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ababab^4-ab^4ab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ababab^2ab-ab^3ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ababab^3-ab^2ab^2ab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3abababab-ab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ababab^4-abab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2abab^4-abab^3ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4ababab^3-b^3ababab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4ababab^2a-b^3ab^2abab^4"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ababab^4-b^4ababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^3abab^3-b^3abab^3ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abab^4-b^4abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababab^4a-b^4abab^2ab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababababab^2ab-bab^2ababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababababababa-bab^2ababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2abab^4-ab^4ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2abab^2ab-ab^3ababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2abababa-ab^3ababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4abab^2ab^4-ab^3ab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4abab^2ab^3a-ab^2ababab^4"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3abab^2ab^3-ab^3ab^2ab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3abab^2abab-ab^3abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ababab^2ab-ab^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2abab^2ab^4-ababababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^3abab^3-b^3abab^3ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2abab^4-b^4ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3abab^2ab^4-b^4ab^2abab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^3ab^2ab^2a-bab^2ab^2abab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2abab-babab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababab^2ab^4-bab^2ababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2ab^2ab^4-ab^3ababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2ab^2ab^2ab-ab^3ab^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2ab^2ababa-ab^3ab^2abab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^3abab^2aba-ab^3abab^2abab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3abab^2abab^3a-ab^2ab^3abab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^3abab^2ab^2a-ab^2abab^2abab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ab^2abab^3a-abab^2abab^2ab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2ab^4-abababab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2abab^2aba-ab^2ab^3ab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4ab^2ab^2ab^3-b^3ab^2ab^2ab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4ab^2ab^2abab-b^3ab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^3abab^2ab-b^3abab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ab^2ab^4-b^4ab^2ab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2abab^2aba-b^4ab^2ab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3abab^2abab^3-b^2ab^3abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3abab^2abab^2a-b^3ab^3abab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^3abab^2ab^2-b^2abab^2abab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^3abab^2aba-b^3abab^2abab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2abab^3-bab^2abab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ababab^2ab-bab^2ababab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2abab^3a-b^2ab^3abab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2abab^2ab-bab^2abab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2ab^3a-b^2ab^2ab^2abab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^3ab^2ab^4-ab^3abab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2ab^3ab^4-ab^3ab^2abab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2ab^3ab^2ab-ab^3ab^3abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ab^2ab^2ab^3a-ab^2ababab^2ab^4"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ab^2ab^2ab^2ab-ab^2abababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ab^2ab^2ababa-ab^2abababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^3ab^2ab^4-ab^2ab^3ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ababababab^3-ab^2ab^2ababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4ab^2ab^3ab^3-b^3ab^3ab^2ab^4a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^4ab^2ab^3ab^2a-b^2ab^2abab^2abab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^3ab^2ab^4-b^4ab^2ab^3ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ababababab^3-b^3ababababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^3ab^2ab^4a-b^3abab^2abab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2abab^2ab-bab^2abababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2abababa-bab^2abababababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abab^2abab^2-b^2abab^2abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^3ab^2ab^3ab-ab^3abababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^3ab^2ab^2aba-ab^3abababababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2abab^2abab^3-ab^3abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3abab^2abab^2ab^3-ab^2ab^2abab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3abab^2abab^2ab^2a-ab^2ab^3ab^2ab^4"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ab^2ab^3ab^4-ab^2abababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ab^2ab^3ab^3a-ab^3ab^2ab^2abab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ab^2ab^3ab^2ab-ab^3ab^2ab^3ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^3ab^2ab^3ab-ab^2ab^3ab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2abababab^3-ab^2ab^2abab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^3ab^2ab^3ab-b^4ab^2ab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ab^2abab^2ab-bab^3ab^2ab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2abababab^3-b^3abababab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3abababab^2ab^3-b^3ab^2abababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2ab^2ab-bab^2ab^2ab^2ab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^3ab^2ab^3ab^2a-ab^3ab^2ab^3ab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2ab^3ab^2ab^3a-ab^2ab^3ab^2ab^3ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2ab^2ab^2ababab-ab^3ab^2ab^2ab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2ababab^3-ab^2ab^2ab^2ab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2abababab^2ab-ab^2ab^2abababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2abababab^2ab^3a-ab^2abababababab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ab^3ab^2ab^3-b^2ab^3ab^2ab^3ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ab^2ab^2abab^2-b^3ab^2ab^2ababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ab^2ababab^3-b^3ab^2ab^2ab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2abababab^2ab-b^3abababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3abababababab^2a-b^3ab^2abababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^3ab^2ab^3ab^3-b^2ab^2ab^3ab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^3ab^2ab^3ab^2a-b^3ab^2ab^3ab^2ab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abababababab-bab^2ab^2abababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababababab^3a-bab^2abababab^2ab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2ab^2abab^2ab^3-ab^2ab^2abababab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4ab^2ababab^2ab^4-ab^3ab^2ab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^3ab^2ab^2abab^3-ab^3abababab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2ab^2abababab^2a-ab^2ab^2abababab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2abababab^2ab^2a-ab^3ab^2ab^2abababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ab^2ab^2ab^2abab^2-ab^2ab^2ababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2ab^2abab^2a-ab^3ab^3ab^2ab^2ab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2ab^2ababab-ab^2ab^2ab^2ab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^3ab^2ab^2abab^2-b^3ab^2ab^2abab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ab^2ab^2ab^2ab^2-b^2ab^2ab^2ab^2ab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ab^2abababab^2-b^2ab^2abababab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ab^2ababababa-b^3ab^2abababab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ababab^2ab^4-b^4ab^2ababab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ababab^2ab^3a-b^3ab^2ab^2ab^2ab^2ab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abababab^2ab^3-b^2ab^2ab^2abababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abababab^2ab^2a-b^3ab^2ab^2abababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2ab^2ab^2ab^2ab^3a-ab^3ab^2ababab^2ab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2ab^2ab^2ab^2ababa-ab^2ab^2ab^2ab^2ababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2ab^2ababababab^2-ab^2ab^2ab^2ab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2abababababab^3-ab^3ab^2ab^2ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ab^2ab^2ab^2abab-b^2ab^2ab^2ab^2ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2abababababab^2-b^2ab^2ab^2ababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2ababab^2a-b^3ab^2ab^2ab^2ab^2abab"),
			},
			complete: true,
			long:     true,
		},
		// G8, Example 4.2.26 Xiu Xingqiang, expected basis from bergman.
		{
			ideal: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "(ababab^3)^2-1"),
			},
			maxiter: 309,
			basis: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^4-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3abab-ab^3aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^3ab-abab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^2ab-b^3ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3abab^2-b^2ab^3aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ababa-b^3ab^3ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^3ab-b^3abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^3-bab^3ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^3ab^2-b^2abab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^3aba-b^3abab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^3-bab^2ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab-babab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^3a-bab^3ab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ab-bab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababababa-bab^2ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^3aba-ab^3abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3abab^3a-abab^3ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^3ab^3a-ababab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2aba-abababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^3ab^3ab-b^3ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2ab-bab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2ab-babababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ababab^3-ab^3ab^3a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3abababab-ab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babababab^3-abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^3a-bab^2ababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab-bab^2ab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abab^2-b^2abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababab^3-b^3ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3abab^2ab^2a-abababab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ababab^2a-ab^2ab^2ab^3"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2abababab-ab^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ababab-b^2abababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababab^2-b^2ab^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2ab^2ab^3-ab^2ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3ab^2abababa-ab^2abababab"),
			},
			complete: true,
		},
		// G9, Example 4.2.26 Xiu Xingqiang, expected basis from bergman.
		{
			ideal: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "(abab^2)^2-1"),
			},
			maxiter: 30,
			basis: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab-b^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2-bab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab-bababa"),
			},
			complete: true,
		},
		// G10, Example 4.2.26 Xiu Xingqiang, expected basis from bergman.
		{
			ideal: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "(ababab^2)^2-1"),
			},
			maxiter: 115,
			basis: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2aba-ab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abab^2a-abab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2a-ababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2abab-ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2ab-abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab-b^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2-bab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababa-b^2ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2-babab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2aba-b^2abab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2a-bab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababab-babababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ababab^2-ab^2ab^2a"),
			},
			complete: true,
		},
		// G11, Example 4.2.26 Xiu Xingqiang, expected basis from bergman.
		{
			ideal: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "(abababab^2)^2-1"),
			},
			maxiter: 325,
			basis: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ababab-ab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2abab-ab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bababab^2ab-abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababa-b^2ab^2ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ababa-b^2ab^2abab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2aba-b^2abab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababab^2a-bab^2ab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababababab-bababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ab^2aba-ab^2ababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2abab^2a-abab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abab^2ab^2a-ababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abababab^2-ab^2ab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2a-abababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab-b^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abab^2-bab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2ab^2-babab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2ab^2-bababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2ab-bab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ab^2abab-babab^2ab^2ababa"),
			},
			complete: true,
		},
		// G12, Example 4.2.26 Xiu Xingqiang, expected basis from bergman.
		{
			ideal: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "(ababab^2abab^2)^2-1"),
			},
			maxiter: 1256,
			basis: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2abab^2aba-ab^2abab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abab^2ab^2aba-ab^2ababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abab^2abab^2a-abab^2abab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2abab^2a-abab^2ababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2abab^2ab^2a-ababab^2abab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2abab^2abab-ab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ababab^2ab-abab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2abab^2ab-abab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abab^2ab-b^2abab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2ab^2ab-b^2ababab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2abab^2-bab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2ababa-b^2ab^2abab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababab^2aba-b^2abab^2ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2abab^2-bab^2ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2ab^2-babab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2aba-b^2abab^2abab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ababab^2a-bab^2ab^2abab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abababab-babababab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2abab^2a-bab^2abab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababababab^2ab-bab^2ababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abab^2ababab^2-ab^2ab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ababab^2abab^2-ab^2abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababababab^2-b^2ababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababab^2abab^2-b^2abab^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2ab^2ab^2ab-bab^2ab^2ab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2abab^2-b^2abab^2ab^2ab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abab^2ab^2ab^2ababa-ab^2ab^2abababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ababab^2ab^2ab^2aba-ab^2abababab^2ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2ababab^2a-abab^2ab^2abababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2abababab^2ab-abab^2ab^2ab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2abababab^2ab^2ab-ababab^2ab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abababab^2aba-b^2abab^2ab^2ab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababab^2ab^2ab^2ab-b^2abababab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababab^2ab^2aba-b^2ababab^2ab^2ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababab^2ababab-bab^2abababab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2ab^2ababab^2-bab^2ab^2abababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2abababab^2a-bab^2ab^2ab^2ababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ab^2ab^2abab^2-bab^2abababab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ababababababab^2-ababab^2ab^2ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababababababab^2a-b^2ababab^2ab^2ababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2ab^2abababab-babababab^2ab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ababab^2ab^2ababab^2a-ab^2ababababababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababab^2ab^2ababab^2-b^2ababababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababababababababab-babababababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abababab^2ab^2abababa-ab^2ab^2ab^2ab^2ab^2ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2ab^2ab^2ab^2ab-abababab^2ab^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2ab^2ab^2aba-b^2abababab^2ab^2ababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2abababab^2-b^2abababab^2ab^2ab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababab^2ab^2ab^2ab^2-b^2ab^2ab^2ab^2abababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababababababababab-bababababababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababab^2ab^2abababab^2-bab^2ab^2ab^2ab^2ab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abababab^2ab^2ab^2ababa-abababab^2ab^2ab^2ababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babababab^2ab^2ab^2ababab-abababab^2ab^2ab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2ab^2abababa-bab^2ab^2ab^2abababababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abababababababab-bababab^2ab^2abababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababab^2ab^2ab^2abab-bababab^2ab^2ab^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababab^2ab^2ab^2abababa-b^2abababab^2ab^2ab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababab^2ab^2ababababab-b^2ab^2ababababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ab^2ab^2ab^2ab^2ababab-ab^2ab^2ab^2ababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2ab^2ab^2ababab^2-abab^2ab^2abababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2ababababababa-ab^2ab^2ab^2ab^2ab^2ababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ababababababab-babab^2ab^2abababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ab^2abababababab^2-b^2ab^2ab^2abababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2abababababab-bab^2ab^2ab^2ab^2abababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ababab^2ab^2ababababababa-ababab^2ab^2ababababababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bababab^2ab^2ababababababab-ababab^2ab^2ababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababab^2ab^2abababababab-babab^2ab^2abababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ab^2abababababababa-b^2ababab^2ab^2abababababab"),
			},
			complete: true,
		},
		// G13, Example 4.2.26 Xiu Xingqiang, expected basis from bergman.
		{
			ideal: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "(ababababab^2ab^2)^2-1"),
			},
			maxiter: 7823,
			basis: []*Polynomial{
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "a^2-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^3-1"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2abababab-ab^2ab^2ab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2ab^2ababab-ab^2ab^2ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bababab^2ab^2abab-ab^2ababab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babababab^2ab^2ab-ababab^2ab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ababababa-b^2ab^2ab^2ab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababababab^2a-bab^2ab^2ab^2ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2abababa-b^2ab^2ab^2ababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ab^2ababa-b^2ab^2ababab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababab^2ab^2aba-b^2ababab^2ab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababababab^2ab^2a-babab^2ab^2ab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ab^2ab^2ababa-ab^2ab^2abababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ab^2ababab^2a-abab^2ab^2ababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ababab^2ab^2a-ababab^2ab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ababababab^2-ab^2ab^2ab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ababab^2ab^2ab^2a-abababab^2ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ababababab^2ab^2-abab^2ab^2ab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2ab^2aba-ab^2ababababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2ab^2ab^2ab^2a-ababababab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2abab-b^2ab^2ababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ababab^2-bab^2ab^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ababab^2ab^2-babab^2ab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababab^2ab^2ab^2-bababab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2ab^2ab^2ab-b^2ababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ab^2ab^2ab^2-babababab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababab^2abababab-babababab^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababababab^2ababab-bababab^2ababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ab^2abab^2ababa-ab^2ab^2abab^2ababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2abab^2ababab^2a-abab^2ab^2abab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ababab^2abab^2a-abab^2abab^2ab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abab^2ab^2ab^2ababa-ab^2ab^2abababab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abab^2ababab^2ab^2a-ababab^2ab^2abab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ababab^2abab^2ab^2a-ababab^2abab^2ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2abab^2aba-ab^2abab^2abababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2abab^2ababab-ab^2ab^2abab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2abab^2ab^2ab^2aba-ab^2abababab^2abab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2ab^2ab^2abab^2a-abab^2abababab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2ab^2abab^2abab-ab^2abab^2ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2abab^2ab^2ab^2a-abababab^2abab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2abab^2ab^2abab-ab^2ababab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2abab^2abababa-ab^2ababababababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bababab^2abab^2ab^2ab-ababab^2abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bababab^2abab^2ababa-abababababababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2abab^2abab-b^2ab^2abab^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abab^2ababab^2-bab^2ab^2abab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abab^2abababa-b^2ab^2ab^2abab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ababab^2abab^2-bab^2abab^2ab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abababab^2aba-b^2abab^2ab^2ab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2ab^2ab^2abab-b^2ab^2abababab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2ababab^2ab^2-babab^2ab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2abababab^2a-bab^2ab^2ab^2abab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababab^2abab^2ab^2-babab^2abab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababab^2abab^2a-bab^2abab^2ab^2ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababab^2ababab-bababab^2abababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abababababababa-babab^2abab^2ababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2ab^2abab^2ab-b^2abab^2abababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2abab^2ababa-b^2ab^2abab^2ababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2ab^2ab^2ab-b^2abababab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2ab^2ababa-b^2ab^2ababab^2abab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abababab^2ab^2a-babab^2ab^2ab^2abab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ab^2ab^2abab^2-bab^2abababab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ab^2abab^2aba-b^2abab^2ababab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2abab^2ab^2ab^2-bababab^2abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2abab^2ab^2aba-b^2ababab^2abab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2abab^2ababab-b^2abababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ababababab^2-b^2ababababab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababab^2abab^2ab^2a-babab^2abab^2ab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababab^2abab^2abab-bababababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababababababab^2a-bababab^2abab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2abab^2abababab^2-ab^2ab^2ab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2abababab^2abab^2-ab^2abab^2ab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ababababababab-abab^2abab^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abab^2abababab^2ab^2-abab^2ab^2ab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abababab^2abab^2ab^2-abab^2abab^2ab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abababab^2abab^2aba-ababababababab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abababababababab^2-ababab^2abab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bababababababab^2ab^2-abababab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2abababab^2-b^2ab^2ababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2abababab^2ab^2-b^2ab^2abababab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ababab^2ababab^2-bab^2ab^2abab^2ab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababab^2ababab^2ab^2-babab^2ab^2abab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2abab^2ab^2abab-b^2ab^2ababab^2ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ab^2abab^2ab^2ab-b^2ababab^2ababab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2ab^2abab^2aba-b^2ababababababababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2ababababab^2-b^2ababababab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababab^2ab^2ab^2ababab-bababab^2ab^2ab^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2abab^2ab^2ab^2ab-bab^2ab^2ab^2abab^2ab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababab^2ababab^2ababab-bababab^2ababab^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2ab^2ab^2ab^2ab^2a-ab^2ab^2abab^2abab^2ab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2ab^2abab^2ab^2ababa-ab^2ab^2ababab^2abababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bababab^2abab^2abab^2ababa-ababab^2ababab^2abababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2ab^2ab^2ab^2ab^2ab-b^2ab^2ab^2abab^2abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ababab^2abababab^2a-bab^2ab^2ab^2abab^2ab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2abab^2abab^2ababab-babababab^2ababab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ababab^2abababab^2-b^2abababab^2ababab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ababab^2ababababa-bababab^2abab^2abab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ab^2abab^2abab^2ab^2ab^2-ab^2ab^2ab^2ab^2ab^2ab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ababab^2ababab^2abababab-ababab^2abab^2abab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abababab^2ababab^2ab^2aba-ababab^2abababababab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2abab^2abab^2ababab^2a-ab^2ababab^2ababab^2ababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bababab^2abababababab^2ab^2-abababab^2ababab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2ababab^2abab^2abab-babab^2abab^2ababab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababab^2abab^2abab^2abab-bababab^2ababab^2ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababab^2ababab^2abababa-babab^2abab^2abab^2ababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2abab^2abab^2ab^2ab^2-b^2ab^2ab^2abab^2abab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2ababab^2abababab^2-b^2ab^2abababababab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ababab^2ababab^2ababab-abab^2abab^2abab^2ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bababab^2abababababab^2abab-ab^2ababab^2abababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2abababababab^2ababa-b^2ab^2ababab^2abababababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ababab^2abababababab^2a-ababab^2abababababab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ababab^2abababababab^2-babab^2abababababab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababab^2ababab^2abab^2abab-babab^2abab^2ababab^2ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ab^2abab^2ab^2abab^2ab^2aba-abab^2ab^2abab^2abab^2abab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2abab^2ab^2abab^2ab^2ab^2a-ab^2ab^2abab^2abab^2abab^2ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2ababab^2abababab^2abab^2-ab^2ababab^2abab^2abab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ab^2abab^2ab^2abab^2ab^2ab-bab^2ab^2abab^2abab^2abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abab^2ab^2abab^2ab^2ab^2ab-bab^2ab^2ab^2abab^2ab^2abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abab^2abab^2abab^2ab^2aba-bab^2ab^2abab^2ab^2abab^2ab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2abab^2abab^2abab^2ab^2a-b^2ab^2ab^2abab^2ab^2abab^2ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2abab^2ababab^2abab^2-b^2abab^2ababab^2abab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ab^2ab^2abab^2ab^2abab^2ab^2a-ab^2ab^2abab^2abab^2abab^2ab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ab^2abab^2abab^2abab^2ab^2ab-ab^2ab^2abab^2ab^2abab^2ab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ababab^2ababababababab^2a-abab^2abab^2abababababab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abab^2ababab^2abababababab^2a-ababab^2abababababab^2abab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2abab^2abab^2abab^2ab^2ab^2-ab^2ab^2ab^2abab^2ab^2abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bababab^2ababababababab^2ababa-ababababababababababababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2ababab^2ababababababab^2-bab^2abab^2abababababab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2ababab^2abababababab^2-babab^2abababababab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2abababababab^2ababa-b^2ab^2ababab^2ababababababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2abab^2abab^2abab^2ababab-bababab^2abab^2abab^2abab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2abababababab^2abab^2aba-b^2abab^2ababab^2abababababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ab^2abab^2abab^2abab^2ababab-abab^2ab^2ababab^2ababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ababababababababababababab^2-ababab^2ababababababab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bababab^2abab^2abab^2abab^2ab^2ab^2-abab^2ababab^2abababababab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2abab^2abab^2ababab^2abab^2-b^2abab^2ababab^2abab^2abab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2abab^2abab^2abab^2abab^2ababab-abab^2abab^2ababababababab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2abab^2ababab^2ababab^2ababab-ab^2abab^2ababab^2ababab^2ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bababab^2abab^2abab^2abab^2abab^2ab^2-ababab^2ababababababab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abab^2abab^2abab^2abab^2ab^2ab-bab^2ab^2abab^2abab^2abab^2abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abab^2abab^2abab^2abab^2ababa-bab^2abab^2abab^2ababab^2ababababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2abab^2abab^2abab^2abab^2ab^2-b^2ab^2abab^2abab^2abab^2abab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2abababababab^2abab^2aba-b^2abababababababababababababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2abab^2abab^2abab^2abab^2ab^2a-b^2ababab^2ababababababab^2abab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ab^2abab^2abab^2abab^2abab^2abab-ab^2abab^2abab^2ababab^2abababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ab^2abab^2abab^2ababab^2abababa-ab^2ab^2abab^2abab^2ababab^2abababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2abab^2abab^2ababab^2abababab^2a-abab^2ab^2abab^2abab^2ababab^2ababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ababab^2ababababababab^2abab^2aba-ababab^2abab^2abab^2abab^2abab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abababab^2ababab^2abab^2abab^2ab^2a-abababab^2ababab^2abab^2abab^2ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2abab^2abab^2ababab^2abab^2aba-ab^2abab^2ababababababab^2abab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2ab^2abab^2abab^2ababab^2abababab-ab^2ab^2abab^2abab^2ababab^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2abab^2abab^2ababab^2abababababa-ab^2ab^2abab^2abab^2abab^2abab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2abab^2abab^2abab^2abab^2ab^2ab^2-ab^2ababab^2ababababababab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babababab^2ababab^2abab^2abab^2ab^2ab-abababab^2ababab^2abab^2abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abab^2abab^2ababab^2abababab^2-bab^2ab^2abab^2abab^2ababab^2abababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abab^2abab^2ababab^2ababababa-b^2ab^2ab^2abab^2abab^2ababab^2ababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababab^2ababababababab^2abab^2ab-babab^2abab^2abab^2abab^2abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2abab^2abab^2ababab^2abababa-b^2ab^2abab^2abab^2ababab^2abababab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2abab^2ababab^2ababababab-b^2ab^2abab^2abab^2abab^2abab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababab^2ababab^2abab^2abab^2ab^2ab^2-babababab^2ababab^2abab^2abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababab^2ababab^2abab^2abab^2ab^2aba-b^2abababab^2ababab^2abab^2abab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2ab^2abab^2abab^2ababab^2ababab^2a-abab^2ab^2abab^2abab^2abab^2abab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2abab^2ababababababab^2abab^2ab^2-ab^2ab^2abab^2abab^2ababab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2ababababababab^2abab^2ab^2-b^2ab^2abab^2ababababababab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2abab^2abab^2abab^2abab^2ababab-b^2abababababababababababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababababababababababababababab^2-bababab^2abab^2abab^2abab^2abab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abab^2ababab^2ababab^2ababab^2abab^2a-abab^2abab^2abab^2ababab^2ababab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2ababab^2ababab^2ababab^2abab^2-bab^2abab^2abab^2ababab^2ababab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2ababab^2ababababababab^2abab-babab^2abab^2ababab^2ababababababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2abab^2ababab^2ababab^2ababa-b^2abab^2ababab^2ababab^2ababab^2abab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ababababababababababababababababa-ababab^2ababab^2ababab^2ababab^2ababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bab^2abab^2abab^2abab^2abab^2abab^2ababab-ab^2ab^2abab^2abab^2abab^2abab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2abab^2abab^2abab^2abab^2abab^2abab-ab^2abab^2abab^2abab^2abab^2abab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bababab^2abab^2abab^2abab^2abab^2abab^2ab-abab^2abab^2abab^2abab^2abab^2abab^2ab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "bababab^2ababab^2ababab^2ababab^2ababab-ababababababababababababababababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2abab^2abab^2abab^2abab^2ababa-b^2ab^2abab^2abab^2abab^2abab^2abab^2abab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2abab^2abab^2abab^2abab^2abab^2aba-b^2abab^2abab^2abab^2abab^2abab^2abab^2ab^2"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ab^2abab^2abab^2abab^2abab^2abab^2abab^2a-abab^2abab^2abab^2abab^2abab^2abab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2abab^2abab^2abab^2abab^2abab^2abab^2ab^2a-ababab^2abab^2abab^2abab^2abab^2abab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abab^2abab^2abab^2abab^2abab^2abab^2-bab^2abab^2abab^2abab^2abab^2abab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ab^2abab^2abab^2abab^2abab^2abab^2ababa-bab^2abab^2ababab^2ababab^2ababab^2ababab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2abab^2abab^2abab^2abab^2abab^2ab^2-babab^2abab^2abab^2abab^2abab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2abab^2ababab^2ababab^2ababab^2ababab^2-bab^2ab^2abab^2abab^2ababab^2ababab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ab^2ababab^2ababab^2ababab^2ababab^2abab^2-bab^2abab^2abab^2ababab^2ababab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2abab^2abab^2ababab^2ababab^2abab-b^2abab^2ababab^2ababab^2ababab^2ababab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2abab^2ababab^2ababab^2abab^2ab-b^2ababab^2ababab^2ababab^2ababab^2abab^2a"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2ababab^2ababab^2ababab^2ababa-bab^2ab^2abab^2abab^2abab^2abab^2abab^2abab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2abab^2abab^2abab^2abab^2abab^2ab^2a-bababab^2ababab^2ababab^2ababab^2abab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "ababab^2ababab^2ababab^2ababab^2abab^2aba-babab^2abab^2abab^2abab^2abab^2abab^2ab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "b^2ababab^2ababab^2ababab^2ababab^2abab^2ab-abab^2abab^2abab^2abab^2abab^2abab^2ab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "babab^2abab^2abab^2abab^2abab^2abab^2ab^2aba-ababab^2ababab^2ababab^2ababab^2abab^2ab"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2ab^2abab^2abab^2abab^2abab^2abab^2abab-bab^2abab^2ababab^2ababab^2ababab^2ababa"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abab^2abab^2abab^2abab^2abab^2abab^2abab^2ab-bab^2abab^2abab^2abab^2abab^2abab^2abab^2aba"),
				parseMust(map[string]Symbol{"a": 2, "b": 1}, Deglex, "abababababababababababababababababababab-babababababababababababababababababababa"),
			},
			complete: true,
			long:     true,
		},
	}

	for i, test := range tests {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			t.Parallel()
			if testing.Short() && test.long {
				t.Skip()
			}
			ideal := make([]*Polynomial, len(test.ideal))
			copy(ideal, test.ideal)

			basis, complete := Buchberger(ideal, test.maxiter)
			if len(basis) != len(test.basis) {
				t.Fatalf("%d %v", len(basis), basis)
			}
			for i := range basis {
				if basis[i].Cmp(test.basis[i]) != 0 {
					t.Errorf("%d %v", i, basis[i])
				}
			}
			if complete != test.complete {
				t.Errorf("got %v want %v", complete, test.complete)
			}
		})
	}
}

func TestAddObstructions(t *testing.T) {
	tests := []struct {
		b   []obstruction
		g   []*Polynomial
		obs []obstruction
	}{
		// OBS(2), Example 5.12, Mora.
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
		// OBS(4), Example 5.12, Mora.
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
		// Example 4.15, Clemens Hofstadler.
		{
			b: []obstruction{
				{i: 0, j: 1, iLeft: Monomial{}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2, 1}},
				{i: 0, j: 1, iLeft: Monomial{}, iRight: Monomial{}, jLeft: Monomial{2}, jRight: Monomial{1}},
				{i: 2, j: 3, iLeft: Monomial{1, 2}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2}},
				{i: 2, j: 3, iLeft: Monomial{}, iRight: Monomial{1}, jLeft: Monomial{}, jRight: Monomial{}},
			},
			g: []*Polynomial{
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 2, 1}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 2}},
					PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 2, 1}},
					PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 1}},
				),
				NewPolynomial(
					Deglex,
					PolynomialTerm{Coefficient: big.NewRat(2, 1), Monomial: Monomial{1}},
				),
			},
			obs: []obstruction{
				{i: 2, j: 3, iLeft: Monomial{1, 2}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2}},
				{i: 2, j: 3, iLeft: Monomial{}, iRight: Monomial{1}, jLeft: Monomial{}, jRight: Monomial{}},
				{i: 0, j: 4, iLeft: Monomial{}, iRight: Monomial{}, jLeft: Monomial{2, 2}, jRight: Monomial{}},
				{i: 2, j: 4, iLeft: Monomial{}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2}},
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
		// OBS(2), Example 5.12, Mora.
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
		// OBS(4), Example 5.12, Mora.
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
		// Example 3.13, Kreuzer.
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
		// Example 3.10, Kreuzer.
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
		obs     []obstruction
		order   Order
		removed []obstruction
	}{
		// OBS(2) Example 5.10, Mora.
		{
			obs: []obstruction{
				{i: 2, j: 2, iLeft: Monomial{2, 1}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{1, 2}},
				{i: 1, j: 2, iLeft: Monomial{}, iRight: Monomial{2}, jLeft: Monomial{1}, jRight: Monomial{}},
				{i: 1, j: 2, iLeft: Monomial{2}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{1}},
			},
			order: Deglex,
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
			order: Deglex,
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
			order: Deglex,
			removed: []obstruction{
				{i: 1, j: 5, iLeft: Monomial{1}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{1}},
				{i: 3, j: 5, iLeft: Monomial{1, 1}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2}},
				{i: 4, j: 5, iLeft: Monomial{}, iRight: Monomial{1, 2}, jLeft: Monomial{2}, jRight: Monomial{}},
			},
		},
		// Example 3.10, Kreuzer.
		{
			obs: []obstruction{
				{i: 1, j: 2, iLeft: Monomial{1, 2}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{}},
				{i: 1, j: 2, iLeft: Monomial{}, iRight: Monomial{1, 2}, jLeft: Monomial{}, jRight: Monomial{}},
				{i: 2, j: 2, iLeft: Monomial{}, iRight: Monomial{1, 2}, jLeft: Monomial{1, 2}, jRight: Monomial{}},
			},
			order: Deglex,
			removed: []obstruction{
				{i: 1, j: 2, iLeft: Monomial{}, iRight: Monomial{1, 2}, jLeft: Monomial{}, jRight: Monomial{}},
			},
		},
		{
			obs: []obstruction{
				{i: 1, j: 1, iLeft: Monomial{1, 2, 3}, iRight: Monomial{}, jLeft: Monomial{1, 2}, jRight: Monomial{2, 1}},
				{i: 1, j: 1, iLeft: Monomial{1, 2}, iRight: Monomial{}, jLeft: Monomial{1, 2}, jRight: Monomial{2, 1}},
			},
			order: Deglex,
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

			remove4b(obs, test.order)
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
		// Example 4.2.3, Xiu Xingqiang.
		{
			o:      obstruction{i: 1, j: 2, iLeft: Monomial{1, 2, 1}, iRight: Monomial{}, jLeft: Monomial{1}, jRight: Monomial{1, 1, 2, 1, 2}},
			shrunk: obstruction{i: 1, j: 2, iLeft: Monomial{2, 1}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{1, 1, 2, 1, 2}},
		},
		// Example 3.5, Kreuzer.
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
		// Example 4.2.3, Xiu Xingqiang.
		{
			obstruction: obstruction{i: 1, j: 3, iLeft: Monomial{1, 2, 1}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{1, 2}},
			im:          Monomial{1, 2, 1, 2},
			jm:          Monomial{1, 2, 1, 1, 2},
			overlap:     true,
		},
		// Example 4.2.3, Xiu Xingqiang.
		{
			obstruction: obstruction{i: 1, j: 2, iLeft: Monomial{1, 2, 1}, iRight: Monomial{}, jLeft: Monomial{1}, jRight: Monomial{1, 1, 2, 1, 2}},
			im:          Monomial{1, 2, 1, 2},
			jm:          Monomial{2},
			overlap:     false,
		},
		// Example 4.2.5, Xiu Xingqiang.
		{
			obstruction: obstruction{i: 2, j: 3, iLeft: Monomial{1, 2, 1, 1}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2, 2}},
			im:          Monomial{2, 2, 2},
			jm:          Monomial{1, 2, 1, 1, 2},
			overlap:     true,
		},
		// Example 4.2.5, Xiu Xingqiang.
		{
			obstruction: obstruction{i: 1, j: 2, iLeft: Monomial{1, 2}, iRight: Monomial{2}, jLeft: Monomial{1, 2, 1, 1}, jRight: Monomial{}},
			im:          Monomial{1, 1, 2, 2},
			jm:          Monomial{2, 2, 2},
			overlap:     true,
		},
		// Case d, example 4.1.8, Xiu Xingqiang.
		{
			obstruction: obstruction{i: 1, j: 2, iLeft: Monomial{}, iRight: Monomial{5, 6, 5, 6, 5, 6, 2}, jLeft: Monomial{4, 6, 5, 6, 5, 6, 5}, jRight: Monomial{}},
			im:          Monomial{4, 6, 5, 6, 5, 6, 5},
			jm:          Monomial{5, 6, 5, 6, 5, 6, 2},
			overlap:     false,
		},
		// Case e, example 4.1.8, Xiu Xingqiang.
		{
			obstruction: obstruction{i: 1, j: 2, iLeft: Monomial{}, iRight: Monomial{6, 5, 6, 5, 6, 5, 6, 2}, jLeft: Monomial{4, 6, 5, 6, 5, 6, 5, 6}, jRight: Monomial{}},
			im:          Monomial{4, 6, 5, 6, 5, 6, 5},
			jm:          Monomial{5, 6, 5, 6, 5, 6, 2},
			overlap:     false,
		},
		// Case f, example 4.1.8, Xiu Xingqiang.
		{
			obstruction: obstruction{i: 1, j: 2, iLeft: Monomial{}, iRight: Monomial{6, 5, 6, 5, 6, 5, 6, 5, 6, 2}, jLeft: Monomial{4, 6, 5, 6, 5, 6, 5, 6, 5, 6}, jRight: Monomial{}},
			im:          Monomial{4, 6, 5, 6, 5, 6, 5},
			jm:          Monomial{5, 6, 5, 6, 5, 6, 2},
			overlap:     false,
		},
		// Case a, Example 2.9, Kreuzer.
		{
			obstruction: obstruction{i: 1, j: 1, iLeft: Monomial{}, iRight: Monomial{1, 1, 1}, jLeft: Monomial{1, 1, 1}, jRight: Monomial{}},
			im:          Monomial{1, 1},
			jm:          Monomial{1, 1},
			overlap:     false,
		},
		// Case a, Example 2.9, Kreuzer.
		{
			obstruction: obstruction{i: 1, j: 1, iLeft: Monomial{}, iRight: Monomial{2, 1, 1}, jLeft: Monomial{1, 1, 2}, jRight: Monomial{}},
			im:          Monomial{1, 1},
			jm:          Monomial{1, 1},
			overlap:     false,
		},
		// Case b, Example 2.9, Kreuzer.
		{
			obstruction: obstruction{i: 1, j: 2, iLeft: Monomial{}, iRight: Monomial{1, 1, 2}, jLeft: Monomial{1, 1, 1}, jRight: Monomial{}},
			im:          Monomial{1, 1},
			jm:          Monomial{1, 2},
			overlap:     false,
		},
		// Case b, Example 2.9, Kreuzer.
		{
			obstruction: obstruction{i: 1, j: 2, iLeft: Monomial{1, 2, 2}, iRight: Monomial{}, jLeft: Monomial{}, jRight: Monomial{2, 1, 1}},
			im:          Monomial{1, 1},
			jm:          Monomial{1, 2},
			overlap:     false,
		},
		// Case c, Example 2.9, Kreuzer.
		{
			obstruction: obstruction{i: 2, j: 2, iLeft: Monomial{}, iRight: Monomial{1, 1, 2}, jLeft: Monomial{1, 2, 1}, jRight: Monomial{}},
			im:          Monomial{1, 2},
			jm:          Monomial{1, 2},
			overlap:     false,
		},
		// Example 3.10, Kreuzer.
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
		// Example 4.1.8, Xiu Xingqiang.
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
		// Example 3.2.2, Xiu Xingqiang.
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
		// Example 3.2.5, Xiu Xingqiang.
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
		// Section 6.4 Simplifying Polynomial Expressions, NCAlgebra.
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
		// Section 6.5 Poynomials and Rules, NCAlgebra.
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

			// Check if quotient * g + remainder == f.
			f := NewPolynomial(test.f.order)
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

func TestPolynomialCmp(t *testing.T) {
	tests := []struct {
		x *Polynomial
		y *Polynomial
		c int
	}{
		{
			x: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 1, 2}},
				PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2}},
			),
			y: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 1}},
				PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1, 2}},
			),
			c: 1,
		},
		{
			x: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(9, 1), Monomial: Monomial{2}},
				PolynomialTerm{Coefficient: big.NewRat(7, 1), Monomial: Monomial{}},
			),
			y: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(3, 1), Monomial: Monomial{1, 2}},
				PolynomialTerm{Coefficient: big.NewRat(2, 1), Monomial: Monomial{2, 1}},
			),
			c: -1,
		},
	}

	for i, test := range tests {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			t.Parallel()
			if c := test.x.Cmp(test.y); c != test.c {
				t.Errorf("%d", c)
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
			z := NewPolynomial(test.z.order)
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
			z := NewPolynomial(test.z.order).Mul(test.x, test.y)
			if z.Cmp(test.z) != 0 {
				t.Errorf("%v", z)
			}
		})
	}
}

func TestPow(t *testing.T) {
	tests := []struct {
		x *Polynomial
		y int
		z *Polynomial
	}{
		{
			x: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1}},
				PolynomialTerm{Coefficient: big.NewRat(-2, 1), Monomial: Monomial{2}},
			),
			y: 2,
			z: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(4, 1), Monomial: Monomial{2, 2}},
				PolynomialTerm{Coefficient: big.NewRat(-2, 1), Monomial: Monomial{2, 1}},
				PolynomialTerm{Coefficient: big.NewRat(-2, 1), Monomial: Monomial{1, 2}},
				PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 1}},
			),
		},
		{
			x: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1}},
				PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2}},
			),
			y: 3,
			z: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2, 2, 2}},
				PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 2, 1}},
				PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{2, 1, 2}},
				PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{2, 1, 1}},
				PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 2, 2}},
				PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1, 2, 1}},
				PolynomialTerm{Coefficient: big.NewRat(-1, 1), Monomial: Monomial{1, 1, 2}},
				PolynomialTerm{Coefficient: big.NewRat(1, 1), Monomial: Monomial{1, 1, 1}},
			),
		},
	}

	for i, test := range tests {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			t.Parallel()
			z := NewPolynomial(Deglex)
			z.Pow(test.x, test.y)
			if z.Cmp(test.z) != 0 {
				t.Errorf("%v", z)
			}
		})
	}
}

func TestMulScalar(t *testing.T) {
	tests := []struct {
		x      *Polynomial
		scalar *big.Rat
		z      *Polynomial
	}{
		{
			x: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(-2, 3), Monomial: Monomial{2, 1, 2}},
				PolynomialTerm{Coefficient: big.NewRat(5, 1), Monomial: Monomial{2}},
			),
			scalar: big.NewRat(-6, 1),
			z: NewPolynomial(
				Deglex,
				PolynomialTerm{Coefficient: big.NewRat(4, 1), Monomial: Monomial{2, 1, 2}},
				PolynomialTerm{Coefficient: big.NewRat(-30, 1), Monomial: Monomial{2}},
			),
		},
	}
	for i, test := range tests {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			t.Parallel()
			z := NewPolynomial(test.z.order)
			z.mulScalar(test.scalar, test.x)
			if z.Cmp(test.z) != 0 {
				t.Errorf("%v", z)
			}

			z = NewPolynomial(Deglex).Set(test.x)
			z.mulScalar(test.scalar, z)
			if z.Cmp(test.z) != 0 {
				t.Errorf("%v", z)
			}
		})
	}
}

func TestHomogeneous(t *testing.T) {
	tests := []struct {
		f           *Polynomial
		homogeneous bool
	}{
		{
			f:           parseMust(map[string]Symbol{"a": 1, "b": 2}, Deglex, "ab"),
			homogeneous: true,
		},
		{
			f:           parseMust(map[string]Symbol{"a": 1, "b": 2}, Deglex, "ab-a"),
			homogeneous: false,
		},
		{
			f:           parseMust(map[string]Symbol{"a": 1, "b": 2}, Deglex, "(a-b)a"),
			homogeneous: true,
		},
	}

	for i, test := range tests {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			t.Parallel()
			if h := homogeneous(test.f); h != test.homogeneous {
				t.Errorf("got %v want %v", h, test.homogeneous)
			}
		})
	}
}

func TestOrder(t *testing.T) {
	tests := []struct {
		words  []Monomial
		order  Order
		sorted []Monomial
	}{
		{
			words:  []Monomial{{2}, {1, 1, 2}, {1, 2}},
			order:  lexicographic,
			sorted: []Monomial{{1, 1, 2}, {1, 2}, {2}},
		},
		{
			words:  []Monomial{{1, 3, 4}, {1, 2, 4}, {1, 3, 6}, {3, 4, 6}, {2, 5, 6}, {1, 5, 6}, {2, 3, 4}, {3, 4, 5}, {1, 4, 5}, {1, 4, 6}, {2, 3, 6}, {2, 4, 5}, {1, 3, 5}, {1, 2, 5}, {2, 3, 5}, {1, 2, 3}, {4, 5, 6}, {3, 5, 6}, {2, 4, 6}, {1, 2, 6}},
			order:  lexicographic,
			sorted: []Monomial{{1, 2, 3}, {1, 2, 4}, {1, 2, 5}, {1, 2, 6}, {1, 3, 4}, {1, 3, 5}, {1, 3, 6}, {1, 4, 5}, {1, 4, 6}, {1, 5, 6}, {2, 3, 4}, {2, 3, 5}, {2, 3, 6}, {2, 4, 5}, {2, 4, 6}, {2, 5, 6}, {3, 4, 5}, {3, 4, 6}, {3, 5, 6}, {4, 5, 6}},
		},
		{
			words:  []Monomial{{2, 3}, {1, 3}, {1, 1}, {3, 3}, {1, 2}, {2, 2}},
			order:  lexicographic,
			sorted: []Monomial{{1, 1}, {1, 2}, {1, 3}, {2, 2}, {2, 3}, {3, 3}},
		},
		{
			words:  []Monomial{{'b', 'e', 'n', 'e', 'f', 'i', 't'}, {'b', 'e'}, {'b', 'e', 'n', 't'}, {'b', 'a', 'r', 'n', 'a', 'c', 'l', 'e'}, {'b', 'e', 'e', 'n'}},
			order:  lexicographic,
			sorted: []Monomial{{'b', 'a', 'r', 'n', 'a', 'c', 'l', 'e'}, {'b', 'e'}, {'b', 'e', 'e', 'n'}, {'b', 'e', 'n', 'e', 'f', 'i', 't'}, {'b', 'e', 'n', 't'}},
		},
		{
			words:  []Monomial{{1, 1}, {1, 1, 2}, {2, 1, 1}, {2, 2}, {1, 2, 2}, {2, 2, 2}, {}, {1, 2}, {1, 1, 1}, {2, 1, 2}, {2}, {1}, {2, 2, 1}, {1, 2, 1}, {2, 1}},
			order:  Deglex,
			sorted: []Monomial{{}, {1}, {2}, {1, 1}, {1, 2}, {2, 1}, {2, 2}, {1, 1, 1}, {1, 1, 2}, {1, 2, 1}, {1, 2, 2}, {2, 1, 1}, {2, 1, 2}, {2, 2, 1}, {2, 2, 2}},
		},
	}
	for i, test := range tests {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			t.Parallel()
			slices.SortFunc(test.words, test.order)
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
		z = NewPolynomial(z.order).Mul(z, y[i])
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

func parseMust(variables map[string]Symbol, order Order, input string) *Polynomial {
	p, err := Parse(variables, order, input)
	if err != nil {
		panic(fmt.Sprintf("%+v", err))
	}
	return p
}
