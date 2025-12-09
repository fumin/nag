package parse

import (
	"bytes"
	"flag"
	"fmt"
	"log"
	"testing"

	"github.com/fumin/nag/parse/scan"
)

func TestParse(t *testing.T) {
	tests := []struct {
		input string
		tree  string
	}{
		{
			input: "ba^3",
			tree:  "(b*(a^3))",
		},
		{
			input: "-ba^3",
			tree:  "(0-(b*(a^3)))",
		},
		{
			input: "(a+b)^4",
			tree:  "((a+b)^4)",
		},
		{
			input: "-12/5a^3((a+cc)b)^2a+7/3ca-3/2b",
			tree:  "(((0-((((12/5)*(a^3))*(((a+(c*c))*b)^2))*a))+(((7/3)*c)*a))-((3/2)*b))",
		},
		{
			input: "5/3b(a+b)^2c+9a",
			tree:  "(((((5/3)*b)*((a+b)^2))*c)+(9*a))",
		},
	}

	for i, test := range tests {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			t.Parallel()
			n, err := Parse(scan.NewScanner(bytes.NewBufferString(test.input)))
			if err != nil {
				t.Fatalf("%+v", err)
			}
			if tree(n) != test.tree {
				t.Errorf("%s", tree(n))
			}
		})
	}
}

func TestMain(m *testing.M) {
	flag.Parse()
	log.SetFlags(log.Lmicroseconds | log.Llongfile | log.LstdFlags)

	m.Run()
}
