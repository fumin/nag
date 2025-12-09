package scan

import (
	"bytes"
	"flag"
	"fmt"
	"log"
	"slices"
	"testing"
)

func TestScanner(t *testing.T) {
	tests := []struct {
		input  string
		tokens []Token
	}{
		{
			input: `101b + 71/2 {b^{\dagger}}^2 (a - bab)^3`,
			tokens: []Token{
				{Type: Int, Text: "101", Location: Location{Line: 0, Column: 0}},
				{Type: Identifier, Text: "b", Location: Location{Line: 0, Column: 3}},
				{Type: Operator, Text: "+", Location: Location{Line: 0, Column: 4}},
				{Type: Int, Text: "71", Location: Location{Line: 0, Column: 5}},
				{Type: Operator, Text: "/", Location: Location{Line: 0, Column: 7}},
				{Type: Int, Text: "2", Location: Location{Line: 0, Column: 8}},
				{Type: Identifier, Text: `{b^{\dagger}}`, Location: Location{Line: 0, Column: 9}},
				{Type: Operator, Text: "^", Location: Location{Line: 0, Column: 22}},
				{Type: Int, Text: "2", Location: Location{Line: 0, Column: 23}},
				{Type: Parenthesis, Text: "(", Location: Location{Line: 0, Column: 24}},
				{Type: Identifier, Text: "a", Location: Location{Line: 0, Column: 25}},
				{Type: Operator, Text: "-", Location: Location{Line: 0, Column: 26}},
				{Type: Identifier, Text: "b", Location: Location{Line: 0, Column: 27}},
				{Type: Identifier, Text: "a", Location: Location{Line: 0, Column: 28}},
				{Type: Identifier, Text: "b", Location: Location{Line: 0, Column: 29}},
				{Type: Parenthesis, Text: ")", Location: Location{Line: 0, Column: 30}},
				{Type: Operator, Text: "^", Location: Location{Line: 0, Column: 31}},
				{Type: Int, Text: "3", Location: Location{Line: 0, Column: 32}},
			},
		},
	}

	for i, test := range tests {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			t.Parallel()
			l := NewScanner(bytes.NewBufferString(test.input))
			var tokens []Token
			for i := l.Next(); i.Type != EOF; i = l.Next() {
				tokens = append(tokens, i)
			}
			if !slices.Equal(tokens, test.tokens) {
				t.Errorf("%v", tokens)
			}
		})
	}
}

func TestMain(m *testing.M) {
	flag.Parse()
	log.SetFlags(log.Lmicroseconds | log.Llongfile | log.LstdFlags)

	m.Run()
}
