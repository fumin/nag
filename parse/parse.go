package parse

import (
	"fmt"

	"github.com/pkg/errors"

	"github.com/fumin/nag/parse/scan"
)

const (
	// AddedLine represents a line that does not exist in the source input, but which is specifically added by the parser.
	AddedLine = -1
)

type Node struct {
	Token  scan.Token
	Parent *Node
	Left   *Node
	Right  *Node
}

func Parse(scanner *scan.Scanner) (*Node, error) {
	// rightMost is the right most node in the current stack.
	// It is either:
	//   * A parenthesis node.
	//   * An operator node whose right child is nil.
	//   * An identifier node.
	rightMost, err := parseFirstToken(scanner)
	if err != nil {
		return nil, errors.Wrap(err, "")
	}

	for {
		tok := scanner.Next()
		var err error
		switch tok.Type {
		case scan.EOF:
			return root(rightMost), nil
		case scan.Parenthesis:
			if tok.Text == ")" {
				return root(rightMost), nil
			}
			rightMost, err = parseParenthesis(rightMost, tok, scanner)
		case scan.Operator:
			rightMost, err = parseOperator(rightMost, tok, nil)
		case scan.Int:
			rightMost, err = parseIdentifier(rightMost, tok)
		case scan.Identifier:
			rightMost, err = parseIdentifier(rightMost, tok)
		default:
			err = errors.Errorf("%d: %s", tok.Location.Column, tok.Text)
		}
		if err != nil {
			return nil, err
		}
	}
}

func parseParenthesis(rightMost *Node, tok scan.Token, scanner *scan.Scanner) (*Node, error) {
	// Parse the expression enclosed within the parentheses.
	expr, err := Parse(scanner)
	if err != nil {
		return nil, errors.Wrap(err, "")
	}
	pNode := &Node{Token: tok}
	setLeft(pNode, expr)

	// Directly insert pNode when right most node is an operator whose right child is nil.
	if rightMost == nil {
		return pNode, nil
	}
	if rightMost.Token.Type == scan.Operator {
		setRight(rightMost, pNode)
		return pNode, nil
	}

	// Perform implicit multiplication when right most node is an identifier.
	mulTok := scan.Token{Type: scan.Operator, Text: "*", Location: scan.Location{Line: AddedLine}}
	if _, err := parseOperator(rightMost, mulTok, pNode); err != nil {
		return nil, errors.Wrap(err, "")
	}
	return pNode, nil
}

func parseOperator(rightMost *Node, tok scan.Token, rightChild *Node) (*Node, error) {
	// Find right most node according to operator precedence.
	for rightMost.Parent != nil {
		if opOrder(tok.Text) > opOrder(rightMost.Parent.Token.Text) {
			break
		}
		rightMost = rightMost.Parent
	}

	op := &Node{Token: tok}
	setRight(rightMost.Parent, op)
	setLeft(op, rightMost)
	setRight(op, rightChild)
	return op, nil
}

func parseIdentifier(rightMost *Node, tok scan.Token) (*Node, error) {
	// Directly insert iNode when right most node is an operator whose right child is nil.
	iNode := &Node{Token: tok}
	if rightMost.Token.Type == scan.Operator {
		setRight(rightMost, iNode)
		return iNode, nil
	}

	// Perform implicit multiplication when right most node is an identifier.
	mulTok := scan.Token{Type: scan.Operator, Text: "*", Location: scan.Location{Line: AddedLine}}
	if _, err := parseOperator(rightMost, mulTok, iNode); err != nil {
		return nil, errors.Wrap(err, "")
	}
	return iNode, nil
}

func parseFirstToken(scanner *scan.Scanner) (*Node, error) {
	tok := scanner.Next()
	switch tok.Type {
	case scan.Parenthesis:
		return parseParenthesis(nil, tok, scanner)
	case scan.Operator:
		rightMost := &Node{Token: tok}
		setLeft(rightMost, &Node{Token: scan.Token{Type: scan.Int, Text: "0", Location: scan.Location{Line: AddedLine}}})
		return rightMost, nil
	case scan.Int:
		fallthrough
	case scan.Identifier:
		return &Node{Token: tok}, nil
	default:
		return nil, errors.Errorf("unknown token %#v", tok)
	}
}

func tree(n *Node) string {
	switch n.Token.Type {
	case scan.Parenthesis:
		return tree(n.Left)
	case scan.Operator:
		return "(" + tree(n.Left) + n.Token.Text + tree(n.Right) + ")"
	case scan.Identifier:
		fallthrough
	case scan.Int:
		return n.Token.Text
	default:
		panic(fmt.Sprintf("%#v", n))
	}
}

func opOrder(op string) int {
	switch op {
	case "+":
		return 0
	case "-":
		return 0
	case "*":
		return 1
	case "/":
		return 2
	case "^":
		return 3
	default:
		panic(fmt.Sprintf("unknown operator \"%s\"", op))
	}
}

func root(n *Node) *Node {
	for n.Parent != nil {
		n = n.Parent
	}
	return n
}

func setLeft(parent, n *Node) {
	if parent == n {
		panic(fmt.Sprintf("%#v", n))
	}
	n.Parent = parent
	parent.Left = n
}

func setRight(parent, n *Node) {
	if parent == n {
		panic(fmt.Sprintf("%#v", n))
	}
	if n != nil {
		n.Parent = parent
	}
	if parent != nil {
		parent.Right = n
	}
}
