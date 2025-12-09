package scan

import (
	"fmt"
	"io"
	"strings"
	"unicode"
	"unicode/utf8"
)

type Type int

const (
	EOF Type = iota
	Error
	Int
	Identifier
	Operator
	Parenthesis
)

type Location struct {
	Line   int
	Column int
}

type Token struct {
	Type     Type
	Text     string
	Location Location
}

const eof = -1

type stateFn func(*Scanner) stateFn

type Scanner struct {
	token Token

	r         io.ByteReader
	input     string   // the line of text being scanned
	start     int      // starting position of this token
	pos       int      // current position of this token
	loc       Location // location of starting position
	done      bool
	lastWidth int
	buf       []byte
}

func NewScanner(r io.ByteReader) *Scanner {
	l := &Scanner{r: r}
	return l
}

func (l *Scanner) Next() Token {
	state := lexAny
	for {
		state = state(l)
		if state == nil {
			return l.token
		}
	}
}

func (l *Scanner) loadLine() {
	l.buf = l.buf[:0]
	for {
		c, err := l.r.ReadByte()
		if err != nil {
			l.done = true
			break
		}
		if c == '\r' {
			continue
		}

		l.buf = append(l.buf, c)
		if c == '\n' {
			break
		}
	}

	if l.start == l.pos {
		l.input = string(l.buf)
		l.start, l.pos = 0, 0
	} else {
		l.input += string(l.buf)
	}
}

func (l *Scanner) readRune() (rune, int) {
	if !l.done && l.pos == len(l.input) {
		l.loadLine()
	}
	if l.pos == len(l.input) {
		return eof, 0
	}
	return utf8.DecodeRuneInString(l.input[l.pos:])
}

func (l *Scanner) peek() rune {
	r, _ := l.readRune()
	return r
}

func (l *Scanner) next() rune {
	var r rune
	r, l.lastWidth = l.readRune()
	l.pos += l.lastWidth
	return r
}

func (l *Scanner) backup() {
	l.pos -= l.lastWidth
}

func (l *Scanner) acceptRun(valid string) {
	for strings.ContainsRune(valid, l.next()) {
	}
	l.backup()
}

func (l *Scanner) emit(t Type) stateFn {
	text := l.input[l.start:l.pos]
	l.token = Token{Type: t, Text: text, Location: l.loc}

	// Update starting position.
	for _, c := range text {
		switch c {
		case '\n':
			l.loc.Line++
			l.loc.Column = 0
		default:
			l.loc.Column++
		}
	}
	l.start = l.pos

	return nil
}

func (l *Scanner) errorf(format string, args ...interface{}) stateFn {
	l.token = Token{
		Type:     Error,
		Text:     fmt.Sprintf(format, args...),
		Location: l.loc,
	}
	l.input = l.input[:0]
	l.start, l.pos = 0, 0
	return nil
}

func lexAny(l *Scanner) stateFn {
	switch r := l.next(); {
	case r == eof:
		l.token = Token{Type: EOF, Text: "EOF"}
		return nil
	case r == '\n':
		l.start = l.pos
		return lexAny
	case isSpace(r):
		return lexSpace
	case strings.ContainsRune("()", r):
		return l.emit(Parenthesis)
	case strings.ContainsRune("+-/^", r):
		return l.emit(Operator)
	case '0' <= r && r <= '9':
		return lexInt
	case r == '{':
		l.backup()
		return lexBracket
	case unicode.IsLetter(r):
		return l.emit(Identifier)
	default:
		return l.errorf("unrecognized character: %#U", r)
	}
}

func lexBracket(l *Scanner) stateFn {
	var stack int
	for {
		r := l.next()
		switch r {
		case '{':
			stack++
		case '}':
			stack--
			if stack == 0 {
				return l.emit(Identifier)
			}
		}
	}
}

func lexInt(l *Scanner) stateFn {
	const digits = "0123456789"
	l.acceptRun(digits)
	return l.emit(Int)
}

func lexSpace(l *Scanner) stateFn {
	for isSpace(l.peek()) {
		l.next()
	}
	l.start = l.pos
	return lexAny
}

func isSpace(r rune) bool {
	return r == ' ' || r == '\t'
}
