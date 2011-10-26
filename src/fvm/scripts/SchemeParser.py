## adapted from pyscheme's parser code


import re

from UserString import UserString

class Symbol(UserString): pass

def pogo(bouncer):
    """A trampoline that bounces a single bouncer.  See:

    http://www.cs.indiana.edu/hyplan/sganz/publications/icfp99/paper.pdf
    """
    try:
        while True:
            if bouncer[0] == 'land':
                return bouncer[1]
            elif bouncer[0] == 'bounce':
                bouncer = bouncer[1](*bouncer[2])
            else:
                traceback.print_exc()
                raise TypeError, "not a bouncer"
    except TypeError:
        traceback.print_exc()
        raise TypeError, "not a bouncer"


def bounce(function, *args):
    """Returns a new trampolined value that continues bouncing on the
    trampoline."""
    return ('bounce', function, args)


def land(value):
    """Returns a new trampolined value that lands off the trampoline."""
    return ('land', value)



"""The End Of File token is a sentinal that terminates a list of tokens."""
EOF_TOKEN = (None, None)


"""Here are the patterns we pay attention when breaking down a string
into tokens."""
PATTERNS = [ ('whitespace', re.compile(r'(\s+)')),
             ('comment', re.compile(r'(;[^\n]*)')),
             ('(', re.compile(r'(\()')),
             (')', re.compile(r'(\))')),
             ('number', re.compile(r'''( [+\-]?    ## optional sign,
                                         (?:       ## followed by some
                                                   ## decimals
                                            \d+\.?\d*[eE]?[+\-]?\d+
                                             | \d+\.\d+
                                            | \d+\.
                                            | \.\d+
                                            | \d+
#                                            |\d+e-\d+
                                            )
                                       )
                                         ''',
                                   re.VERBOSE)),
             ('symbol',
              re.compile(r'''([_a-zA-Z:\+\=\?\!\@\#\$\%\^\&\*\-\/\.\>\<]
                              [\w\+\=\?\!\@\#\$\%\^\&\*\-\/\.\>\<]*)''',
                                   re.VERBOSE)),
             ('string', re.compile(r'''
                                      "
                                      (([^\"] | \\")*)
                                      "
                                      ''',
                                   re.VERBOSE)),
             ]

broken_number_re1 = re.compile( r'(\d\.)(\s|\))')

def tokenize(s_in):
    """Given a string 's', return a list of its tokens.  A token can
    be one of the following types listed in the PATTERNS above, and
    each token is a 2-tuple: (type, content)."""

    ## fix numbers like 4. to 4.0 so that the tokenizer works correctly
    s = broken_number_re1.sub(r'\g<1>0\g<2>',s_in)

    tokens = []
    while 1:
        should_continue = 0
        for tokenType, regex in PATTERNS:
            match_obj = regex.match(s)
            if match_obj:
                should_continue = 1
                tokens.append( (tokenType, match_obj.group(1)) )
                s = s[match_obj.span()[1] :]
        if should_continue == 0:
            break
    tokens.append(EOF_TOKEN)
    return filter(lambda x: x[0] not in ('whitespace', 'comment'),
                  tokens)


class ParserError(Exception):
    """Our personalized exception class."""
    pass


def peek(tokens):
    """Take a quick glance at the first token in our tokens list."""
    if len(tokens) == 0:
        raise ParserError, "While peeking: ran out of tokens."
    return tokens[0]


def eat(tokens, tokenType):
    """Digest the first token in our tokens list, making sure that we're
    biting on the right tokenType of thing."""
    if len(tokens) == 0:
        raise ParserError, "While trying to eat %s: ran out of tokens." % \
              (repr(tokenType),)
    if tokens[0][0] != tokenType:
        raise ParserError, "While trying to eat %s: got token %s instead." % \
                            (repr(tokenType), repr(tokens[0]))
    return tokens.pop(0)


def identity_cont(val):
    return land(val)


def parseSingleExpression(tokens, cont=identity_cont):
    """Returns a single Expression, given a sequence of tokens.
    Raises a ParserException if our tokens haven't been exhausted."""
    def c(expression):
        eat(tokens, None)
        return bounce(cont, expression)

    return parseExpression(tokens, c)


def parseExpression(tokens, cont=identity_cont):
    """Returns an Expression, given a sequence of tokens.
    An expression is made up of one of the following things:
        o  A quoted expression
        o  An atom (like a number or symbol or string)
        o  A list.
    This procedure tries to take care of all these potentials."""
    look_ahead_type = peek(tokens)[0]

    if look_ahead_type == '(':
        return parseList(tokens, cont)
    elif look_ahead_type in ('number', 'symbol', 'string'):
        return parseAtom(tokens, cont)
    else:
        raise ParserError, "While parsing Expression: no alternatives."


def parseAtom(tokens, cont):
    """Returns an Atom, given a sequence of tokens.
    An atom is either a number, a symbol, or a string."""
    if peek(tokens)[0] == 'number':
        return bounce(cont, toNumber(eat(tokens, 'number')[1]))
    elif peek(tokens)[0] == 'symbol':
        return bounce(cont, Symbol(eat(tokens, 'symbol')[1]))
    elif peek(tokens)[0] == 'string':
        return bounce(cont, eat(tokens, 'string')[1])
    else:
        raise ParserError, "While parsing Atom: no alternatives."


def toNumber(s):
    """Tries to convert string 's' into a number."""
    try:
        return int(s)
    except ValueError:
        return float(s)

NIL = []


class ConsPair(list): pass


def cons(head, rest):
    """Returns the concatentation of head with rest."""
    return ConsPair([head, rest])


def parseList(tokens, cont):
    """Parses a parenthesized list expression."""
    eat(tokens, "(")
    def c_expressionsEaten(val):
        eat(tokens, ")")
        return bounce(cont, val)
    return bounce(parseExpressionStar, tokens, c_expressionsEaten)


def parseExpressionStar(tokens, cont):
    """Tries to eat as many expressions as it can see."""
    START_EXPR_SET = ( '(', 'number', 'symbol', 'string')
    if peek(tokens) == ('symbol', '.'):
        ## Handle dotted pairs
        eat(tokens, "symbol")
        return bounce(parseExpression, tokens, cont)
    elif peek(tokens)[0] in START_EXPR_SET:
        def c_first_eaten(firstVal):
            def c_rest_eaten(restVal):
                return bounce(cont, cons(firstVal, restVal))
            return bounce(parseExpressionStar, tokens, c_rest_eaten)
        return bounce(parseExpression, tokens, c_first_eaten)
    else:
        return bounce(cont, NIL)


def parse(s):
    """Parse a single string.  This is just a convenience function."""
    return pogo(parseSingleExpression(tokenize(s),
                                      identity_cont))


