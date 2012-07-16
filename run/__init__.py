import builtins
import operator
import functools
import importlib

from .shell import dg

# Choose a function based on the number of arguments.
varary  = lambda *fs: lambda *xs: fs[len(xs) - 1](*xs)


builtins.__dict__.update({
    # Runtime counterparts of some stuff in `Compiler.builtins`.
    '$': lambda f, *xs: f(*xs)
  , ':': lambda f, *xs: f(*xs)
  , ',': lambda a, *xs: (a,) + xs

  , '<':  operator.lt
  , '<=': operator.le
  , '==': operator.eq
  , '!=': operator.ne
  , '>':  operator.gt
  , '>=': operator.ge
  , 'is': operator.is_
  , 'in': lambda a, b: a in b

  , 'not': operator.not_
  , '~':  operator.invert
  , '+':  varary(operator.pos, operator.add)
  , '-':  varary(operator.neg, operator.sub)
  , '*':  operator.mul
  , '**': operator.pow
  , '/':  operator.truediv
  , '//': operator.floordiv
  , '%':  operator.mod
  , '!!': operator.getitem
  , '&':  operator.and_
  , '^':  operator.xor
  , '|':  operator.or_
  , '<<': operator.lshift
  , '>>': operator.rshift

    # Useful stuff.
  , 'import': importlib.import_module
  , 'foldl': functools.reduce
  , '~:': functools.partial
  , '...': Ellipsis
})