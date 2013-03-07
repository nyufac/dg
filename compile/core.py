import dis
import types
import collections

from .. import parse

__all__ = ['INFIX_LEFT', 'INFIX_RIGHT', 'PREFIX', 'scanvars', 'compile', 'Code']


### VARIABLE SCANNER
#
# Instead of actually compiling code, this part looks for whatever variables
# it references.
#

# scanvars :: (StructMixIn, set Link) -> (set Link, set Link, set Link, set Link)
#
# Scan an AST for variable names. The returned sets contain:
#   1. global variables
#   2. local variables
#   3. cell variables
#   4. free variables
#
def scanvars(code, cell, nolocals, rightbind=False):

    if isinstance(code, parse.tree.Link):

        return (
          (set(), set(), set(), {code}) if code in cell else
          ({code}, set(), set(), set())
        )

    if isinstance(code, parse.tree.Constant):

        return set(), set(), set(), set()

    if isinstance(code, parse.tree.Expression):

        f, *args = code
        g =  INFIX_RIGHT[f] if f.infix and not f.closed and rightbind \
        else INFIX_LEFT[f]  if f.infix and not f.closed and len(args) == 1 \
        else PREFIX[f] if isinstance(f, parse.tree.Link) else Code.nativecall

        return getattr(g, 'scanvars', joinvars)(code, cell, nolocals)

    raise TypeError('not an AST output structure: {!r}'.format(code))


# joinvars :: ([StructMixIn], set Link) -> (set Link, set Link, set Link, set Link)
#
# Like `scanvars`, but for multiple items in the same expression.
#
def joinvars(code, cell, nolocals):

    globals = set()
    locals  = set()
    ncell   = set()
    free    = set()

    for item in code:

        a, b, c, d = scanvars(item, cell, nolocals)
        globals |= a
        globals -= b | c | d
        locals  |= b
        locals  -= c | d
        ncell   |= c
        ncell   -= d
        free    |= d

    return globals, locals, ncell, free


# joinvars1 :: ([StructMixIn], set Link) -> (set Link, set Link, set Link, set Link)
#
# Like `joinvars`, but ignores the first item, which is generally the function name.
# Useful for compile-time functions::
#
#     (PREFIX !! '+') = bind binary_op 'BINARY_ADD' 'INPLACE_ADD'
#     (PREFIX !! '+').scanvars = compile.joinvars1
#
def joinvars1(code, cell, nolocals):

    return joinvars(code[1:], cell, nolocals)


# joinvarsX :: (*StructMixIn, Maybe bool) -> typeof joinvars
#
# A wrapper for `joinvars` and `joinvars1` that adds a few more items.
# Useful when a compile-time function uses variables.
#
def joinvarsX(*add, f=joinvars):

    add = list(add)
    return lambda code, cell, nolocals: f(add + list(code), cell, nolocals)

### VARIABLE SCANNER--


# data Code = Code (str, int, int, int, set Link, set Link, set Link, set Link)
#           where bytecode :: [(int, int)]
#                 lineno   :: int
#                 depth    :: int
#                 depthmax :: int
#                 lnomap   :: dict int int
#
# A mutable `code` object.
#
class Code:

    OPTIMIZED = 1
    NEWLOCALS = 2
    VARARGS   = 4
    VARKWARGS = 8
    NESTED    = 16
    GENERATOR = 32
    NOFREE    = 64

    OPCODE_STACK_DELTA = {
     'NOP':          0,
     'EXTENDED_ARG': 0,

     'POP_TOP':    -1,
     'ROT_TWO':     0,
     'ROT_THREE':   0,
     'DUP_TOP':     1,
     'DUP_TOP_TWO': 2,

     'BINARY_ADD':           -1,
     'BINARY_AND':           -1,
     'BINARY_FLOOR_DIVIDE':  -1,
     'BINARY_LSHIFT':        -1,
     'BINARY_MODULO':        -1,
     'BINARY_MULTIPLY':      -1,
     'BINARY_OR':            -1,
     'BINARY_POWER':         -1,
     'BINARY_RSHIFT':        -1,
     'BINARY_SUBSCR':        -1,
     'BINARY_SUBTRACT':      -1,
     'BINARY_TRUE_DIVIDE':   -1,
     'BINARY_XOR':           -1,
     'COMPARE_OP':           -1,
     'INPLACE_ADD':          -1,
     'INPLACE_AND':          -1,
     'INPLACE_FLOOR_DIVIDE': -1,
     'INPLACE_LSHIFT':       -1,
     'INPLACE_MODULO':       -1,
     'INPLACE_MULTIPLY':     -1,
     'INPLACE_OR':           -1,
     'INPLACE_POWER':        -1,
     'INPLACE_RSHIFT':       -1,
     'INPLACE_SUBTRACT':     -1,
     'INPLACE_TRUE_DIVIDE':  -1,
     'INPLACE_XOR':          -1,

     'UNARY_INVERT':   0,
     'UNARY_NEGATIVE': 0,
     'UNARY_NOT':      0,
     'UNARY_POSITIVE': 0,

     'LOAD_ATTR':        0,
     'LOAD_BUILD_CLASS': 1,
     'LOAD_CLOSURE':     1,
     'LOAD_CONST':       1,
     'LOAD_DEREF':       1,
     'LOAD_FAST':        1,
     'LOAD_GLOBAL':      1,
     'LOAD_NAME':        1,

     'STORE_ATTR':   -2,
     'STORE_DEREF':  -1,
     'STORE_FAST':   -1,
     'STORE_GLOBAL': -1,
     'STORE_LOCALS': -1,
     'STORE_NAME':   -1,
     'STORE_SUBSCR': -3,

     'DELETE_ATTR':   -1,
     'DELETE_DEREF':   0,
     'DELETE_FAST':    0,
     'DELETE_GLOBAL':  0,
     'DELETE_NAME':    0,
     'DELETE_SUBSCR': -2,

     'FOR_ITER': 1,
     'GET_ITER': 0,

     'IMPORT_FROM':  1,
     'IMPORT_NAME': -1,
     'IMPORT_STAR': -1,

     'JUMP_ABSOLUTE':         0,
     'JUMP_FORWARD':          0,
     'JUMP_IF_FALSE_OR_POP': -1,
     'JUMP_IF_TRUE_OR_POP':  -1,
     'POP_JUMP_IF_FALSE':    -1,
     'POP_JUMP_IF_TRUE':     -1,

     'BUILD_MAP':    1,
     'LIST_APPEND': -1,
     'MAP_ADD':     -2,
     'SET_ADD':     -1,
     'STORE_MAP':   -2,

     'PRINT_EXPR':   -1,
     'RETURN_VALUE': -1,
     'YIELD_FROM':   -1,
     'YIELD_VALUE':  -1,

     'SETUP_EXCEPT':  0,
     'SETUP_FINALLY': 0,
     'SETUP_LOOP':    0,
     'SETUP_WITH':    2,
    }

    def __init__(self, name, argc, kwargc, flags, globals, locals, cell, free):

        super().__init__()
        self.name  = name
        self.flags = flags

        self.argc   = argc
        self.kwargc = kwargc

        self.filename = '<generated>'
        self.lineno   = 0
        self.depth    = 0
        self.depthmax = 0
        self.bytecode = []
        self.lnomap   = {}

        # Just in case.
        # Run with `-O` (or was it `-OO`?) to ignore these checks.
        assert len(locals) >= argc + kwargc + bool(flags & self.VARARGS) + bool(flags & self.VARKWARGS)
        assert bool(flags & self.NOFREE) == (not cell)
        assert bool(flags & self.NESTED) == bool(free)

        assert not (free & cell)
        assert not (free & locals)
        assert not (cell & locals)
        assert not (free & globals)
        assert not (cell & globals)
        assert not (locals & globals)

        self.consts   = {}
        self.names    = dict(zip(globals, range(len(globals))))
        self.varnames = dict(zip(locals,  range(len(locals))))
        self.cellvars = dict(zip(cell,    range(len(cell))))
        self.freevars = dict(zip(free,    range(len(cell), len(cell) + len(free))))

    # error :: (type < Exception, *, **) -> _|_
    #
    # An expression equivalent to `raise`.
    #
    def error(self, T, *args, **kwargs):

        raise T(*args, **kwargs)

    # appendcode :: (str, Maybe (Either int DelayedComputation), Maybe int) -> ()
    #
    # Append a single opcode given its name (sometimes with an argument.)
    #
    def appendcode(self, name, argument=None, stackdelta=None):

        self.depth += stackdelta if stackdelta is not None \
                 else self.OPCODE_STACK_DELTA[name] if name in self.OPCODE_STACK_DELTA \
                 else self.error(AssertionError, '{!r} requires manual delta calculation'.format(name))
        self.depthmax = max(self.depth, self.depthmax)
        self.bytecode.append((dis.opmap[name], argument))

    for opname in dis.opmap:

        locals()[opname] = lambda self, argument=None, stackdelta=None, name=opname: \
          self.appendcode(name, argument, stackdelta)

    # appendname :: str -> int
    #
    # Add something to the `names` map.
    #
    def appendname(self, name):

        assert isinstance(name, parse.tree.Link), 'not a valid name'
        v = self.names[name] = self.names.get(name, len(self.names))
        return v

    # bakedopcode :: (int, int) -> iter bytes
    #
    # Compile a single instruction.
    #
    def bakedopcode(self, code, arg):

        if code >= dis.HAVE_ARGUMENT and arg > 0xFFFF:

            yield from self.bakedopcode(dis.opmap['EXTENDED_ARG'], arg // 0xFFFF)

        yield code

        if code >= dis.HAVE_ARGUMENT:

            yield arg %  256
            yield arg // 256

    @property
    # bakedcode :: iter bytes
    #
    # Compiled bytecode string.
    #
    def bakedcode(self):

        for code, arg in self.bytecode:

            yield from self.bakedopcode(code, arg)

    @property
    # immutable :: code
    #
    # Compiled code object.
    #
    def immutable(self):

        return types.CodeType(
            self.argc,
            self.kwargc,
            len(self.varnames),
            self.depthmax,
            self.flags,
            bytes(self.bakedcode),
            tuple(x for _, x in sorted(self.consts,   key=self.consts.__getitem__)),
            tuple(              sorted(self.names,    key=self.names.__getitem__)),
            tuple(              sorted(self.varnames, key=self.varnames.__getitem__)),
            self.filename, self.name,
            self.lineno, b'',
            tuple(              sorted(self.freevars, key=self.freevars.__getitem__)),
            tuple(              sorted(self.cellvars, key=self.cellvars.__getitem__)),
        )

    # push :: StructMixIn -> ()
    #
    # Append a sequence of opcodes that pushes the result of an expression
    # onto the stack.
    #
    def push(self, item):

        mdepth = self.depth

        if isinstance(item, parse.tree.Link):

            self.LOAD_FAST (self.varnames[item]) if item in self.varnames else \
            self.LOAD_DEREF(self.freevars[item]) if item in self.freevars else \
            self.LOAD_DEREF(self.freevars[item]) if item in self.cellvars else \
            self.LOAD_NAME (self.names[item])    if item in self.names    else \
            self.error(AssertionError, 'variable scanner error')

        elif isinstance(item, parse.tree.Constant):

            val = item.value
            self.consts[type(val), val] = a = self.consts.get((type(val), val), len(self.consts))
            self.LOAD_CONST(a)

        else:

            self.call(*item)

        assert mdepth + 1 == self.depth, 'stack depth calculation error: {} -> {}'.format(mdepth, self.depth)

    # call :: *StructMixIn -> ()
    #
    # Call a (possibly built-in) function.
    #
    def call(self, f, *args, rightbind=False):

        return self.call(*args, rightbind=True) if f.infix and f == '' \
          else INFIX_RIGHT[f](self, f, *args) if f.infix and not f.closed and rightbind \
          else INFIX_LEFT[f](self, f, *args)  if f.infix and not f.closed and len(args) == 1 \
          else PREFIX[f](self, f, *args) if isinstance(f, parse.tree.Link) else self.nativecall(f, *args)

    # NOTE this is NOT the same as `joinvars1` because of `rightbind`.
    call.scanvars = lambda code, cell, nolocals: scanvars(parse.tree.Expression(code[1:]), cell, nolocals, rightbind=True)

    # infixbindl :: (StructMixIn, StructMixIn) -> ()
    #
    # Default implementation of a left infix bind.
    #
    def infixbindl(self, op, arg):

        self.push(parse.tree.Link('bind'))
        self.push(op)
        self.push(arg)
        self.CALL_FUNCTION(2, -2)

    infixbindl.scanvars = joinvarsX(parse.tree.Link('bind'))

    # infixbindr :: (StructMixIn, StructMixIn+) -> ()
    #
    # Default implementation of right infix bind.
    #
    def infixbindr(self, op, *args):

        self.push(parse.tree.Link('bind'))
        self.push(parse.tree.Link('flip'))
        self.push(op)
        self.CALL_FUNCTION(1, -1)
        self.push(*args) if len(args) == 1 else self.call(*args)
        self.CALL_FUNCTION(2, -2)

    infixbindr.scanvars = joinvarsX(parse.tree.Link('bind'), parse.tree.Link('flip'))

    # nativecall :: (StructMixIn, *StructMixIn) -> ()
    #
    # Call a non-builtin (i.e. defined in runtime) function.
    #
    def nativecall(self, f, *args):

        self.push(f)
        self.nativecall_args(args, 0, f.infix and not f.closed)

    # nativecall_args :: ([StructMixIn], int, bool) -> ()
    #
    # Call a non-builtin function that is already on top of the stack.
    #
    def nativecall_args(self, args, preloaded, infix=False):

        defs  = args, None, None, {}, (), ()
        args, _, _, kwargs, vararg, varkwarg = defs if infix else parse.syntax.argspec(args, definition=False)

        for arg in args:

            self.push(arg)

        for kw, value in kwargs.items():

            self.push(parse.tree.Constant(str(kw)))
            self.push(value)

        for item in (vararg + varkwarg):

            self.push(item)

        (
          self.CALL_FUNCTION_VAR_KW if vararg and varkwarg else
          self.CALL_FUNCTION_VAR    if vararg else
          self.CALL_FUNCTION_KW     if varkwarg else
          self.CALL_FUNCTION
        )(
            len(args) + 256 * len(kwargs) + preloaded,
           -len(args) -   2 * len(kwargs) - preloaded - len(vararg + varkwarg)
        )

    # store_top :: StructMixIn -> ()
    #
    # Append a sequence of opcodes that stores the top value in a variable.
    #
    def store_top(self, var):

        type, var, args = var._store_top_scanner_cache

        if type == const.AT.UNPACK:

            ln, star = args
            op  = 'UNPACK_SEQUENCE' if star < 0 else 'UNPACK_EX'
            arg = ln                if star < 0 else star + 256 * (ln - star - 1)
            self.appendcode(op, arg, ln - 1)

            for item in var:

                self.store_top(item)

        elif type == const.AT.ATTR:

            self.push(args)
            self.STORE_ATTR(var)

        elif type == const.AT.ITEM:

            self.push(args)
            self.push(var)
            self.STORE_SUBSCR()

        else:

            self.STORE_FAST (self.varnames[item]) if item in self.varnames else \
            self.STORE_DEREF(self.freevars[item]) if item in self.freevars else \
            self.STORE_DEREF(self.freevars[item]) if item in self.cellvars else \
            self.STORE_NAME (self.names[item])    if item in self.names    else \
            self.error(AssertionError, 'variable scanner error')

    ### ESSENTIAL BUILT-INS
    #
    # type Builtin = (StructMixIn, StructMixIn+) -> ()
    #

    # chain :: Builtin
    #
    # Push & pop each expression, except for the last one.
    #
    def chain(self, _, *args):

        for not_last, expr in enumerate(args, 1 - len(args)):

            self.push(expr)
            self.POP_TOP() if not_last else None

    chain.scanvars = joinvars1

    # getattr :: Builtin
    #
    # Retrieve an attribute of some object.
    #
    def getattr(self, f, obj, *args):

        len(args) > 0 or parse.syntax.error('which attribute?', obj)
        len(args) < 2 or parse.syntax.error('one at a time',    args[1])
        isinstance(args[0], parse.tree.Link) or parse.syntax.error('not an attribute', args[0])

        self.push(obj)
        self.LOAD_ATTR(self.appendname(args[0]))

    getattr.scanvars = lambda code, cell, nolocals: scanvars(code[1], cell, nolocals)

# _ :: dict StructMixIn ((Code, *StructMixIn) -> ())
#
# Map macro name to its implementation.
#
# INFIX_LEFT[R]  : `arg R`     where `R` is infix
# INFIX_RIGHT[R] : `R arg ...` where `R` is infix
# PREFIX[R]      : `R arg ...` or `arg1 R arg2`
#
# ---
#
# Each macro implementation may reimplement the variable scanning algorithm
# by supplying a `scanvars`-like function as its `scanvars` attribute.
#
# For example, `switch (a = b) (c = d)` may want to return joined results of
# `scanvars(a)`, `scanvars(b)`, `scanvars(c)`, and `scanvars(d)`, while the
# default behavior would be to do `scanvars(a = b)` and `scanvars(c = d)`.
#
INFIX_LEFT  = collections.defaultdict(lambda: Code.infixbindl)
INFIX_RIGHT = collections.defaultdict(lambda: Code.infixbindr)
PREFIX      = collections.defaultdict(lambda: Code.nativecall)

PREFIX['']   = Code.call
PREFIX['\n'] = Code.chain
PREFIX['.']  = Code.getattr
