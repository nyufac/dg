import io
import dis
import types
import collections

from .. import parse

__all__ = ['INFIX_LEFT', 'INFIX_RIGHT', 'PREFIX', 'scanvars', 'compile', 'Code']


# scanvars :: (StructMixIn, set Link) -> (set Link, set Link, set Link, set Link)
#
# Scan an AST for variable names. The returned sets contain:
#   1. global variables
#   2. local variables
#   3. cell variables
#   4. free variables
#
def scanvars(code, cell, rightbind=False):

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
        else PREFIX[f]

        if hasattr(g, 'scanvars'):

            return g.scanvars(code, cell)

        globals = set()
        locals  = set()
        ncell   = set()
        free    = set()

        for item in code:

            a, b, c, d = scanvars(item, cell)
            globals |= a
            globals -= b | c | d
            locals  |= b
            locals  -= c | d
            ncell   |= c
            ncell   -= d
            free    |= d

        return globals, locals, ncell, free

    raise TypeError('not an AST output structure: {!r}'.format(code))


# compile :: (Either StructMixIn TextIOBase str, Code, Maybe str) -> Code
#
# Compile a parse tree into a mutable code object.
#
def compile(code, target, name='<unknown>'):

    if isinstance(code, parse.tree.StructMixIn):

        target.push(code)
        target.RETURN_VALUE()
        return target

    elif isinstance(code, io.TextIOBase):

        return compile(parse.fd(code, name), target, name)

    elif isinstance(code, str):

        return compile(parse.it(code, name), target, name)

    else:

        raise TypeError('not compilable: {!r}'.format(code))


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

    def __init__(self, name, argc, kwargc, flags, globals, locals, cell, free):

        super().__init__()
        self.name  = name
        self.flags = flags

        self.argc   = argc
        self.kwargc = kwargc

        self.filename = '<generated>'
        self.lineno   = 0
        self.depth    = 10
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

        #raise NotImplementedError

        #stackdelta = self.OPCODE_STACK_DELTA[name] if stackdelta is None else stackdelta
        self.bytecode.append((dis.opmap[name], argument))

    @property
    # bakedcode :: iter bytes
    #
    # Compiled bytecode string.
    #
    def bakedcode(self):

        for code, arg in self.bytecode:

            yield code

            if code >= dis.HAVE_ARGUMENT:

                yield arg %  256
                yield arg // 256

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

        oldstacksize = self.depth

        if isinstance(item, parse.tree.Link):

            self.appendcode('LOAD_FAST',  self.varnames[item]) if item in self.varnames else \
            self.appendcode('LOAD_DEREF', self.freevars[item]) if item in self.freevars else \
            self.appendcode('LOAD_DEREF', self.freevars[item]) if item in self.cellvars else \
            self.appendcode('LOAD_NAME',  self.names[item])    if item in self.names    else \
            self.error(AssertionError, 'variable scanner error')

        elif isinstance(item, parse.tree.Constant):

            val = item.value
            self.consts[type(val), val] = a = self.consts.get((type(val), val), len(self.consts))
            self.appendcode('LOAD_CONST', a)

        else:

            self.call(*item)

        #assert oldstacksize + 1 == self.depth, 'stack depth calculation error'

    # call :: *StructMixIn -> ()
    #
    # Call a (possibly built-in) function.
    #
    def call(self, f, *args, rightbind=False):

        return self.call(*args, rightbind=True) if f.infix and f == '' \
          else INFIX_RIGHT[f](self, f, *args) if f.infix and not f.closed and rightbind \
          else INFIX_LEFT[f](self, f, *args)  if f.infix and not f.closed and len(args) == 1 \
          else PREFIX[f](self, f, *args)

    call.scanvars = lambda code, cell: scanvars(parse.tree.Expression(code[1:]), cell, rightbind=True)

    # infixbindl :: (StructMixIn, StructMixIn) -> ()
    #
    # Default implementation of a left infix bind.
    #
    def infixbindl(self, op, arg):

        raise NotImplementedError

    # infixbindr :: (StructMixIn, StructMixIn+) -> ()
    #
    # Default implementation of right infix bind.
    #
    def infixbindr(self, op, *args):

        raise NotImplementedError

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

    for opname in dis.opmap:

        locals()[opname] = lambda self, argument=None, stackdelta=None, name=opname: \
          self.appendcode(name, argument, stackdelta)

# _ :: dict StructMixIn ((Code, *StructMixIn) -> ())
#
# Map macro name to its implementation.
#
# INFIX_LEFT[R]  : `arg R`     where `R` is infix
# INFIX_RIGHT[R] : `R arg ...` where `R` is infix
# PREFIX[R]      : `R arg ...`
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

PREFIX[''] = Code.call
