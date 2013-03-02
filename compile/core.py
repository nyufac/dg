from io import TextIOBase

from .. import parse


# scanvars :: (StructMixIn, set Link) -> (set Link, set Link, set Link, set Link, set (type a, a))
#
# Scan an AST for variable names. The returned sets contain:
#   1. global variables
#   2. local variables
#   3. cell variables
#   4. free variables
#   5. constant (type, value) pairs
#
def scanvars(code, cell):

    if isinstance(code, parse.tree.Link):

        return (
          (set(), set(), set(), {code}, set()) if code in cell else
          ({code}, set(), set(), set(), set())
        )

    if isinstance(code, parse.tree.Constant):

        return set(), set(), set(), set(), {(type(code.value), code.value)}

    globals  = set()
    locals   = set()
    closures = set()

    raise NotImplementedError


# compile :: (Either StructMixIn TextIOBase str, Code, Maybe str) -> Code
#
# Compile a parse tree into a mutable code object.
#
def compile(code, target, name='<unknown>'):

    if isinstance(code, parse.tree.StructMixIn):

        target.RETURN_VALUE(code)
        return target

    elif isinstance(code, TextIOBase):

        return compile(parse.fd(code, name), target, name)

    elif isinstance(code, str):

        return compile(parse.it(code, name), target, name)

    else:

        raise TypeError('not compilable: {!r}'.format(code))


# data Code = Code (str, int, int, int, set Link, set Link, set Link, set Link, set (type a, a))
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

    def __init__(self, name, argc, kwargc, flags, globals, locals, cell, free, consts):

        super().__init__()
        self.name  = name
        self.flags = flags

        self.argc   = argc
        self.kwargc = kwargc

        self.lineno   = 0
        self.depth    = 0
        self.depthmax = 0
        self.bytecode = []
        self.lnomap   = {}

        # Just in case.
        # Run with `-O` (or was it `-OO`?) to ignore these checks.
        assert len(locals) >= argc + kwargc + bool(flags & self.VARARGS) + bool(flags & VARKWARGS)
        assert bool(flags & self.NOFREE) == (not cell)
        assert bool(flags & self.NESTED) == bool(free)

        assert not (free & cell)
        assert not (free & locals)
        assert not (cell & locals)
        assert not (free & globals)
        assert not (cell & globals)
        assert not (locals & globals)

        self.consts   = dict(zip(consts,  range(len(consts))))
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

        raise NotImplementedError

        stackdelta = self.OPCODE_STACK_DELTA[name] if stackdelta is None else stackdelta
        ...

    # push :: StructMixIn -> ()
    #
    # Append a sequence of opcodes that pushes the result of an expression
    # onto the stack.
    #
    def push(self, item):

        if isinstance(item, parse.tree.Link):

            self.appendcode('LOAD_FAST',  self.varnames[item]) if item in self.varnames else \
            self.appendcode('LOAD_DEREF', self.freevars[item]) if item in self.freevars else \
            self.appendcode('LOAD_DEREF', self.freevars[item]) if item in self.cellvars else \
            self.appendcode('LOAD_NAME',  self.names[item])    if item in self.names    else \
            self.error(AssertionError, 'variable scanner error')

        elif isinstance(item, parse.tree.Constant):

            self.appendcode('LOAD_CONST', self.consts[type(item.value), item.value])

        else:

            raise NotImplementedError
