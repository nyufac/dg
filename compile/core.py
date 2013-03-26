import dis
import types
import posixpath
import itertools

from .. import parse
from .. import const

__all__ = ['INFIX_LEFT', 'INFIX_RIGHT', 'PREFIX', 'scanvars', 'compile', 'Code']


# throw :: BaseException -> _|_
#
# An expression equivalent to `raise`.
#
def throw(e): raise e

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
def scanvars(code, cell, nolocals, queue):

    assert isinstance(code, parse.tree.StructMixIn)
    return (
        (set(),  set(), set(), {code}) if isinstance(code, parse.tree.Link) and code in cell else
        ({code}, set(), set(), set())  if isinstance(code, parse.tree.Link)                  else
        (set(),  set(), set(), set())  if isinstance(code, parse.tree.Constant)              else
        Code.call.scanvars(code, cell, nolocals, queue, rightbind=False)
    )


# joinsets :: [(set, set, set, set)] -> (set, set, set, set)
#
# Combine the results of scanning multiple expressions.
#
def joinsets(sets):

    globals, locals, cell, free = map(set.union, *sets)
    cell    -= free
    locals  -= free | cell
    globals -= free | cell | locals
    return globals, locals, cell, free


# joinvars :: ([StructMixIn], set Link) -> (set Link, set Link, set Link, set Link)
#
# Like `scanvars`, but for multiple items in the same expression.
#
def joinvars(code, cell, nolocals, queue, scanf=scanvars):

    return joinsets(scanf(item, cell, nolocals, queue) for item in code)


# joinvars1 :: ([StructMixIn], set Link) -> (set Link, set Link, set Link, set Link)
#
# Like `joinvars`, but ignores the first item, which is generally the function name.
#
def joinvars1(code, cell, nolocals, queue):

    return joinvars(code[1:], cell, nolocals, queue)


# joinvarsX :: (*StructMixIn, Maybe bool) -> typeof joinvars
#
# A wrapper for `joinvars` and `joinvars1` that adds a few more items.
# Useful when a compile-time function uses variables.
#
def joinvarsX(*add, f=joinvars):

    add = list(add)
    return lambda code, cell, nolocals, queue: f(add + list(code), cell, nolocals, queue)


# binary_scanner :: (str, str) -> typeof scanvars -> typeof scanvars
#
# Restrict a built-in function to two arguments by raising `SyntaxError`s.
#
def binary_scanner(one='this needs 2 arguments', three='too many arguments for {!r}'):

    def wrapper(method):

        def wrapped(code, cell, nolocals, queue):

            f, lhs, *rhs = code
            len(rhs) > 0 or parse.syntax.error(one  .format(f), f)
            len(rhs) < 2 or parse.syntax.error(three.format(f), rhs[1])
            return method(code, cell, nolocals, queue)

        return wrapped

    return wrapper


# scanner_of :: ((StructMixIn, StructMixIn+) -> ()) -> typeof scanvars -> NoneType
#
# Assign a scanner to a built-in function.
#
def scanner_of(f):

    return lambda x: setattr(f, 'scanvars', x)

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

    # `CodeType` flags. See CPython source code for information.
    OPTIMIZED = 1
    NEWLOCALS = 2
    VARARGS   = 4
    VARKWARGS = 8
    NESTED    = 16
    GENERATOR = 32
    NOFREE    = 64

    # How many items an opcode will push onto the stack.
    # Note that some opcodes are missing from this list, as they change
    # behavior depending on the argument (e.g. `END_FINALLY` may pop
    # either one or three objects.)
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
     'INPLACE_ADD':          -1,  # NOTE these opcodes, despite having `INPLACE`
     'INPLACE_AND':          -1,  #   in their names, still require `STORE_TOP`
     'INPLACE_FLOOR_DIVIDE': -1,  #   to be added afterwards. The only
     'INPLACE_LSHIFT':       -1,  #   difference from `BINARY` opcodes
     'INPLACE_MODULO':       -1,  #   is the method used (e.g. `__iadd__`
     'INPLACE_MULTIPLY':     -1,  #   instead of `__add__`.)
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
        self.flags = flags | self.NOFREE * (not cell) | self.NESTED * bool(free)
        self.noldname = flags & self.OPTIMIZED

        self.argc   = argc
        self.kwargc = kwargc

        self.filename = '<generated>'
        self.lineno   = 0
        self.depth    = 0
        self.depthmax = 0
        self.bytecode = []
        self.lnomap   = {}

        # Run with `-O` to ignore these checks, though that's not recommended.
        assert len(locals) >= argc + kwargc + bool(flags & self.VARARGS) + bool(flags & self.VARKWARGS)
        assert bool(self.flags & self.NOFREE) == (not cell)
        assert bool(self.flags & self.NESTED) == bool(free)

        assert not (free & cell)
        assert not (free & locals)
        assert not (cell & locals)
        assert not (free & globals)
        assert not (cell & globals)
        assert not (locals & globals)

        self.consts   = {}
        self.names    = {k: i for i, k in enumerate(globals)}
        self.varnames = {k: i for i, k in enumerate(locals)}
        self.cellvars = {k: i for i, k in enumerate(cell)}
        self.freevars = {k: i + len(cell) for i, k in enumerate(free)}

    # appendcode :: (str, Maybe (Either int DelayedComputation), Maybe int) -> ()
    #
    # Append a single opcode given its name (sometimes with an argument.)
    #
    def appendcode(self, name, argument=None, stackdelta=None):

        assert stackdelta is not None or name in self.OPCODE_STACK_DELTA
        self.depth += stackdelta if stackdelta is not None else self.OPCODE_STACK_DELTA[name]
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

    ## HIGH-LEVEL OPCODE GENERATION

    # push :: StructMixIn -> ()
    #
    # Append a sequence of opcodes that pushes the result of an expression
    # onto the stack.
    #
    def push(self, item):

        assert isinstance(item, parse.tree.StructMixIn)
        mdepth = self.depth

        if isinstance(item, parse.tree.Link):

            self.LOAD_DEREF (self.freevars[item]) if item in self.freevars else \
            self.LOAD_DEREF (self.cellvars[item]) if item in self.cellvars else \
            self.LOAD_FAST  (self.varnames[item]) if item in self.varnames else \
            self.LOAD_GLOBAL(self.names[item])    if self.noldname         else \
            self.LOAD_NAME  (self.names[item])    if item in self.names    else \
            throw(AssertionError('variable scanner error'))

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

        return self.call(*args, rightbind=True) if f.infix and f == '' else f._f(self, f, *args)

    @scanner_of(call)
    def _(code, cell, nolocals, queue, rightbind=True):

        f, *args = code
        f._f = INFIX_RIGHT.get(f, Code.infixbindr) if f.infix and not f.closed and rightbind \
          else INFIX_LEFT.get(f, Code.infixbindl)  if f.infix and not f.closed and len(args) == 1 \
          else PREFIX.get(f, Code.nativecall) if isinstance(f, parse.tree.Link) else Code.nativecall

        code = args if f.infix and f == '' else code
        return getattr(f._f, 'scanvars', joinvars)(code, cell, nolocals, queue)

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
            self.STORE_ATTR(self.appendname(var))

        elif type == const.AT.ITEM:

            self.push(args)
            self.push(var)
            self.STORE_SUBSCR()

        else:

            self.STORE_DEREF(self.freevars[var]) if var in self.freevars else \
            self.STORE_DEREF(self.cellvars[var]) if var in self.cellvars else \
            self.STORE_FAST (self.varnames[var]) if var in self.varnames else \
            self.STORE_NAME (self.names[var])    if var in self.names and not self.noldname else \
            throw(AssertionError('variable scanner error'))

    @scanner_of(store_top)
    # NOTE since `store_top` is not a real built-in, you
    #   should call this manually::
    #
    #     Code.store_top.scanvars(varname, cell, nolocals)
    #
    def _(code, cell, nolocals, queue):

        type, var, args = code._store_top_scanner_cache = parse.syntax.assignment_target(code)

        return (
            joinvars([args, var], cell, nolocals, queue)                         if type == const.AT.ITEM   else
            joinvars(var,  cell, nolocals, queue, scanf=Code.store_top.scanvars) if type == const.AT.UNPACK else
            scanvars(args, cell, nolocals, queue)                                if type == const.AT.ATTR   else
            (set(), set(), set(), {var}) if var in cell else
            ({var}, set(), set(), set()) if nolocals else
            (set(), {var}, set(), set())
        )


    # make_function :: (Code, [StructMixIn], dict Link StructMixIn) -> ()
    #
    # Push a function object onto the stack.
    #
    def make_function(self, target, defs, kwdefs):

        for k, v in kwdefs.items():

            self.push(parse.tree.Constant(str(k)))
            self.push(v)

        for a in defs:

            self.push(a)

        if target.freevars:

            for v in target.freevars:

                self.LOAD_CLOSURE(self.freevars[v] if v in self.freevars else self.cellvars[v])

            self.BUILD_TUPLE(len(target.freevars), 1 - len(target.freevars))

        self.push(parse.tree.Constant(target.immutable))
        self.push(parse.tree.Constant('<FIXME:qualname>'))
        (self.MAKE_CLOSURE if target.freevars else self.MAKE_FUNCTION)(
                len(defs) + 256 * len(kwdefs),
            1 - len(defs) -   2 * len(kwdefs) - bool(target.freevars) - 2
        )

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
    def getattr(self, _, obj, attr):

        self.push(obj)
        self.LOAD_ATTR(self.appendname(attr))

    @scanner_of(getattr)
    @binary_scanner('which attribute?', 'one at a time')
    def _(code, cell, nolocals, queue):

        _, lhs, rhs = code
        isinstance(rhs, parse.tree.Link) or parse.syntax.error('not an attribute', rhs)
        return scanvars(lhs, cell, nolocals, queue)

    # store :: Builtin
    #
    # Store the result of an expression in `var`.
    #
    def store(self, _, var, expr):

        #e = self.assigned_to
        #self.assigned_to = var
        self.push(expr)
        #self.assigned_to = e
        self.DUP_TOP()
        self.store_top(var)

    @scanner_of(store)
    @binary_scanner('store what?', 'one at a time')
    def _(code, cell, nolocals, queue):

        _, var, expr = code
        return joinsets([
            scanvars(expr, cell, nolocals, queue),
            Code.store_top.scanvars(var, cell, nolocals, queue)
        ])

    # function :: Builtin
    #
    # Create a function given an argument specification and a body.
    #
    def function(self, f, argspec, body):

        code = Code('<FIXME:name>', argspec._argc, argspec._kwargc, Code.OPTIMIZED | Code.NEWLOCALS, *body._vars)
        # Argument names should be in the correct order.
        code.varnames = {k: i + len(argspec._argnames) for i, k in enumerate(code.varnames.keys() - set(argspec._argnames))}
        code.varnames.update({k: i for i, k in enumerate(argspec._argnames)})

        for name, pattern in argspec._targets.items():

            code.LOAD_FAST(code.varnames[name])
            code.store_top(pattern)

        code.push(body)
        code.RETURN_VALUE()

        self.make_function(code, argspec._defs, argspec._kwdefs)

    @scanner_of(function)
    @binary_scanner('this function should do what?', 'one body per function')
    def _(code, cell, nolocals, queue):

        _, argspec, body = code
        argspec._argspec_cache = parse.syntax.argspec(argspec, definition=True)
        args, kwargs, defs, kwdefs, varargs, varkwargs = argspec._argspec_cache

        argspec._argc     = len(args)
        argspec._kwargc   = len(kwargs)
        argspec._flags    = Code.VARARGS * len(varargs) | Code.VARKWARGS | len(varkwargs)
        argspec._defs     = defs
        argspec._kwdefs   = kwdefs
        argspec._argnames = []
        argspec._targets  = {}

        for i, arg in enumerate(itertools.chain(args, kwargs, varargs, varkwargs)):

            if not isinstance(arg, parse.tree.Link) or arg in argspec._argnames:

                arg, expr = 'pattern-{}'.format(i), arg
                argspec._targets[arg] = expr

            argspec._argnames.append(arg)

        def scanbody(_, locals, cell, free):

            new_queue = []
            body._vars = joinsets([
                (set(), set(argspec._argnames), set(), set()),
                scanvars(body, locals | free, False, new_queue)
            ])

            for f in new_queue:

                f(*body._vars)

            locals -= body._vars[-1]
            cell   |= body._vars[-1] - free

        queue.append(scanbody)
        return set(), set(), set(), set()

    # import_ :: Builtin
    #
    # Import a module given a POSIX-style path.
    # Normally, this function works as `from ... import ...` in Python;
    # however, supplying "qualified" keyword after a module name makes it
    # work like a simple `import`. (This only affects paths with multiple slashes.)
    #
    def import_(self, _, name, *qualified):

        path, parent = name._data
        self.push(parse.tree.Constant(parent))

        if qualified or len(path) == 1:

            self.push(parse.tree.Constant(None))
            self.IMPORT_NAME(self.appendname(parse.tree.Link('.'.join(path))))

        else:

            *mod, var = path
            self.push(parse.tree.Constant((str(var),)))
            self.IMPORT_NAME(self.appendname(parse.tree.Link('.'.join(mod))))
            self.IMPORT_FROM(self.appendname(var))
            self.ROT_TWO()
            self.POP_TOP()

        self.DUP_TOP()
        self.store_top(path[0 - len(qualified)])


    @scanner_of(import_)
    def _(code, cell, nolocals, queue):

        _, name, *qualified = code
        if not isinstance(name, parse.tree.Constant): parse.syntax.error('not a module name', name)
        if not isinstance(name.value, str):           parse.syntax.error('not a module name', name)
        if len(qualified) >  1:                                 parse.syntax.error('excess arguments', qualified[1])
        if len(qualified) == 1 and qualified[0] != 'qualified': parse.syntax.error('only `qualified` is allowed here', qualified[0])

        path   = [parse.tree.Link(x) for x in posixpath.normpath(name.value).split(posixpath.sep)]
        parent = 1

        while path and path[0] == '':               parent = path.pop(0) or 0
        while path and path[0] == posixpath.curdir:          path.pop(0)
        while path and path[0] == posixpath.pardir: parent = path.pop(0) and parent + 1
        if not path: parse.syntax.error('invalid module name', name)
        if qualified and parent: parse.syntax.error('cannot perform qualified relative imports', qualified[0])

        name._data = path, parent
        return Code.store_top.scanvars(path[0 - len(qualified)], cell, nolocals, queue)


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
INFIX_LEFT  = {}
INFIX_RIGHT = {}
PREFIX      = {
    '':       Code.call
  , '.':      Code.getattr
  , '=':      Code.store
  , '\n':     Code.chain
  , '->':     Code.function
  , 'import': Code.import_
}
