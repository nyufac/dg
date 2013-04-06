import os
import imp
import marshal

from .  import core
from .. import parse


# ast :: StructMixIn -> code
# fd  :: TextIOBase  -> code
# it  :: str         -> code
#
# Compile something.
#
def ast(code, name='<module>'):

    queue = []
    a, b, c, d = core.scanvars(code, set(), True, queue)

    for f in queue:

        f(a, b, c, d)

    target = core.Code(name, 0, 0, 0, a, b, c, d)
    target.push(code)
    target.RETURN_VALUE()
    return target.immutable


def fd(fd, name='<stream>'):

    return ast(parse.fd(fd, name))


def it(text, name='<string>'):

    return ast(parse.it(text, name))


for f in [
       'shortcuts.dg'
#  , 'conditionals.dg',      'unary.dg'
#  ,       'binary.dg', 'comparison.dg'
#  ,      'inherit.dg',     'switch.dg', 'where.dg'
#  ,        'loops.dg',     'unsafe.dg',  'with.dg', 'yield.dg'
#  ,      'imphook.dg', 'functional.dg'
]:
    f = os.path.join(__path__[0], f)
    q = imp.cache_from_source(f)

    try:

        c = os.stat(q).st_mtime > os.stat(f).st_mtime and marshal.load(open(q, 'rb'))

    except Exception:

        c = None

    if not c:

        c = fd(open(f))
        os.makedirs(os.path.dirname(q), exist_ok=True)
        marshal.dump(c, open(q, 'wb'))

    eval(c, {'__package__': __package__})
