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
    a, b, c, d = core.scanvars(code, set(), nolocals=True, queue=queue)

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
