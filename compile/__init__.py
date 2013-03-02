from .  import core
from .. import parse


# ast :: StructMixIn -> code
# fd  :: TextIOBase  -> code
# it  :: str         -> code
#
# Compile something.
#
def ast(code, name='<module>'):

    target = core.Code('<unknown>', 0, 0, core.Code.NOFREE, *core.scanvars(code, set(), nolocals=True))
    target.push(code)
    target.RETURN_VALUE()
    return target.immutable


def fd(fd, name='<stream>'):

    return ast(parse.fd(fd, name))


def it(text, name='<string>'):

    return ast(parse.it(text, name))
