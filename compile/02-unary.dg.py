..compile = import


compile.r.builtins !! 'return' = (self, a) -> self.opcode: 'RETURN_VALUE' None a delta: 1
compile.r.builtins !! 'yield'  = (self, a) -> self.opcode: 'YIELD_VALUE'       a delta: 1
compile.r.builtins !! 'not'    = (self, a) -> self.opcode: 'UNARY_NOT'         a delta: 1
compile.r.builtins !! '~'      = (self, a) -> self.opcode: 'UNARY_INVERT'      a delta: 1