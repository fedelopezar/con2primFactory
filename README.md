# con2primFactory

Suite of conservative to primitive algorithms, agnostic of the plasma model.

See `idealFluid/` for an instantiation and test example.

## Current Branch: AsterX

In this branch, guided by GPU compatibility, virtual functions of
the parent class have been abandoned.
Instead, we define `con2prim` algorithms as 
templates functions for the different instances of
the parent class.

Also, the code is simplified and portable to the Einstein Toolkit.
For instance, `std::vectors` are moved to `C` 
arrays.

Some code design is compromised, as now
child classes must be defined in `con2primFactory.hxx`,
and explicit instances of the template algorithms must
be defined in them same files.

