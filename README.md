# con2primFactory

Suite of conservative to primitive algorithms, agnostic of the plasma model.

See AsterX_Flat_idealFluid_Magnetized/ for an instantiation and test example.

## Current Branch: AsterX

In this branch, guided by GPU compatibility, virtual functions of
the parent class have been abandoned.

Instead, we define `con2prim` algorithms as 
templates functions for the different instances of
the parent class.

Some code design is compromised, as now
child classes must be defined in `con2primFactory.hxx`,
and explicit instances of the template algorithms must
be defined in them same files.

