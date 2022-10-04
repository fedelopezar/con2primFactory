# con2primFactory

## Branch: AsterX

In this branch, virtual functions of parent class 
have been abandoned.

Instead, we define `con2prim` algorithms as templates
for the different instances of the parent class.
This is crucial for running on GPUs.

Some code design is compromised, as now
child classes must be defined in `con2primFactory.hxx`,
and explicit instances of the template algorithms must
be defined in them same files.

