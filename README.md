# con2primFactory

Suite of conservative to primitive algorithms, 
agnostic of the plasma model.

We are guided by the abstract factory pattern: 
The parent class `con2primFactory` serves as a
template for the model-dependent member functions
that the user must provide&mdash;these are virtual
functions to be overwritten.

See `AsterX_Flat_idealFluid_Magnetized/` for an 
instantiation and test example.

Issue: Virtual functions conflict with GPU performance. For an alternate design, see branch
`AsterX`.