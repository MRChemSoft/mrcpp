------------------
MultiResolutionAnalysis
------------------

The MultiResolutionAnalysis (MRA) class contains the methods to project objects onto the spatial grid.
That is, to combine different functions and operators in mathematical operations, they need to be compatible; 
they must be defined on the same computational domain and constructed using the same polynomial basis 
(order and type). This information constitutes an MRA, which needs to be defined and passed as argument to 
all function and operator constructors, and only functions and operators with compatible MRAs can be 
combined in subsequent calculations.

.. doxygenclass:: mrcpp::MultiResolutionAnalysis
   :members:
   :protected-members:
   :private-members:

