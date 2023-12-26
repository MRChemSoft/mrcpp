------------------
CornerOperatorTree
------------------

This is an introduction to the CornerOperatorTree class. We write a small overarching summary of the class where we define the
algorithm/equation/structure reasoning for having this class or where it fits with the rest of the code.

This class inherits from OperatorTree and represents the non-standard form with
matrices :math:`A, B, C` having negligible matrix elements around their diagonals.
Note that only these three matrices can be considered narrow banded at the corners,
and so only for them we redefine the notion of a band width.
A band of one of them is the distance between the first non-negligible diagonal and the main diagonal.


.. doxygenclass:: mrcpp::CornerOperatorTree
   :members:
   :protected-members:
   :private-members:

