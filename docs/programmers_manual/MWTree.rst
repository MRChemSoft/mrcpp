------------------
MWTree
------------------

This is an introduction to the mwtree class. We write a small overarching summary of the class where we define the
algorithm/equation/structure reasoning for having this class or where it fits with the rest of the code.

.. graphviz::

   digraph {
       "MWTree" -> "FunctionTree"
       "MWTree" -> "OperatorTree"
   }

.. doxygenclass:: mrcpp::MWTree
   :members:
   :protected-members:
   :private-members:
