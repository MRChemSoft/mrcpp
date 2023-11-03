---------------------
Apply for complex valued functions
---------------------

Application of a complex convolution.


.. doxygenfunction:: mrcpp::apply
(
    double prec, ComplexObject< FunctionTree<D> > &out,
    ComplexObject< ConvolutionOperator<D> > &oper, ComplexObject< FunctionTree<D> > &inp,
    int maxIter, bool absPrec
)