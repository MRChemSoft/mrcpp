

-----------
MWOperators
-----------

The MW operators discussed in this chapter is available to the application
program by including:

.. code-block:: cpp

    #include "MRCPP/MWOperators"


ConvolutionOperator
-------------------

.. NOTE::

    The convolution operators have separate precision parameters for their
    construction and application. The ``build_prec`` argument to the operator
    constructors will affect e.g. the number of terms in the separated
    representations of the Poisson/Helmholtz approximations, as well as the
    operator bandwidth. The ``apply_prec`` argument to the apply function relates
    only to the adaptive construction of the output function, based on a wavelet
    norm error estimate.

.. doxygenclass:: mrcpp::IdentityConvolution
   :members:

.. doxygenclass:: mrcpp::DerivativeConvolution
   :members:

.. doxygenclass:: mrcpp::PoissonOperator
   :members:

.. doxygenclass:: mrcpp::HelmholtzOperator
   :members:

.. doxygenfunction:: mrcpp::apply(double, FunctionTree<D>&, ConvolutionOperator<D>&, FunctionTree<D>&, int, bool)
.. doxygenfunction:: mrcpp::apply(double, FunctionTree<D>&, ConvolutionOperator<D>&, FunctionTree<D>&, FunctionTreeVector<D>&, int, bool)

DerivativeOperators
-------------------

.. NOTE::

    The derivative operators have clearly defined requirements on the output
    grid structure, based on the grid of the input function. This means that
    there is no real grid adaptivity, and thus no precision parameter is needed
    for the application of such an operator.

.. doxygenclass:: mrcpp::ABGVOperator
   :members:

.. doxygenclass:: mrcpp::PHOperator
   :members:

.. doxygenclass:: mrcpp::BSOperator
   :members:

.. doxygenfunction:: mrcpp::apply(FunctionTree<D>&, DerivativeOperator<D>&, FunctionTree<D>&, int)
.. doxygenfunction:: mrcpp::divergence(FunctionTree<D>&, DerivativeOperator<D>&, FunctionTreeVector<D>&)
.. doxygenfunction:: mrcpp::gradient(DerivativeOperator<D>&, FunctionTree<D>&)


Examples
--------

PoissonOperator
+++++++++++++++

The electrostatic potential :math:`g` arising from a charge distribution
:math:`f` are related through the Poisson equation

.. math:: -\nabla^2 g(r) = f(r)

This equation can be solved with respect to the potential by inverting the
differential operator into the Green's function integral convolution operator

.. math:: g(r) =  \int \frac{1}{4\pi\|r-r'\|} f(r') dr'

This operator is available in the MW representation, and can be solved with
arbitrary (finite) precision in linear complexity with respect to system size.
Given an arbitrary charge dirtribution ``f_tree`` in the MW representation, the
potential is computed in the following way:

.. code-block:: cpp

    double apply_prec;                              // Precision for operator application
    double build_prec;                              // Precision for operator construction

    mrcpp::PoissonOperator P(MRA, build_prec);      // MW representation of Poisson operator
    mrcpp::FunctionTree<3> f_tree(MRA);             // Input function
    mrcpp::FunctionTree<3> g_tree(MRA);             // Output function

    mrcpp::apply(apply_prec, g_tree, P, f_tree);    // Apply operator adaptively

The Coulomb self-interaction energy can now be computed as the dot product:

.. code-block:: cpp

    double E = mrcpp::dot(g_tree, f_tree);

HelmholtzOperator
+++++++++++++++++

The Helmholtz operator is a generalization of the Poisson operator and is given
as the integral convolution

.. math:: g(r) =  \int \frac{e^{-\mu\|r-r'\|}}{4\pi\|r-r'\|} f(r') dr'

The operator is the inverse of the shifted Laplacian

.. math:: \big[-\nabla^2 + \mu^2 \big] g(r) = f(r)

and appears e.g. when solving the SCF equations. The construction and
application is similar to the Poisson operator, with an extra argument for the
:math:`\mu` parameter

.. code-block:: cpp

    double apply_prec;                              // Precision for operator application
    double build_prec;                              // Precision for operator construction
    double mu;                                      // Must be a positive real number

    mrcpp::HelmholtzOperator H(MRA, mu, build_prec);// MW representation of Helmholtz operator
    mrcpp::FunctionTree<3> f_tree(MRA);             // Input function
    mrcpp::FunctionTree<3> g_tree(MRA);             // Output function

    mrcpp::apply(apply_prec, g_tree, H, f_tree);    // Apply operator adaptively


ABGVOperator
++++++++++++

The ABGV (Alpert, Beylkin, Gines, Vozovoi) derivative operator is initialized
with two parameters :math:`a` and :math:`b` accounting for the boundary
conditions between adjacent nodes, see `Alpert et al.
<http://www.sciencedirect.com/science/article/pii/S0021999102971603>`_

.. code-block:: cpp

    double a = 0.0, b = 0.0;                        // Boundary conditions for operator
    mrcpp::ABGVOperator<3> D(MRA, a, b);            // MW derivative operator
    mrcpp::FunctionTree<3> f(MRA);                  // Input function
    mrcpp::FunctionTree<3> f_x(MRA);                // Output function
    mrcpp::FunctionTree<3> f_y(MRA);                // Output function
    mrcpp::FunctionTree<3> f_z(MRA);                // Output function

    mrcpp::apply(f_x, D, f, 0);                     // Operator application in x direction
    mrcpp::apply(f_y, D, f, 1);                     // Operator application in y direction
    mrcpp::apply(f_z, D, f, 2);                     // Operator application in z direction

The tree structure of the output function will depend on the choice of
parameters :math:`a` and :math:`b`: if both are zero, the output grid will be
identical to the input grid; otherwise the grid will be widened by one node (on
each side) in the direction of application.


PHOperator
++++++++++

The PH derivative operator is based on the noise reducing derivative of `Pavel Holoborodko
<http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/smooth-low-noise-differentiators/>`_.
This operator is also available as a direct second derivative.


.. code-block:: cpp

    mrcpp::PHOperator<3> D1(MRA, 1);                // MW 1st derivative operator
    mrcpp::PHOperator<3> D2(MRA, 2);                // MW 2nd derivative operator
    mrcpp::FunctionTree<3> f(MRA);                  // Input function
    mrcpp::FunctionTree<3> f_x(MRA);                // Output function
    mrcpp::FunctionTree<3> f_xx(MRA);               // Output function

    mrcpp::apply(f_x, D1, f, 0);                    // Operator application in x direction
    mrcpp::apply(f_xx, D2, f, 0);                   // Operator application in x direction


Special thanks to Prof. Robert J. Harrison (Stony Brook University) for sharing the
operator coefficients.

BSOperator
++++++++++

The BS derivative operator is based on a pre-projection onto B-splines in order
to remove the discontinuities in the MW basis, see `Anderson et al.
<https://www.sciencedirect.com/science/article/pii/S2590055219300496>`_
This operator is also available as a direct second and third derivative.


.. code-block:: cpp

    mrcpp::BSOperator<3> D1(MRA, 1);                // MW 1st derivative operator
    mrcpp::BSOperator<3> D2(MRA, 2);                // MW 2nd derivative operator
    mrcpp::BSOperator<3> D3(MRA, 3);                // MW 3nd derivative operator
    mrcpp::FunctionTree<3> f(MRA);                  // Input function
    mrcpp::FunctionTree<3> f_x(MRA);                // Output function
    mrcpp::FunctionTree<3> f_yy(MRA);               // Output function
    mrcpp::FunctionTree<3> f_zzz(MRA);              // Output function

    mrcpp::apply(f_x, D1, f, 0);                    // Operator application in x direction
    mrcpp::apply(f_yy, D2, f, 1);                   // Operator application in x direction
    mrcpp::apply(f_zzz, D3, f, 2);                  // Operator application in x direction

