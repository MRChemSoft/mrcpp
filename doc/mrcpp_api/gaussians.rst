
---------
Gaussians
---------

MRCPP provides some simple features for analytic Gaussian functions. These
are meant to be used as a starting point for MW computations, they are
*not* meant for heavy analytical computation, like GTO basis sets. The
Gaussian features are available by including:

.. code-block:: cpp

    #include "MRCPP/Gaussians"


Available functions
-------------------

evalf
  Evaluate function in a point.

calcOverlap
  Compute overlap with another Gaussian.

differentiate
  Compute analytical derivative.

getSquareNorm
  Returns the squared L2 norm.

normalize
  Rescale the function by its norm.

mult
  Multiply two Gaussian functions.

calcCoulombEnergy
  Compute analytical electrostatic energy between two Gaussian charge
  distributions.


GaussFunc
---------

A ``GaussFunc`` is a simple D-dimensional Gaussian function with a Cartesian
monomial in front, e.g. in 3D:


.. math:: f(r) = \alpha (x-x_0)^a (y-y_0)^b (z-z_0)^c e^{-\beta \|r-r_0\|^2}

.. code-block:: cpp

    double alpha, beta;
    std::array<int, 3> pow = {a, b, c};
    std::array<double, 3> pos = {x_0, y_0, z_0};
    mrcpp::GaussFunc<3> gauss(beta, alpha, pos, pow);

    double E = gauss.calcCoulombEnergy(gauss);              // Analytical energy

This Gaussian function can be used to build an empty grid based on the position
and exponent. The grid will then be refined close to the center of the Gaussian,
with deeper refinement for higher exponents (steeper function):

.. code-block:: cpp

    mrcpp::FunctionTree<3> g_tree(MRA);
    mrcpp::build_grid(g_tree, gauss);                       // Build empty grid
    mrcpp::project(prec, g_tree, gauss);                    // Project Gaussian

GaussPoly
---------

``GaussPoly`` is a generalization of the ``GaussFunc``, where there is an
arbitrary polynomial in front of the exponential

.. math:: f(r) = \alpha P(r-r_0) e^{-\beta \|r-r_0\|^2}

for instance making the GaussPoly:

.. math:: f(r) = \alpha (a_x + b_x x + c_x x^2) (a_y + b_y y + c_y y^2) (a_z + b_z z + c_z z^2)e^{-\beta \|r-r_0\|^2}



.. code-block:: cpp

    auto gauss_poly = GaussPoly<D>(beta, alpha, pos, pow);

    // Create polynomial in x, y and z direction
    auto pol_x = Polynomial(2); // 2 is the order of the polynomial
    pol_x.getCoefs() << a_x, b_x, c_x;
    auto pol_y = Polynomial(2);
    pol_y.getCoefs() << a_y, b_y, c_y;
    auto pol_z = Polynomial(2);
    pol_z.getCoefs() << a_z, b_z, c_z;

    // Add polynomials to gauss_poly
    guass_poly.setPoly(0, pol_x);
    guass_poly.setPoly(1, pol_y);
    guass_poly.setPoly(2, pol_z);

GaussExp
--------

A ``GaussExp`` is a collection of Gaussians in the form

.. math:: G(r) = \sum_i c_i g_i(r)

where :math:`g_i` can be either ``GaussFunc`` or ``GaussPoly``

.. math:: g_i(r) =  \alpha_i P_i(r-r_i)e^{-\beta_i\|r-r_i\|^2}

Individual Gaussian functions can be appended to the ``GaussExp`` and treated as
a single function:

.. code-block:: cpp

    mrcpp::GaussExp<3> g_exp;                               // Empty Gaussian expansion
    for (int i = 0; i < N; i++) {
        double alpha_i, beta_i;                             // Individual parameters
        std::array<int, 3> pow_i;                           // Individual parameters
        std::array<double, 3> pos_i;                        // Individual parameters
        mrcpp::GaussFunc<3> gauss_i(beta_i, alpha_i, pos_i, pow_i);
        g_exp.append(gauss_i);                              // Append Gaussian to expansion
    }
    mrcpp::project(prec, tree, g_exp);                      // Project full expansion
