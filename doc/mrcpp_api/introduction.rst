------------
Introduction
------------

The main features of MRCPP are the numerical multiwavelet (MW) representations
of functions and operators. Two integral convolution operators are implemented
(the Poisson and Helmholtz operators), as well as the partial derivative and
arithmetic operators. In addition to the numerical representations there are a
limited number of analytic functions that are usually used as starting point
for the numerical computations. Also, MRCPP provides three convenience classes
(Timer, Printer and Plotter) that can be made available to the application
program.

The API consists of seven include files which will be discussed in more detail
below:

MRCPP/MWFunctions
  Provides features for representation and manipulation of real-valued
  scalar functions in a MW basis, including projection of analytic function,
  numerical integration and scalar products, as well as arithmetic operations
  and function mappings.

MRCPP/MWOperators
  Provides features for representation and application of MW operators.
  Currently there are three operators available: Poisson, Helmholtz and
  Cartesian derivative.

MRCPP/Gaussians
  Provides some simple features for analytical Gaussian functions, useful e.g.
  to generate initial guesses for MW computations.

MRCPP/Parallel
  Provides some simple MPI features for MRCPP, in particular the possibility to
  send complete MW function representations between MPI processes.

MRCPP/Printer
  Provides simple (parallel safe) printing options. All MRCPP internal printing
  is done with this class, and the printer must be initialized in order to get
  any printed output, otherwise MRCPP will run silently.

MRCPP/Plotter
  Provides options to generate data files for plotting of MW function
  representations. These include line plots, surface plots and cube plots, as
  well as grid visualization using geomview.

MRCPP/Timer
  Provides an accurate timer for the wall clock in parallel computations.


Analytic functions
------------------

The general way of defining an analytic function in MRCPP is to use lambdas
(or std::function), which provide lightweight functions that can be used on
the fly. However, some analytic functions, like Gaussians, are of special
importance and have been explicitly implemented with additional functionality
(see Gaussian chapter).

In order to be accepted by the MW projector (see MWFunctions chapter), the
lambda must have the following signature:

.. code-block:: cpp

    auto f = [] (const double *r) -> double;

e.i. it must take a ``double`` pointer defining a set of Cartesian coordinates,
and return a ``double``. For instance, the electrostatic potential from a point
nuclear charge :math:`Z` (in atomic units) is

.. math:: f(r) = \frac{Z}{r}

which can be written as the lambda function

.. code-block:: cpp

    double Z = 1.0;                                 // Hydrogen nuclear charge
    auto f = [Z] (const double *r) -> double {
        double R = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        return Z/R;
    }

Note that the function signature must be *exactly* as given above, which means
that any additional arguments (such as :math:`Z` in this case) must be given in
the capture list (square brackets), see e.g. `cppreference.com 
<http://en.cppreference.com/w/cpp/language/lambda>`_ for more
details on lambda functions and how to use the capture list.

