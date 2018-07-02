-----------
MWFunctions
-----------

Everything that is discussed in the following chapter is available to the
application program by including:

.. code-block:: cpp

    #include "MRCPP/MWFunctions"

Multiwavelet (MW) representations of real-valued scalar functions are in MRCPP
called ``FunctionTrees``. These are in principle available in any dimension
using the template parameter D (in practice D=1,2,3). There are several
different ways of constructing MW functions (computing the expansion
coefficients in the MW basis):

* Projection of analytic function
* Arithmetic operations
* Mapping of function values
* Application of MW operator

The first three will be discribed in this chapter, while the last one
regarding operators will be the topic of the next chapter.

The interface for constructing MW representations in MRCPP has a dual focus:
on the one hand we want a simple, intuitive way of producing adaptive numerical
approximations with guaranteed precision that does not require detailed
knowledge of the internals of the MW code and with a minimal number of
parameters that have to be set. On the other hand we want the possibility for
more detailed control of the construction and refinement of the numerical grid
where such control is possible and even necessary. In the latter case it is
important to be able to reuse the existing grids in e.g. iterative algorithms
without excessive allocation/deallocation of memory.


MultiResolution Analysis (MRA)
------------------------------

In order to combine different functions and operators in mathematical
operations, they need to be compatible. That is, they must be defined on the
same computational domain and constructed using the same polynomial basis
(order and type). This information constitutes an MRA, which needs to be
defined and passed as argument to all function and operator constructors,
and only functions and operators with compatible MRAs can be combined in
subsequent calculations. An MRA is defined in two steps, first the
computational domain is given by a ``BoundingBox`` (D is the dimension)

.. code-block:: cpp

    int n;                                          // Root scale defines box size 2^{-n}
    int l[D];                                       // Translation of first box
    int nb[D];                                      // Number of boxes
    BoundingBox<D> world(n, l, nb);

which is combined with a ``ScalingBasis`` to give an MRA

.. code-block:: cpp

    int N;                                          // Maximum refinement scale 2^{-N}
    int k;                                          // Polynomial order
    ScalingBasis basis(k);                          // Legendre or Interpolating basis
    MultiResolutionAnalysis<D> MRA(world, basis, N);

Two types of ``ScalingBasis`` are supported (``LegendreBasis`` and
``InterpolatingBasis``), and they are both available at orders
:math:`k=1,2,\dots,40` (note that some operators are constructed using
intermediates of order :math:`2k`, so in that case the maximum order available
is :math:`k=20`).


FunctionTree
------------

Constructing a full grown ``FunctionTree`` involves a number of steps,
including setting up a memory allocator, constructing root nodes according
to the given MRA, building an adaptive tree structure and computing MW
coefficients. The ``FunctionTree`` constructor does only half of these steps:

.. code-block:: cpp

    FunctionTree<D> tree(MRA);

It takes an MRA argument, which defines the computational domain and scaling
basis (these are fixed parameters that cannot be changed after construction).
The tree is initialized with a memory allocator and a set of root nodes, but
it does not compute any coefficients and the function is initially *undefined*.
An undefined ``FunctionTree`` will have a well defined tree structure (at the
very least the root nodes of the given MRA, but possibly with additional
refinement, as discussed below) and its MW coefficient will be allocated but
uninitialized, and its square norm will be negative (minus one).


Creating defined FunctionTrees
++++++++++++++++++++++++++++++

The following functions will *define* MW coefficients where there are none, and
thus *require* that the output ``FunctionTree`` is in an *undefined* state.
All functions marked with adaptive grid will use the same building algorithm:

1. Start with an initial guess for the grid
2. Compute the MW coefficients for the output function on the current grid
3. Refine the grid where necessary based on the local wavelet norm
4. Iterate points 2 and 3 until the grid is converged

With a *negative* precision argument, the grid will be *fixed*, e.i. it will
not be refined beyond the initial grid. There is also an argument to limit the
number of *extra* refinement levels beyond the initial grid, in which the
adaptive refinement will stop, even if the local precision requirement is not
met.

setZero
  Set the MW coefficients to zero, fixed grid.

copy_func
  Copy existing function into a new tree, fixed grid.

project
  Project an analytic function onto the MW basis, adaptively grid.

add
  Add existing functions, adaptive grid.

multiply
  Multiply existing functions, adaptive grid.

square
  Multiply an existing function with itself, adaptive grid.

power
  Raise an existing function to a given power, adaptive grid.

map
  Map an existing function through an analytic function, adaptive grid.



Creating undefined FunctionTrees
++++++++++++++++++++++++++++++++

The grid of a ``FunctionTree`` can also be constructed *without* computing any
MW coefficients:

build_grid
  Build an empty grid based on information from an analytic function, e.g.
  position and exponent of Gaussian, or based on the structure of another grid.

copy_grid
  Build an empty grid that is identical to that of an existing function.

clear_grid
  Clear MW coefficients of an existing function.

clear
  Clear MW coefficients and remove all grid refinement.


Changing FunctionTrees
++++++++++++++++++++++

There are also a number of in-place operations that *change* the MW
coefficients of a given defined ``FunctionTree``:

rescale
  Multiply the function with a scalar, fixed grid.

normalize
  Rescale the function by its norm, fixed grid.

add
  Add an existing function, fixed grid.

multiply
  Multiply an existing function, fixed grid.

square
  Multiply an existing function with itself, fixed grid.

power
  Raise an existing function to a given power, fixed grid.

crop
  Truncate the wavelet expansion accoring to a new precision threshold.

refine_grid
  Refine grid at most one level based on precision and the local wavelet norm,
  or on the structure of another grid. The grid changes, but the `represented`
  function remains the same.

All changing operations *require* that the ``FunctionTree`` is in a
*defined* state.


File I/O
++++++++

saveTree
  Write function to file.

loadTree
  Read function from file.


Extracting data
+++++++++++++++

Given a ``FunctionTree`` that is a *well defined* function representation, the
following data can be extracted:

getSquareNorm
  Returns the squared L2 norm of the function.

integrate
  Returns the integral of the function over the entire computational domain.

evalf
  Returns the function value in a given point.

dot
  Returns the dot product of two functions over the entire computational domain.


FunctionTreeVector
------------------

The ``FunctionTreeVector`` is simply an alias for a ``std::vector`` of tuples
containing a numerical coefficient and a ``FunctionTree`` pointer.
Elements can be appended to the vector using the ``std::make_tuple``, elements
are obtained with the ``get_func`` and ``get_coef`` functions:

.. code-block:: cpp
    
    FunctionTreeVector<D> tree_vec;
    tree_vec.push_back(std::make_tuple(2.0, &tree_a)); // Push back pointer to FunctionTree
    double coef = get_coef(tree_vec, 0);               // Get coefficient of first entry
    FunctionTree<3> &tree = get_func(tree_vec, 0);     // Get function of first entry
    clear(tree_vec, false);                            // Bool argument for tree destruction

Clearing the vector means removing all its elements, and the ``bool`` argument
tells if the elements should be properly deallocated (default ``false``).


Examples
--------

Building empty grids
++++++++++++++++++++

Sometimes it is useful to construct an empty grid based on some available
information of the function that is about to be represented. This can be e.g.
that you want to copy the grid of an existing ``FunctionTree`` or that an
analytic function has more or less known grid requirements (like Gaussians).
Sometimes it is even necessary to force the grid refinement beyond the coarsest
scales in order for the adaptive refining algorithm to detect a wavelet
"signal" that allows it to do its job properly (this happens for narrow
Gaussians where none of the initial quadrature points hits a function value
significantly different from zero).

The simplest way to build an empty grid is to copy the grid from an existing
tree (assume that ``f_tree`` has been properly built so that it contains more
than just root nodes)

.. code-block:: cpp

    mrcpp::FunctionTree<D> f_tree(MRA);                     // Input tree
    mrcpp::FunctionTree<D> g_tree(MRA);                     // Output tree

    mrcpp::project(prec, f_tree, f_func);                   // Build adaptive grid for f_tree
    mrcpp::copy_grid(g_tree, f_tree);                       // Copy grid from f_tree to g_tree

Passing an analytic function as argument to the generator will build a grid
based on some predefined information of the function (if there is any,
otherwise it will do nothing)

.. code-block:: cpp

    mrcpp::RepresentableFunction<D> func;                   // Analytic function
    mrcpp::FunctionTree<D> tree(MRA);                       // Output tree
    mrcpp::build_grid(tree, func);                          // Build grid based on f_func

The lambda analytic functions do `not` provide such information, this must be
explicitly implemented as a ``RepresentableFunction`` sub-class (see MRCPP
programmer's guide for details).

Actually, the effect of the ``build_grid`` is to *extend* the existing grid
with any missing nodes relative to the input. There is also a version of
``build_grid`` taking a ``FunctionTree`` argument. Its effect is very similar to the
``copy_grid`` above, with the only difference that now the output grid is
*extended* with the missing nodes (e.i. the nodes that are already there are
*not* removed first). This means that we can build the union of two grids by
successive applications of ``build_grid``

.. code-block:: cpp

    mrcpp::FunctionTree<D> f_tree(MRA);             // Construct empty grid of root nodes
    mrcpp::build_grid(f_tree, g_tree);              // Extend f with missing nodes relative to g
    mrcpp::build_grid(f_tree, h_tree);              // Extend f with missing nodes relative to h

In contrast, doing the same with ``copy_grid`` would clear the ``f_tree`` grid in
between, and you would *only* get a (identical) copy of the last ``h_tree`` grid,
with no memory of the ``g_tree`` grid that was once there. One can also make the
grids of two functions equal to their union

.. code-block:: cpp

    mrcpp::build_grid(f_tree, g_tree);              // Extend f with missing nodes relative to g
    mrcpp::build_grid(g_tree, f_tree);              // Extend g with missing nodes relative to f

The union grid of several trees can be constructed in one go using a
``FunctionTreeVector``

.. code-block:: cpp

    mrcpp::FunctionTreeVector<D> inp_vec;
    inp_vec.push_back(std::make_tuple(1.0, tree_1));
    inp_vec.push_back(std::make_tuple(1.0, tree_2));
    inp_vec.push_back(std::make_tuple(1.0, tree_3));

    mrcpp::FunctionTree<D> f_tree(MRA);
    mrcpp::build_grid(f_tree, inp_vec);


Projection
++++++++++

The ``project`` function takes an analytic D-dimensional scalar function (which
can be defined as a lambda function or one of the explicitly implemented
sub-classes of the ``RepresentableFunction`` base class in MRCPP) and projects
it with the given precision onto the MRA defined by the ``FunctionTree``.
E.g. a unit charge Gaussian is projected in the following way (the MRA must
be initialized as above)

.. code-block:: cpp

    // Defining an analytic function
    double beta = 10.0;
    double alpha = std::pow(beta/pi, 3.0/2.0);
    auto func = [alpha, beta] (const double *r) -> double {
        double R = std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        return alpha*std::exp(-beta*R*R);
    }

    double prec = 1.0e-5;
    FunctionTree<3> tree(MRA);
    mrcpp::project(prec, tree, func);

This projection will start at the default initial grid (only the root nodes of
the given MRA), and adaptively build the full grid. Alternatively, the grid can
be estimated *a priori* if the analytical function has some known features, such
as for Gaussians:


.. code-block:: cpp

    double prec;                                            // Precision of the projection
    int max_iter;                                           // Maximum levels of refinement

    mrcpp::GaussFunc<D> func;                               // Analytic Gaussian function
    mrcpp::FunctionTree<D> tree(MRA);                       // Output tree

    mrcpp::build_grid(tree, func);                          // Empty grid from analytic function
    mrcpp::project(prec, tree, func, max_iter);             // Starts projecting from given grid

This will first produce an empty grid suited for representing the analytic
function ``func`` (this is meant as a way to make sure that the projection
starts on a grid where the function is actually visible, as for very narrow
Gaussians, it's `not` meant to be a good approximation of the final grid) and
then perform the projection on the given numerical grid. With a negative
``prec`` (or ``max_iter = 0``) the projection will be performed strictly on the
given initial grid, with no further refinements.


Addition
++++++++

Arithmetic operations in the MW representation are performed using the
``FunctionTreeVector``, and the general sum :math:`f = \sum_i c_i f_i(x)`
is done in the following way

.. code-block:: cpp

    double a, b, c;                                         // Addition parameters
    mrcpp::FunctionTree<D> a_tree(MRA);                     // Input function
    mrcpp::FunctionTree<D> b_tree(MRA);                     // Input function
    mrcpp::FunctionTree<D> c_tree(MRA);                     // Input function

    mrcpp::FunctionTreeVector<D> inp_vec;                   // Vector to hold input functions
    inp_vec.push_back(std::make_tuple(a, &a_tree));         // Append to vector
    inp_vec.push_back(std::make_tuple(b, &b_tree));         // Append to vector
    inp_vec.push_back(std::make_tuple(c, &c_tree));         // Append to vector

    mrcpp::FunctionTree<D> f_tree(MRA);                     // Output function
    mrcpp::add(prec, f_tree, inp_vec);                      // Adaptive addition

The default initial grid is again only the root nodes, and a positive ``prec``
is required to build an adaptive tree structure for the result. The special
case of adding two functions can be done directly without initializing a
``FunctionTreeVector``

.. code-block:: cpp

    mrcpp::FunctionTree<D> f_tree(MRA);
    mrcpp::add(prec, f_tree, a, a_tree, b, b_tree);

Addition of two functions is usually done on their (fixed) union grid

.. code-block:: cpp

    mrcpp::FunctionTree<D> f_tree(MRA);                     // Construct empty root grid
    mrcpp::build_grid(f_tree, a_tree);                      // Copy grid of g
    mrcpp::build_grid(f_tree, b_tree);                      // Copy grid of h
    mrcpp::add(-1.0, f_tree, a, a_tree, b, b_tree);         // Add functions on fixed grid

Note that in the case of addition there is no extra information to be gained
by going beyond the finest refinement levels of the input functions, so the
union grid summation is simply the best you can do, and adding a positive
``prec`` will not make a difference. There are situations where you want to
use a `smaller` grid, though, e.g. when performing a unitary transformation
among a set of ``FunctionTrees``. In this case you usually don't want to
construct `all` the output functions on the union grid of `all` the input
functions, and this can be done by adding the functions adaptively starting
from root nodes.

If you have a summation over several functions but want to perform the
addition on the grid given by the `first` input function, you first copy the
wanted grid and then perform the operation on that grid

.. code-block:: cpp

    mrcpp::FunctionTreeVector<D> inp_vec;
    inp_vec.push_back(std::make_tuple(a, a_tree));
    inp_vec.push_back(std::make_tuple(b, b_tree));
    inp_vec.push_back(std::make_tuple(c, c_tree));

    mrcpp::FunctionTree<D> f_tree(MRA);                     // Construct empty root grid
    mrcpp::copy_grid(f_tree, get_func(inp_vec, 0));         // Copy grid of first input function
    mrcpp::add(-1.0, f_tree, inp_vec);                      // Add functions on fixed grid

Here you can of course also add a positive ``prec`` to the addition and the
resulting function will be built adaptively starting from the given initial
grid.


Multiplication
++++++++++++++

The multiplication follows the exact same syntax as the addition, where the
product :math:`f = \prod_i c_i f_i(x)` is done in the following way

.. code-block:: cpp

    double a, b, c;                                         // Multiplication parameters
    mrcpp::FunctionTree<D> a_tree(MRA);                     // Input function
    mrcpp::FunctionTree<D> b_tree(MRA);                     // Input function
    mrcpp::FunctionTree<D> c_tree(MRA);                     // Input function

    mrcpp::FunctionTreeVector<D> inp_vec;                   // Vector to hold input functions
    inp_vec.push_back(std::make_tuple(a, &a_tree));         // Append to vector
    inp_vec.push_back(std::make_tuple(b, &b_tree));         // Append to vector
    inp_vec.push_back(std::make_tuple(c, &c_tree));         // Append to vector

    mrcpp::FunctionTree<D> f_tree(MRA);                     // Output function
    mrcpp::multipy(prec, f_tree, inp_vec);                  // Adaptive multiplication

In the special case of multiplying two functions the coefficients are collected
into one argument

.. code-block:: cpp

    mrcpp::FunctionTree<D> f_tree(MRA);
    mrcpp::multiply(prec, f_tree, a*b, a_tree, b_tree);

For multiplications, there might be a loss of accuracy if
the product is restricted to the union grid. The reason for this is that the
product will contain signals of higher frequency than each of the input
functions, which require a higher grid refinement for accurate representation.
By specifying a positive ``prec`` you will allow the grid to adapt to the higher
frequencies, but it is usually a good idea to restrict to one extra refinement
level beyond the union grid (by setting ``max_iter=1``) as the grids are not
guaranteed to converge for such local operations (like arithmetics, derivatives
and function mappings)

.. code-block:: cpp

    mrcpp::FunctionTree<D> f_tree(MRA);                     // Construct empty root grid
    mrcpp::build_grid(f_tree, a_tree);                      // Copy grid of a
    mrcpp::build_grid(f_tree, b_tree);                      // Copy grid of b
    mrcpp::multiply(prec, f_tree, a*b, a_tree, b_tree, 1);  // Allow 1 extra refinement


Re-using grids
++++++++++++++

Given a ``FunctionTree`` that is a valid function representation, we can clear
its MW expansion coefficients as well as its grid refinement

.. code-block:: cpp

    mrcpp::FunctionTree<D> tree(MRA);                       // tree is an undefined function
    mrcpp::project(prec, tree, f_func);                     // tree represents analytic function f
    tree.clear();                                           // tree is an undefined function 
    mrcpp::project(prec, tree, f_func);                     // tree represents analytic function g

This action will leave the ``FunctionTree`` in the same state as after
construction (undefined function, only root nodes), and its coefficients can
now be re-computed.

In certain situations it might be desireable to separate the actions of
computing MW coefficients and refining the grid. For this we can use the
``refine_grid``, which will adaptively refine the grid one level (based on
the wavelet norm and the given precision) and project the existing function
representation onto the new finer grid

.. code-block:: cpp

    mrcpp::refine_grid(tree, prec);

E.i., this will *not* change the function that is represented in ``tree``, but
it *might* increase its grid size. The same effect can be made using another
``FunctionTree`` argument instead of the precision parameter

.. code-block:: cpp

    mrcpp::refine_grid(tree_out, tree_in);

which will *extend* the grid of ``tree_out`` in the same way as ``build_grid``
as shown above, but it will *keep* the function representation in ``tree_out``.

This functionality can be combined with ``clear_grid`` to make a "manual"
adaptive building algorithm. One example where this might be useful is in
iterative algorithms where you want to fix the grid size for all calculations
within one cycle and then relax the grid in the end in preparation for the next
iteration. The following is equivalent to the adaptive projection above
(``refine_grid`` returns the number of new nodes that were created in the
process)

.. code-block:: cpp

    int n_nodes = 1;
    while (n_nodes > 0) {
        mrcpp::project(-1.0, tree, func);                   // Project f on fixed grid
        n_nodes = mrcpp::refine_grid(tree, prec);           // Refine grid based on prec
        if (n_nodes > 0) mrcpp::clear_grid(tree);           // Clear grid for next iteration
    }


