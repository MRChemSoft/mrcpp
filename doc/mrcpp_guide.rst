GaussFunc
---------

A very important analytic function is the Cartesian Gaussian

.. math:: f(r) = \alpha (x-x_0)^a (y-y_0)^b (z-z_0)^c e^{-\beta \|r-r_0\|^2}

which is initialized in the following way (in 3D)

.. code-block:: cpp

    double alpha, beta;
    double pos[3] = {x_0, y_0, z_0};
    int pow[3] = {a, b, c};
    GaussFunc<3> f_func(beta, alpha, pos, pow);

GaussPoly
---------

``GaussPoly`` is a generalization of the ``GaussFunc``, where there is an
arbitrary polynomial in front of the exponential

.. math:: f(r) = \alpha P(r-r_0) e^{-\beta \|r-r_0\|^2}

GaussExp
--------

A ``GaussExp`` is a collection of Gaussians in the form

.. math:: G(r) = \sum_i c_i g_i(r)

where :math:`g_i` can be either ``GaussFunc`` or ``GaussPoly``

.. math:: g_i(r) =  \alpha_i P_i(r-r_i)e^{-\beta_i\|r-r_i\|^2}


TreeCalculator
--------------

This class operates on the node level, computing MW coefficients based on the
proper input data (analytic functions in the case of projection,
``FunctionTrees`` in the case of operators). The ``TreeCalculator`` is hidden within the
``TreeBuilder``, and is not part of its interface. There is one calculator for each
of the MW-types of ``TreeBuilder``:

* ProjectionCalculator
* AdditionCalculator
* MultiplicationCalculator
* OperationCalculator

TreeAdaptor
-----------

Like the ``TreeCalculator``, this class operates on the node level, but instead of
computing coefficients, it decides whether each node needs to be split into
:math:`2^D` children nodes. There can be different reasons for splitting nodes,
the most important being to reduce the wavelet norm of the representation.
There are three different ``TreeAdaptors``:

* WaveletAdaptor
* AnalyticAdaptor
* CopyAdaptor

where the ``WaveletAdaptor`` tests the wavelet norm, the
``AnalyticAdaptor`` use some known information of an analytic function, and the
``CopyAdaptor`` will copy the node structure of another tree.


Advanced initialization
-----------------------

The default projector (``prec`` and ``max_iter`` negative) will
simply do the projection on the initial grid with no further grid refinement.
By specifying a positive ``prec`` the grid will automatically be adapted to
represent the function to the given precision, based on the wavelet norm of
the representation. You can also allow the grid to be refined only a certain number
of iterations beyond the initial grid by specifying a positive ``max_iter``
(this will of course not guarantee the accuracy of the representation).

