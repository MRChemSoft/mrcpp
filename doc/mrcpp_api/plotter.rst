
-------
Plotter
-------

MRCPP comes with its own plotter class which can be used by the host
program to generate data files for visualization using e.g.
`gnuplot <http://www.gnuplot.info/>`_,
`paraview <http://www.paraview.org/>`_,
`blob <https://github.com/densities/blob/>`_ and
`geomview <http://www.geomview.org/>`_.
These features are available by including:

.. code-block:: cpp

    #include "MRCPP/Plotter"


.. doxygenclass:: mrcpp::Plotter
    :members:

.. NOTE::

    When plotting a ``FunctionTree``, only the *scaling* part of the
    leaf nodes will be evaluated, which means that the function
    values will not be fully accurate. This is done to allow a
    fast and ``const`` function evaluation that can be done in
    OMP parallel. If you want to include also the *final* wavelet
    corrections to your function values, you'll have to manually
    extend the MW grid by one level before plotting using
    ``mrcpp::refine_grid(tree, 1)``.


Examples
--------

A parametric line plot of a three-dimensional function along the z axis [-1, 1]:

.. code-block:: cpp

    mrcpp::FunctionTree<3> f_tree(MRA);                  // Function to be plotted

    int nPts = 1000;                                     // Number of points
    mrcpp::Coord<3> o{ 0.0, 0.0,-1.0};                   // Origin vector
    mrcpp::Coord<3> a{ 0.0, 0.0, 2.0};                   // Boundary vector

    mrcpp::Plotter<3> plot(o);                           // Plotter of 3D functions
    plot.setRange(a);                                    // Set plot range
    plot.linePlot({nPts}, f_tree, "f_tree");             // Write to file f_tree.line


A surface plot of a three-dimensional function in the x=[-2,2], y=[-1,1], z=0 plane:

.. code-block:: cpp

    int aPts = 2000;                                     // Number of points in a
    int bPts = 1000;                                     // Number of points in b
    mrcpp::Coord<3> o{-2.0,-1.0, 0.0};                   // Origin vector
    mrcpp::Coord<3> a{ 4.0, 0.0, 0.0};                   // Boundary vector A
    mrcpp::Coord<3> b{ 0.0, 2.0, 0.0};                   // Boundary vector B

    mrcpp::Plotter<3> plot(o);                           // Plotter of 3D functions
    plot.setRange(a, b);                                 // Set plot range
    plot.surfPlot({aPts, bPts}, f_tree, "f_tree");       // Write to file f_tree.surf

A cube plot of a three-dimensional function in the volume x=[-2,2], y=[-1,1], z=[0,2]:

.. code-block:: cpp

    int aPts = 200;                                      // Number of points in a
    int bPts = 100;                                      // Number of points in b
    int cPts = 100;                                      // Number of points in c
    mrcpp::Coord<3> o{-2.0,-1.0, 0.0};                   // Origin vector
    mrcpp::Coord<3> a{ 4.0, 0.0, 0.0};                   // Boundary vector A
    mrcpp::Coord<3> b{ 0.0, 2.0, 0.0};                   // Boundary vector B
    mrcpp::Coord<3> b{ 0.0, 0.0, 2.0};                   // Boundary vector C

    mrcpp::Plotter<3> plot(o);                           // Plotter of 3D functions
    plot.setRange(a, b, c);                              // Set plot range
    plot.cubePlot({aPts, bPts, cPts}, f_tree, "f_tree"); // Write to file f_tree.cube

A grid plot of a three-dimensional FunctionTree:

.. code-block:: cpp

    mrcpp::Plotter<3> plot;                              // Plotter of 3D functions
    plot.gridPlot(f_tree, "f_tree");                     // Write to file f_tree.grid
