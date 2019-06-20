
-------
Plotter
-------

MRCPP comes with its own plotter class which can be used by the application
program to generate data files for visualization using e.g.
`gnuplot <http://www.gnuplot.info/>`_,
`blob <https://github.com/densities/blob/>`_ and
`geomview <http://www.geomview.org/>`_.
These features are available by including:

.. code-block:: cpp

    #include "MRCPP/Plotter"

This class will generate an equidistant grid in one (line), two (surf)
or three (cube) dimensions, and subsequently evaluate the function on
this grid.

The grid is generated from the vectors A, B and C in relation to the origin O:
 - a linePlot will plot the line spanned by A, starting from O
 - a surfPlot will plot the area spanned by A and B, starting from O
 - a cubePlot will plot the volume spanned by A, B and C, starting from O

The vectors A, B and C does not necessarily have to be orthogonal.


Available functions
-------------------

setOrigin
  Set the point of origin O for the plot.

setRange
  Set boundary vectors A, B and C for the plot.

linePlot
  Plot a D-dimensional function along the vector A, starting from O.

surfPlot
  Plot a D-dimensional function on the area defined by the two boundary
  vectors A and B, starting from O.

cubePlot
  Plot a D-dimensional function on the volume defined by the three boundary
  vectors A, B and C, starting from O.

gridPlot
  Produce a file for visualizing the multiresolution grid using geomview.


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
    plot.surfPlot({aPts, bPts}, f_tree, "f_tree");       // Write to file f_tree.line

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
    plot.surfPlot({aPts, bPts, cPts}, f_tree, "f_tree"); // Write to file f_tree.line

A grid plot of a three-dimensional FunctionTree:

.. code-block:: cpp

    mrcpp::Plotter<3> plot;                              // Plotter of 3D functions
    plot.gridPlot(f_tree, "f_tree");                     // Write to file f_tree.line
