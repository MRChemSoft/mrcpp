
-------
Plotter
-------

MRCPP comes with its own plotter class which can be used by the application
program to generate data files for visualization using e.g.
`gnuplot <http://www.gnuplot.info/>`_ and
`geomview <http://www.geomview.org/>`_.
These features are available by including:

.. code-block:: cpp

    #include "MRCPP/Plotter"


Available functions
-------------------

setNPoints
  Set number of equidistant points to be plotted within range.

setRange
  Set boundary points for the plot.

linePlot
  Plot a D-dimensional function on the line between the boundary points.

surfPlot
  Plot a D-dimensional function on a rectangular surface in the x-y plane
  defined by the boundary points (the z coordinate must be constant).

cubePlot
  Plot a D-dimensional function on a cube defined by the boundary points.

gridPlot
  Produce a file for visualizing the multiresolution grid using geomview.


Examples
--------

.. code-block:: cpp

    mrcpp::FunctionTree<3> tree(MRA);               // Function to be plotted

    int nPts = 1000;                                // Number of points
    double a[3] = {-1.0,-1.0, 0.0};                 // Start point of plot
    double b[3] = { 1.0, 1.0, 0.0};                 // End point of plot

    mrcpp::Plotter<3> plot;                         // Plotter of 3D functions
    plot.setNPoints(nPts);                          // Set number of points
    plot.setRange(a, b);                            // Set plot range

    plot.linePlot(tree, "func_tree");               // Write to file func_tree.line
    plot.surfPlot(tree, "func_tree");               // Write to file func_tree.surf
    plot.cubePlot(tree, "func_tree");               // Write to file func_tree.cube
    plot.gridPlot(tree, "func_tree");               // Write to file func_tree.grid

