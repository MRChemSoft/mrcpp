-------
Printer
-------

MRCPP comes with a printer class to handle standard output:

.. code-block:: cpp

    #include "MRCPP/Printer"

The main purpose of this class is to provide (or suppress) any internal printing
in MRCPP routines that might be useful for debugging. Also, it provides a sane
printing environment for parallel computations that can also be used by the
host program. By using the printing routines of this class, as opposed to the
standard ``std::cout``, only the master thread in a OpenMP region will provide
any output while all other threads remain silent. Similarly, when running a
host program in MPI parallel, the ``mrcpp::Printer`` provides three different
options for handling printed output (see examples below):

* Only master rank prints to screen (stdout)
* All ranks prints to screen (stdout)
* All ranks prints to individual files

If you want only the master rank to print to an output file, this can be
achieved by redirecting the output from the first option to a file
(``./program >file.out``).


.. doxygenclass:: mrcpp::Printer
    :members:

Functions
---------

Some convenience functions for printing output is provided within the
``mrcpp::print`` namespace. These functions use the data of the
``mrcpp::Printer`` class to provide pretty output of a few standard data types.

.. doxygenfunction:: mrcpp::print::environment
.. doxygenfunction:: mrcpp::print::separator
.. doxygenfunction:: mrcpp::print::header
.. doxygenfunction:: mrcpp::print::footer
.. doxygenfunction:: mrcpp::print::tree(int, const std::string&, const MWTree<D>&, const Timer&)
.. doxygenfunction:: mrcpp::print::tree(int, const std::string&, int, int, double)
.. doxygenfunction:: mrcpp::print::time
.. doxygenfunction:: mrcpp::print::memory

Macros
------

The following macros should replace the regular calls to ``std::cout``:

.. doxygendefine:: println
.. doxygendefine:: printout

The following macros will print a message along with information on where you
are in the code (file name, line number and function name). Only macros that
end with ``_ABORT`` will kill the program, all other will continue to run after
the message is printed:

.. doxygendefine:: MSG_INFO
.. doxygendefine:: MSG_WARN
.. doxygendefine:: MSG_ERROR
.. doxygendefine:: MSG_ABORT
.. doxygendefine:: INVALID_ARG_ABORT
.. doxygendefine:: NOT_IMPLEMENTED_ABORT
.. doxygendefine:: NOT_REACHED_ABORT
.. doxygendefine:: NEEDS_TESTING
.. doxygendefine:: NEEDS_FIX

Examples
--------

Using the print level to adjust the amount of output:

.. code-block:: cpp

    int level = 10;
    mrcpp::Printer::init(level);            // Initialize printer with printlevel 10

    println( 0, "This is printlevel  0");   // This will be displayed at printlevel 10
    println(10, "This is printlevel 10");   // This will be displayed at printlevel 10
    println(11, "This is printlevel 11");   // This will NOT be displayed at printlevel 10


Using headers and footers to get pretty output:

.. code-block:: cpp

    using namespace mrcpp;

    Timer timer;                            // Start timer
    project(prec, tree, func);              // Project function
    double integral = tree.integrate();     // Integrate function
    timer.stop();                           // Stop timer
    
    print::header(0, "Projecting analytic function");
    print::tree(0, "Projected function", tree, timer);
    print::value(0, "Integrated function", integral, "(au)");
    print::footer(0, timer);


This will produce the following output::

    ============================================================
                    Projecting analytic function
    ------------------------------------------------------------
     Projected function         520 nds       16 MB    0.09 sec
     Integrated function               (au)  9.999999999992e-01
    ------------------------------------------------------------
                      Wall time: 9.32703e-02 sec
    ============================================================
    

As mentioned above, when running in MPI parallel there are three different ways
of handling printed output (master to stdout, all to stdout or all to files).
These can be chosen by adding appropriate arguments to ``init``. The default
setting will in a parallel environment have all MPI ranks printing to screen,
but by adding MPI info to the printer, we can separate the output of the
different ranks:

.. code-block:: cpp


    int level = 10;
    int wrank, wsize;
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);  // Get my rank
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);  // Get total number of ranks

    // All ranks will print to screen
    mrcpp::Printer::init(level);

    // Only master rank will print to screen
    mrcpp::Printer::init(level, wrank, wsize);

    // All ranks will print to separate files called filename-<rank>.out
    mrcpp::Printer::init(level, wrank, wsize, "filename");

