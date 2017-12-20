
-------
Printer
-------

MRCPP comes with its own printer class which can be used by the application
program. These features are available by including:

.. code-block:: cpp

    #include "MRCPP/Printer"


Available functions
-------------------

All Printer methods (except getters and setters) takes an integer print level
as first argument. When initialized with a given print level, only print
statements with a *lower* print level will be displayed. All internal printing
in MRCPP is at print level 10 or higher, so there is some flexibility left for
adjusting the print volume in the application program. By default the print
level is negative, which means that **no** output will be printed unless the
printer is initialized.

init
  Initialize the MRCPP printer with print level and options for parallel
  printing.

printEnvironment
  Print information on how MRCPP was built. Version number and git revision,
  which linear algebra library is used, and parallel options.

printSeparator
  Print a full line of a single character.

printHeader
  Print a headline.

printFooter
  Print a footline, including a timer.

printDouble
  Print a variable.

printTree
  Print a tree size and time.

printTime
  Print a timer.

setPrecision
  Set new precision for floating point values. Returns the *old* precision.

setPrintLevel
  Set new print level. Returns the *old* print level.

getPrecision
  Get the current precision.

getPrintLevel
  Get the current print level.


Available macros
----------------

printout
  Print a string with a given print level, no newline.

println
  Print a string with a given print level, with newline.


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

    Printer::printHeader(0, "Projecting function");

    Timer timer;                            // Start timer
    project(prec, tree, func);              // Project function
    double integral = tree.integrate();     // Integrate function
    timer.stop();                           // Stop timer

    Printer::printDouble(0, "The function integrates to", integral);
    Printer::printFooter(0, timer);


This will produce the following output::

    ============================================================
                         Projecting function
    ------------------------------------------------------------
     The function integrates to              9.999999999992e-01
    ------------------------------------------------------------
                     Wall time: 1.15362e-01 sec
    ============================================================
    

When running an application program in MPI parallel there are three different
ways of handling printed output:

* Only master rank prints to screen (stdout)
* All ranks prints to screen (stdout)
* All ranks prints to individual files

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

