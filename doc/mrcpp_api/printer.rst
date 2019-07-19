
-------
Printer
-------

MRCPP comes with a printer class to handle standard output which is available
by including:

.. code-block:: cpp

    #include "MRCPP/Printer"

The main purpose of this class is to provide (or suppress) any internal printing
in MRCPP routines that might be useful for debugging. Also, it provides a sane
printing environment for parallel computations that can also be used by the
application program. E.g. by using the printing routines of this class, as
opposed to the standard ``std::cout``, only the master thread in a OpenMP region
will provide any output while all other threads remain silent. Similarly, when
running an application program in MPI parallel the ``mrcpp::Printer`` provides
three different options for handling printed output (see examples below):

* Only master rank prints to screen (stdout)
* All ranks prints to screen (stdout)
* All ranks prints to individual files

If you want only the master rank to print to an output file, this can be
achieved by redirecting the output from the first option to a file
(``./program >file.out``).

The main functionality of the MRCPP printer is given by the two macros
``printout`` and ``println`` (see below for details) which should replace the
regular calls to ``std::cout``. In addition there are a few helper functions
and macros for convenience output.


Available functions
-------------------

All ``mrcpp::print`` functions takes an integer print level as first argument.
When the ``mrcpp::Printer`` is initialized with a given print level, only print
statements with a *lower* print level will be displayed. All internal printing
in MRCPP is at print level 10 or higher, so there is some flexibility left
(levels 0 through 9) for adjusting the print volume within the application
program.

init
  Initialize the ``mrcpp::Printer`` with print level and options for parallel
  printing. No output is provided by the `mrcpp::print` routines unless the
  ``Printer`` is initialized.

setScientific
  Use scientific notation (e.g. 1.0e02) for floating point numbers.

setFixed
  Use fixed notation (e.g. 100.0) for floating point numbers.

setPrecision
  Set new precision for floating point values, e.i. number of digits after the
  comma. Returns the *old* precision.

getPrecision
  Get the current precision for floating point values.

setPrintLevel
  Change the print level from this point onwards. Returns the *old* print level.

getPrintLevel
  Get the current print level.

setWidth
  Set the line width (in number of characters) of the helper function output.
  Default is 60 characters.

getWidth
  Get the current line width (in number of characters).


Helper functions
----------------

Some convenience functions for printing output is provided within the
``mrcpp::print`` namespace. These functions use the data of the
``mrcpp::Printer`` class to provide pretty output of a few standard data types.

print::environment
  Print information on how MRCPP was built. Version number and git revision,
  which linear algebra library is used, and parallel options.

print::separator
  Print a full line of a single character.

print::header
  Print a headline.

print::footer
  Print a footline, including a timer.

print::value
  Print a floating point value, including unit.

print::tree
  Print a tree size (nodes and memory) and time.

print::memory
  Print the current memory usage of the this process.


Available macros
----------------

The following macros should replace the regular calls to ``std::cout``:

printout
  Print a string with a given print level, no newline.

println
  Print a string with a given print level, with newline.


The following macros will print a message along with information on where you
are in the code (file name, line number and function name). Only macros that
end with ``_ABORT`` will kill the program, all other will continue to run after
the message is printed.

MSG_INFO
  Prints a "Info: " message, **no** abort.

MSG_WARN
  Prints a "Warning: " message, **no** abort.

MSG_ERROR
  Prints a "Error: " message, **no** abort.

MSG_ABORT
  Prints a "Error: " message and aborts the program.

INVALID_ARG_ABORT
  You have passed an invalid argument to a function. Aborts program.

NOT_IMPLEMENTED_ABORT
  You have reached a function that is yet to be implemeted. Aborts program.

NOT_REACHED_ABORT
  You have reached a point in the code that should not be reached, probably due
  to a bug or inconsistency. Aborts program.

NEEDS_TESTING
  You have reached an experimental part of the code, results cannot be trusted.

NEEDS_FIX
  You have hit a known bug that is yet to be fixed, results cannot be trusted.

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

