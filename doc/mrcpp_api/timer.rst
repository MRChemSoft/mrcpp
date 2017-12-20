
-----
Timer
-----

MRCPP comes with its own timer class which can be used by the application
program:

.. code-block:: cpp

    #include "MRCPP/Timer"


Available functions
-------------------

start
  Start clock from zero.

resume
  Resume clock from previous stopped time.

stop
  Stop clock.

getWallTime
  Get previous stopped time (in sec).


Examples
--------

The timer records wall (human) time, not CPU user time. The clock will by
default start immediately after construction, and should be explicitly stopped
before reading the time:

.. code-block:: cpp

    mrcpp::Timer timer;                         // This will start the timer
    mrcpp::project(prec, tree, func);           // Do some work
    timer.stop();                               // This will stop the timer
    double t = timer.getWallTime();             // Get time spent projecting


The timer can also be started explicitly at a later stage *after* construction,
and can be combined with the ``mrcpp::Printer`` class for convenient printing:

.. code-block:: cpp

    mrcpp::Timer timer(false);                  // This will not start the timer
    timer.start();                              // This will start the timer
    mrcpp::project(prec, tree, func);           // Do some work
    timer.stop();                               // This will stop the timer
    mrcpp::Printer::printTimer(0, "Time spent projecting function", timer);

