
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

elapsed
  Get current time (in sec).


Examples
--------

The timer records wall (human) time, not CPU user time. The clock will by
default start immediately after construction, and will keep running until
explicitly stopped. The elapsed time can be evaluated while clock is running.

.. code-block:: cpp

    mrcpp::Timer timer;                         // This will start the timer
    mrcpp::project(prec, tree, func);           // Do some work
    double t = timer.elapsed();                 // Get time since clock started while still running


The timer can also be started explicitly at a later stage *after* construction,
as well as explicitly stopped after the work is done. Then the `elapsed()`
function will return the time spent between `start()` and `stop()`:

.. code-block:: cpp

    mrcpp::Timer timer(false);                  // This will not start the timer
    timer.start();                              // This will start the timer
    mrcpp::project(prec, tree, func);           // Do some work
    timer.stop();                               // This will stop the timer
    double t = timer.elapsed();                 // Get time spent between start and stop

