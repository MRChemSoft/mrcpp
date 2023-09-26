.. MRChem documentation master file, created by
   sphinx-quickstart on Tue Jan 26 15:03:29 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

----------
Clang-tidy
----------

To ensure modern coding conventions are followed developers
are encouraged to run clang-tidy on the code. Ensure clang-tidy
is installed. Then to display available checkers run::

    $ clang-tidy --list-checks -checks='*' | grep "modernize"

This will generate a list looking like this::

    $ modernize-avoid-bind
    $ modernize-deprecated-headers
    $ modernize-loop-convert
    $ modernize-make-shared
    $ modernize-make-unique
    $ modernize-pass-by-value
    $ modernize-raw-string-literal
    $ modernize-redundant-void-arg
    $ modernize-replace-auto-ptr
    $ modernize-replace-random-shuffle
    $ modernize-return-braced-init-list
    $ modernize-shrink-to-fit
    $ modernize-unary-static-assert
    $ modernize-use-auto
    $ modernize-use-bool-literals
    $ modernize-use-default-member-init
    $ modernize-use-emplace
    $ modernize-use-equals-default
    $ modernize-use-equals-delete
    $ modernize-use-noexcept
    $ modernize-use-nullptr
    $ modernize-use-override
    $ modernize-use-transparent-functors
    $ modernize-use-using


To run any of these modernization's on the code, go to your build directory.
From there run the command::

    $ run-clang-tidy -header-filter='.*' -checks='-*,modernize-your-modernization' -fix
