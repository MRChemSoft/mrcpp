#.rst:
#
# Enables creation of static/shared library.
#
# autocmake.yml configuration::
#
#   docopt: "--static Create only the static library [default: False]."
#   define: "'-DSTATIC_LIBRARY_ONLY={0}'.format(arguments['--static'])"

option_with_print(STATIC_LIBRARY_ONLY "Create the static library only" OFF)
option_with_print(SHARED_LIBRARY_ONLY "Create the shared library only" OFF)
