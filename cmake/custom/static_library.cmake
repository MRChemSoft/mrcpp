#.rst:
#
# Enables creation of static/shared library.
#
# autocmake.yml configuration::
#
#   docopt: "--static Create only the static library [default: False]."
#   define: "'-DSTATIC_LIBRARY_ONLY={0}'.format(arguments['--static'])"


option_with_print(
  NAME
    STATIC_LIBRARY_ONLY
  MESSAGE
    "Create the static library only"
  DEFAULT
    OFF
  )
option_with_print(
  NAME
    SHARED_LIBRARY_ONLY
  MESSAGE
    "Create the shared library only"
  DEFAULT
    OFF
  )
