#.rst:
#
# Enables creation of static/shared library.
#
# autocmake.yml configuration::
#
#   docopt: "--static Create only the static library [default: False]."
#   define: "'-DBUILD_STATIC_LIBS={0}'.format(arguments['--static'])"

option_with_print(BUILD_STATIC_LIBS "Create the static library" OFF)
# By default build shared library
set(BUILD_SHARED_LIBS ON)
if(BUILD_STATIC_LIBS)
  set(BUILD_SHARED_LIBS OFF)
endif()
