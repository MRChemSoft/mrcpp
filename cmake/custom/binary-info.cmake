# This file is a CMake script, not supposed to be included anywhere else.
# Based on: https://github.com/dev-cafe/cmake-cookbook/tree/master/chapter-06/recipe-07

# in case Git is not available, we default to "unknown"
set(_git_last_commit_hash "unknown")
set(_git_last_commit_author "unknown")
set(_git_last_commit_date "unknown")
set(_git_branch "unknown")
set(_git_describe "unknown")

# find Git and if available set GIT_HASH variable
find_package(Git QUIET)
if(GIT_FOUND)
  # Get hash of last commit
  execute_process(
    COMMAND ${GIT_EXECUTABLE} --no-pager show -s --pretty=format:%H -n 1
    OUTPUT_VARIABLE _git_last_commit_hash
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
    )
  # Get author of last commit
  execute_process(
    COMMAND ${GIT_EXECUTABLE} --no-pager show -s --pretty=format:%aN -n 1
    OUTPUT_VARIABLE _git_last_commit_author
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
    )
  # Get date of last commit
  execute_process(
    COMMAND ${GIT_EXECUTABLE} --no-pager show -s --pretty=format:%ad -n 1
    OUTPUT_VARIABLE _git_last_commit_date
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
    )
  # Get branch name
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
    OUTPUT_VARIABLE _git_branch
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
    )
  # Get commit description
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --abbrev=20 --long --always --dirty --match "NOT A TAG"
    OUTPUT_VARIABLE _git_describe
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
    )
endif()

message(STATUS "Git branch: ${_git_branch}")
message(STATUS "Git last commit hash: ${_git_describe}")
message(STATUS "Git last commit author: ${_git_last_commit_author}")
message(STATUS "Git last commit date: ${_git_last_commit_date}")

cmake_host_system_information(RESULT _hostname QUERY HOSTNAME)

execute_process(
  COMMAND whoami
  OUTPUT_VARIABLE _username
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
  )

execute_process(
  COMMAND date +"%Y-%m-%d %H:%M:%S"
  OUTPUT_VARIABLE _configuration_time
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
  )
string(REGEX REPLACE "\"" "" _configuration_time ${_configuration_time})

file(READ "${VERSION_FILE}" PROGRAM_VERSION)
string(STRIP "${PROGRAM_VERSION}" PROGRAM_VERSION)

# generate file
configure_file(
  ${INPUT_DIR}/version.h.in
  ${TARGET_DIR}/version.h
  @ONLY
  )
