set(GIT_REVISION "undefined")
find_package(Git QUIET)
if(GIT_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-list --abbrev-commit --max-count=1 HEAD
    OUTPUT_VARIABLE _git_revision
    ERROR_QUIET
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )
  if(_git_revision)
    string(STRIP ${_git_revision} GIT_REVISION)
  endif()
endif()
