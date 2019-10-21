set(GIT_REVISION "undefined")
find_package(Git)

if(GIT_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-list --abbrev-commit --max-count=1 HEAD
    OUTPUT_VARIABLE _git_revision
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )
  string(STRIP ${_git_revision} _git_revision)
  if(_git_revision)
    set(GIT_REVISION ${_git_revision})
  endif()
endif()
