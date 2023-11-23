add_library(modules "")

target_sources(modules
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/modules.f90
    ${CMAKE_CURRENT_LIST_DIR}/io.f90
    ${CMAKE_CURRENT_LIST_DIR}/godunov.f90
    ${CMAKE_CURRENT_LIST_DIR}/vanleer.f90
  )

