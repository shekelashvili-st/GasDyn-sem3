add_executable(test_godunov test_godunov.f90)
target_link_libraries(test_godunov modules)

add_test(
  NAME
    test_godunov
  COMMAND
    $<TARGET_FILE:test_godunov>
  )