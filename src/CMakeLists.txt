add_executable(prog main.f90)

add_subdirectory(modules)

target_link_libraries(prog 
	PRIVATE
	  modules
	)


