
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

FIND_PACKAGE(deal.II 8.4 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()


SET(CLEAN_UP_FILES *.tex *.vtk *.vtu *.msh *.dofs test_file.prm)

# Set the name of the project and target:
SET(TARGET "main")
#SET(TARGET "test")


FILE(GLOB TARGET_SRC "./source/*.cpp")

INCLUDE_DIRECTORIES("include")

# Declare all source files the target consists of:
SET(TARGET_SRC
	./source/main/main.cpp 				# main version
#	./tests/Poisson_test.cpp			# Mixed FEM and LDG test
#	./tests/IMEX_LDG_test.cpp			# Mixed FEM and LDG-IMEX test
#	./tests/DD_Poisson_test.cpp		# Coupled Mixed FEM-LDG test
#	./tests/interface_test.cpp		# Coupled interface test
	${TARGET_SRC}
	${TEST_SRC}
	)

# Usually, you will not need to modify anything beyond this point...

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})
DEAL_II_INVOKE_AUTOPILOT()
