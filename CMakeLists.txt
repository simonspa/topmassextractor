########################################################
# cmake file for CMSPixelDecoder
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
########################################################

SET(GCC_COMPILE_FLAGS "-Wall -Wextra -g -Wno-deprecated -pedantic -Wno-long-long")
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COMPILE_FLAGS}" )

SET(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -fno-inline -fdiagnostics-show-option -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wmissing-include-dirs -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wsign-promo -Wstrict-null-sentinel -Wswitch-default -Wundef" CACHE STRING "Debug options." FORCE )
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -Wall"  CACHE STRING "Relwithdebinfo options." FORCE )

# Set the correct build type and allow command line options:
IF(!CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Relwithdebinfo)
  #SET(CMAKE_BUILD_TYPE Debug)
ENDIF(!CMAKE_BUILD_TYPE)

SET(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake )

MESSAGE( STATUS "Build Type: " ${CMAKE_BUILD_TYPE} )
MESSAGE( STATUS "Installing to: " ${PROJECT_SOURCE_DIR} )

FIND_PACKAGE(ROOT REQUIRED)

INCLUDE_DIRECTORIES("${CMAKE_CURRENT_BINARY_DIR}" "${ROOT_INCLUDE_DIR}" utils extractor main)

# Preparing the build of the shared decoder library
SET(extractor_sources 
  "extractor/extractor.cc"
  "extractor/extractor_diffxs.cc"
  "extractor/extractor_yields.cc"
  "utils/plotter.cc"
  "utils/latex.cc"
  "utils/helpers.cc")

ADD_LIBRARY(extractor SHARED ${extractor_sources})
ADD_EXECUTABLE (extract main/extract.cc main/extract_yield.cc main/extract_diffxs.cc)
MESSAGE( STATUS "Building extrator library." )

TARGET_LINK_LIBRARIES(extractor ${ROOT_LIBRARIES})
TARGET_LINK_LIBRARIES (extract extractor)

INSTALL(TARGETS extractor
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

INSTALL(TARGETS extract
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

ADD_EXECUTABLE (compare tools/powheg_vs_data.C utils/plotter.cc)
TARGET_LINK_LIBRARIES(compare ${ROOT_LIBRARIES})

INSTALL(TARGETS compare
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

ADD_EXECUTABLE (comptool tools/comptool.cc utils/plotter.cc utils/helpers.cc)
TARGET_LINK_LIBRARIES(comptool ${ROOT_LIBRARIES})

INSTALL(TARGETS comptool
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
