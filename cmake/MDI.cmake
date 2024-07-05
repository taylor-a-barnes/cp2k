# download/ build MDI library
# always build static library with -fpic
# support cross-compilation and ninja-build
include(ExternalProject)
ExternalProject_Add(mdi_build
  URL     "https://github.com/MolSSI-MDI/MDI_Library/archive/v1.4.28.tar.gz"
  URL_MD5 "ba02e0376267bcaa95e1ebbb9673de74"
  CMAKE_ARGS ${CMAKE_REQUEST_PIC}
  -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
  -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
  -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
  -DCMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}
  -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
  -Dlanguage=Fortran
  -Dlibtype=SHARED
  -Dmpi=OFF
  -Dpython_plugins=OFF
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_BYPRODUCTS "<BINARY_DIR>/MDI_Library/libmdi.a"
  )

# where is the compiled library?
ExternalProject_get_property(mdi_build BINARY_DIR)
set(MDI_BINARY_DIR "${BINARY_DIR}/MDI_Library")
# workaround for older CMake versions
file(MAKE_DIRECTORY ${MDI_BINARY_DIR})

# check if found
#    get_target_property(
#      CP2K_SCALAPACK_LINK_LIBRARIES cp2k::BLAS::SCI::scalapack_link
#      INTERFACE_LINK_LIBRARIES)
#find_package_handle_standard_args(MDI
#                                  REQUIRED_VARS CP2K_MDI_LINK_LIBRARIES)
set(CP2K_MDI_LINK_LIBRARIES "${MDI_BINARY_DIR}/libmdi.so")

# add target to link against
if(NOT TARGET cp2k::MDI::mdi)
  add_library(cp2k::MDI::mdi INTERFACE IMPORTED)
endif()

# create imported target for the MDI library
#add_library(mdi-lib UNKNOWN IMPORTED)
add_dependencies(cp2k::MDI::mdi mdi_build)
set_target_properties(cp2k::MDI::mdi PROPERTIES
  IMPORTED_LOCATION "${MDI_BINARY_DIR}/libmdi.a"
  INTERFACE_INCLUDE_DIRECTORIES ${MDI_BINARY_DIR}
  )

set_property(TARGET cp2k::MDI::mdi
             PROPERTY INTERFACE_LINK_LIBRARIES ${CP2K_MDI_LINK_LIBRARIES})
mark_as_advanced(CP2K_MDI_LINK_LIBRARIES)
