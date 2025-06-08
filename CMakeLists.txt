cmake_minimum_required(VERSION 3.20.3)

project(FEM)

# PkgConfig 
find_package(PkgConfig)
if (PKG_CONFIG_FOUND)
    pkg_check_modules(PETSC PETSc)
endif()
#pkg_check_modules(PETSC PETSc)
if (PETSC_FOUND)
    list(APPEND COMPILE_OPTIONS ${PETSC_CFLAGS})

    include_directories(${PETSC_INCLUDE_DIRS})
    
    set(LINK_FLAGS "${LINK_FLAGS} ${PETSC_LDFLAGS}")   
    
    list(APPEND LIBRARIES ${PETSC_LINK_LIBRARIES})   
    
    set(CMAKE_REQUIRED_FLAGS ${PETSC_CFLAGS})   
    
    set(CMAKE_REQUIRED_INCLUDES "${PETSC_INCLUDE_DIRS}")
endif()

MESSAGE(STATUS "PETSC_INCLUDE_DIRS: ${PETSC_INCLUDE_DIRS}")
MESSAGE(STATUS "PETSC_LIB: ${PETSC_LINK_LIBRARIES}")
MESSAGE(STATUS "PETSC_PATH: ${PETSC_PATH}")

set(CMAKE_CXX_FLAGS "-std=c++11 -O3")

#set(CMAKE_CXX_COMPILER "/Users/dual/Documents/projects/petsc/arch-darwin-cxx-debug/bin/mpicxx")
#add_compile_options(-Wl,-rpath,/Users/dual/Documents/projects/petsc/arch-darwin-cxx-debug/lib -L/Users/dual/Documents/projects/petsc/arch-darwin-cxx-debug/lib -lpetsc -lf2clapack -lf2cblas -lc++)

file(GLOB SOURCES "src/*")
add_executable(FEM3D ${SOURCES})

include_directories(
    #/Users/dual/Documents/projects/petsc/arch-darwin-cxx-debug/include
    /Users/dual/Documents/projects/petsc/arch-darwin-cxx-opt/include
)

target_include_directories(FEM3D PUBLIC ${PETSC_INCLUDE_DIRS})
#TARGET_LINK_LIBRARIES(FEM ${PETSC_LINK_LIBRARIES})
#TARGET_LINK_LIBRARIES(FEM /Users/dual/Documents/projects/petsc/arch-darwin-cxx-debug/lib/libmpi.dylib)
#TARGET_LINK_LIBRARIES(FEM /Users/dual/Documents/projects/petsc/arch-darwin-cxx-debug/lib/libpetsc.dylib)
TARGET_LINK_LIBRARIES(FEM3D /Users/dual/Documents/projects/petsc/arch-darwin-cxx-opt/lib/libmpi.dylib)
TARGET_LINK_LIBRARIES(FEM3D /Users/dual/Documents/projects/petsc/arch-darwin-cxx-opt/lib/libpetsc.dylib)

install(TARGETS FEM3D DESTINATION bin)


