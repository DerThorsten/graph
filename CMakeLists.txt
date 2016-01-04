cmake_minimum_required(VERSION 2.6)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)
project(graph)
include_directories(include)
file(GLOB headers
    include/andres/*.hxx
    include/andres/ilp/*.hxx
    include/andres/graph/*.hxx
    include/andres/graph/multicut/*.hxx
    include/andres/graph/multicut-lifted/*.hxx
    include/andres/graph/twocut-lifted/*.hxx
    include/andres/graph/hdf5/*.hxx
    include/andres/graph/detail/*.hxx
)
enable_testing()

##############################################################################
# Doxygen
##############################################################################
find_package(Doxygen)

##############################################################################
# OpenMP
##############################################################################
#find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

##############################################################################
# HDF5
##############################################################################
find_package(HDF5 1.8.0)

##############################################################################
# Gurobi
##############################################################################
find_package(GUROBI)
if(GUROBI_FOUND)
    include_directories(${GUROBI_INCLUDE_DIR})
endif()


##############################################################################
# Boost Python vigra vigra support 
##############################################################################
OPTION(BUILD_PYTHON "build python bindings" ON)
if(BUILD_PYTHON)
    #--------------------------------------------------------------
    # vigra
    #--------------------------------------------------------------
    FIND_PACKAGE(VIGRA REQUIRED)
    include_directories(${VIGRA_INCLUDE_DIR})

    #--------------------------------------------------------------
    # Boost
    #--------------------------------------------------------------
    find_package(Boost 1.47 REQUIRED)
    find_package(Boost COMPONENTS python REQUIRED)

    #--------------------------------------------------------------
    # Python Lib
    #--------------------------------------------------------------
    find_package (PythonLibs REQUIRED)
    include_directories(
        ${PYTHON_INCLUDE_DIRS}
        ${PYTHON_NUMPY_INCLUDE_DIR}
        ${Boost_INCLUDE_DIR}
        ${Boost_PYTHON_INCLUDE_DIR}
    )


    add_library( _graph SHARED  
        src/andres/graph/python/graph_module.cxx
        src/andres/graph/python/grid_graph.cxx
        src/andres/graph/python/graph.cxx
        src/andres/graph/python/lifted_mc_model.cxx
        src/andres/graph/python/lifted_mc.cxx
        # lifted mc solvers
        src/andres/graph/python/lifted_mc_parallel_solver.cxx
    )
    set_target_properties(_graph PROPERTIES PREFIX "")

    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/python/graph)


    add_custom_target(copy_python__init__ ALL#POST_BUILD
        COMMAND cmake -E copy_if_different 
        ${CMAKE_SOURCE_DIR}/src/andres/graph/python/module/__init__.py 
        ${CMAKE_BINARY_DIR}/python/graph/__init__.py)

    ADD_CUSTOM_COMMAND(
            TARGET _graph
            POST_BUILD
            COMMAND ${CMAKE_COMMAND} ARGS -E copy_if_different $<TARGET_FILE:_graph> ${CMAKE_BINARY_DIR}/python/graph/
            COMMENT "Copying target _graph to temporary module directory")

    #DEPENDS _graph )

    #add_custom_target(copy_python_test POST_BUILD
    #COMMAND cmake -E copy_if_different ${CMAKE_SOURCE_DIR}/test_graph.py ${CMAKE_BINARY_DIR}/test_graph.py
    #DEPENDS _graph )


    #--------------------------------------------------------------
    # Link libraries
    #--------------------------------------------------------------
    find_library(RT rt)
    target_link_libraries( _graph
        ${Boost_PYTHON_LIBRARIES}
        ${VIGRA_IMPEX_LIBRARY}
        ${VIGRA_NUMPY_CORE_LIBRARY}
    )




endif()




##############################################################################
# C++11 support
##############################################################################
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Some functionality will not be available.")
endif()

##############################################################################
# targets: graph
##############################################################################
add_executable(test-graph-graph src/andres/graph/unit-test/graph.cxx ${headers})
add_test(test-graph-graph test-graph-graph)

add_executable(test-graph-digraph src/andres/graph/unit-test/digraph.cxx ${headers})
add_test(test-graph-digraph test-graph-digraph)

add_executable(test-graph-shortest-paths src/andres/graph/unit-test/shortest-paths.cxx ${headers})
add_test(test-graph-shortest-paths test-graph-shortest-paths)

add_executable(test-graph-paths src/andres/graph/unit-test/paths.cxx ${headers})
add_test(test-graph-paths test-graph-paths)

add_executable(test-graph-components src/andres/graph/unit-test/components.cxx ${headers})
add_test(test-graph-components test-graph-components)

add_executable(test-graph-max-flow src/andres/graph/unit-test/max-flow.cxx ${headers})
add_test(test-graph-max-flow test-graph-max-flow)

if(COMPILER_SUPPORTS_CXX0X OR COMPILER_SUPPORTS_CXX11)
    add_executable(test-graph-complete-graph src/andres/graph/unit-test/graph-complete.cxx ${headers})
    add_test(test-graph-complete-graph test-graph-complete-graph)

    add_executable(test-graph-grid-graph src/andres/graph/unit-test/graph-grid.cxx ${headers})
    add_test(test-graph-grid-graph test-graph-grid-graph)

    add_executable(test-graph-dfs src/andres/graph/unit-test/dfs.cxx ${headers})
    add_test(test-graph-dfs test-graph-dfs)

    add_executable(test-graph-bfs src/andres/graph/unit-test/bfs.cxx ${headers})
    add_test(test-graph-bfs test-graph-bfs)

    add_executable(test-graph-cut-vertices src/andres/graph/unit-test/cut-vertices.cxx ${headers})
    add_test(test-graph-cut-vertices test-graph-cut-vertices)

    add_executable(test-graph-bridges src/andres/graph/unit-test/bridges.cxx ${headers})
    add_test(test-graph-bridges test-graph-bridges)

    add_executable(test-graph-minimum-spanning-tree src/andres/graph/unit-test/minimum-spanning-tree.cxx ${headers})
    add_test(test-graph-minimum-spanning-tree test-graph-minimum-spanning-tree)

    add_executable(test-graph-lifting src/andres/graph/unit-test/lifting.cxx ${headers})
    add_test(test-graph-lifting test-graph-lifting)

    add_executable(test-graph-multicut-lifted-kl src/andres/graph/unit-test/multicut-lifted/kernighan-lin.cxx ${headers})
    add_test(test-graph-multicut-lifted-kl test-graph-multicut-lifted-kl)

    add_executable(test-graph-multicut-lifted-greedy-additive src/andres/graph/unit-test/multicut-lifted/greedy-additive.cxx ${headers})
    add_test(test-graph-multicut-lifted-greedy-additive test-graph-multicut-lifted-greedy-additive)

    add_executable(test-graph-multicut-lifted-fusion-move-mc src/andres/graph/unit-test/multicut-lifted/fusion_move_mc.cxx ${headers})
    add_test(test-graph-multicut-lifted-fusion-move-mc test-graph-multicut-lifted-fusion-move-mc)


    add_executable(test-graph-pq src/andres/graph/unit-test/pq.cxx ${headers})
    add_test(test-graph-pq test-graph-pq)

    if(HDF5_FOUND)
        include_directories(PUBLIC ${HDF5_INCLUDE_DIRS})
        add_executable(test-hdf5 src/andres/graph/unit-test/hdf5.cxx ${headers})
        target_link_libraries(test-hdf5 ${HDF5_LIBRARIES})
        add_test(test-hdf5 test-hdf5)
    endif()
endif(COMPILER_SUPPORTS_CXX0X OR COMPILER_SUPPORTS_CXX11)

if(GUROBI_FOUND)
    add_executable(test-graph-multicut-ilp src/andres/graph/unit-test/multicut/ilp.cxx ${headers})
    target_link_libraries(test-graph-multicut-ilp ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CPP_LIBRARY} ${GUROBI_LIBRARY})
    add_test(test-graph-multicut-ilp test-graph-multicut-ilp)
endif()

if(DOXYGEN_FOUND)
    configure_file("${graph_SOURCE_DIR}/doxygen/doxyfile-graph.in" "${graph_BINARY_DIR}/doxyfile-graph" @ONLY IMMEDIATE)
    add_custom_target(doc-graph ALL COMMAND ${DOXYGEN} "${graph_BINARY_DIR}/doxyfile-graph")
endif()
