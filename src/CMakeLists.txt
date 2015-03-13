set ( LIBRARY_FORMAT_SRCS
    LibraryFormat.cpp
)

set ( MINCE_SRCS
    Mince.cpp
    FindPartition.cpp
    Decoder.cpp
    BucketModel.cpp
    bloom.c
    misc.c
)

set ( CANONICALIZE_SRCS
    Canonicalize.cpp
)

set ( CANONICAL_COMPARE_SRCS 
    CanonicalCompare.cpp
)

set ( MERGE_PAIRS_SRCS
    MergePairs.cpp
)

include_directories( 
${GAT_SOURCE_DIR}/include
${GAT_SOURCE_DIR}/external
${GAT_SOURCE_DIR}/external/cereal/include
${GAT_SOURCE_DIR}/external/install/include
${GAT_SOURCE_DIR}/external/install/include/jellyfish-2.1.3
${ZLIB_INCLUDE_DIR}
${TBB_INCLUDE_DIRS}
${Boost_INCLUDE_DIRS}
)

link_directories( 
${GAT_SOURCE_DIR}/lib
${GAT_SOURCE_DIR}/external/install/lib
${Boost_LIBRARY_DIRS}
${TBB_LIBRARY_DIRS}
)

message("TBB_LIBRARIES = ${TBB_LIBRARIES}")
message("Boost_LIBRARIES = ${Boost_LIBRARIES}")

# Set the RPATH 
if (APPLE)
    ## This DOES NOT do what I / any one sane, expects.  Setting the 
    ## linker path on OSX is messed up.  Just tell the user to use 
    ## DYLD_FALLBACK_LIBRARY_PATH for now
    set(CMAKE_INSTALL_RPATH "@executable_path/../lib:@executable_path/../../lib:@executable_path/")
else()
    set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib:$ORIGIN/../../lib:$ORIGIN/")
endif()

set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)

# Build the library 
add_library(minceutils STATIC MinceUtils.cpp LibraryFormat.cpp bitfile.cpp)
add_dependencies(minceutils externalProjectTarget)

# Build the mince executable
add_executable(mince ${MINCE_SRCS})
add_dependencies(mince minceutils)

# Build the canonicalize executable
add_executable(canonicalize ${CANONICALIZE_SRCS})

# Build the mergepairs executable
add_executable(mergepairs ${MERGE_PAIRS_SRCS})

# Build the canoncmp executable
add_executable(canoncmp ${CANONICAL_COMPARE_SRCS})

# Link the executable
target_link_libraries(mince
    minceutils
    ${Boost_LIBRARIES} 
    ${ZLIB_LIBRARY} 
    ${GAT_SOURCE_DIR}/external/install/lib/libjellyfish-2.0.a 
    pthread 
    m
    tbb
    tbbmalloc
    ${NON_APPLECLANG_LIBS}
)

# Link the executable
target_link_libraries(canonicalize
    minceutils
    ${Boost_LIBRARIES} 
    ${ZLIB_LIBRARY} 
    ${GAT_SOURCE_DIR}/external/install/lib/libjellyfish-2.0.a 
    pthread 
    m
    tbb
    tbbmalloc
    ${NON_APPLECLANG_LIBS}
)

# Link the executable
target_link_libraries(canoncmp
    minceutils
    ${Boost_LIBRARIES} 
    ${ZLIB_LIBRARY} 
    ${GAT_SOURCE_DIR}/external/install/lib/libjellyfish-2.0.a 
    pthread 
    m
    tbb 
    tbbmalloc
    ${NON_APPLECLANG_LIBS}
)


# Link the executable
target_link_libraries(mergepairs
    minceutils
    ${Boost_LIBRARIES} 
    ${ZLIB_LIBRARY} 
    ${GAT_SOURCE_DIR}/external/install/lib/libjellyfish-2.0.a 
    pthread 
    m
    tbb
    tbbmalloc
    ${NON_APPLECLANG_LIBS}
)

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(CMAKE_INSTALL_PREFIX
    "${GAT_SOURCE_DIR}" CACHE PATH "Default install prefix" FORCE
    )
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

set(INSTALL_LIB_DIR lib )
set(INSTALL_BIN_DIR bin )
set(INSTALL_INCLUDE_DIR include )

install(DIRECTORY 
        ${GAT_SOURCE_DIR}/external/install/lib/ 
        DESTINATION ${INSTALL_LIB_DIR}
    )

install(TARGETS mince 
                RUNTIME DESTINATION bin 
                LIBRARY DESTINATION lib
                ARCHIVE DESTINATION lib
        )

    #set(POST_INSTALL_SCRIPT ${GAT_SOURCE_DIR}/cmake/PostInstall.cmake)

install(
    CODE 
    "
    execute_process(COMMAND \"${CMAKE_COMMAND}\"
                            -DCMAKE_SYSTEM_NAME=${CMAKE_SYSTEM_NAME}
                            -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} 
                            -P \"${POST_INSTALL_SCRIPT}\")
    "
)

include(InstallRequiredSystemLibraries)
add_test( NAME simple_test COMMAND ${CMAKE_COMMAND} -DTOPLEVEL_DIR=${GAT_SOURCE_DIR} -P ${GAT_SOURCE_DIR}/cmake/SimpleTest.cmake )


