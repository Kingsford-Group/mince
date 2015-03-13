cmake_minimum_required (VERSION 2.8)

enable_testing()

project (Mince)

set(CPACK_PACKAGE_VERSION "0.6.0")
set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "6")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_GENERATOR "TGZ")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Mince - Practical and efficient compression boosting for sequencing reads")
set(CPACK_PACKAGE_NAME 
  "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME
  "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}-Source")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")


set (WARNING_IGNORE_FLAGS "-Wno-deprecated-register") 
set (BOOST_CXX_FLAGS "-Wno-deprecated-register -std=c++11")
## Prefer static to dynamic libraries
SET(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})

## Set the standard required compile flags
set (CMAKE_CXX_FLAGS "-funroll-loops -fPIC -O3 -DHAVE_ANSI_TERM -DHAVE_SSTREAM -DHAVE_CONFIG_H -Wall -std=c++11 -Wno-unused-local-typedefs")

##
# OSX is strange (some might say, stupid in this regard).  Deal with it's quirkines here.
##
if (APPLE)
    # To allow ourselves to build a dynamic library, we have to tell the compiler
    # that, yes, the symbols will be around at runtime.
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -undefined dynamic_lookup")
    set (LIBMINCE_LINKER_FLAGS "-all_load")
    # In order to "think different", we also have to use non-standard suffixes
    # for our shared libraries
    set(SHARED_LIB_EXTENSION "dylib")
else()
    # We're in sane linux world
   set (SHARED_LIB_EXTENSION "so")
   set (LIBMINCE_LINKER_FLAGS "")
endif()

set( BOOST_EXTRA_FLAGS "--layout=tagged" )
## this get's set differently below if we 
## are on clang & apple
set (NON_APPLECLANG_LIBS gomp rt)

##
# Compiler-specific C++11 activation.
# http://stackoverflow.com/questions/10984442/how-to-detect-c11-support-of-a-compiler-with-cmake
##
##
# First take care of what to do if we have gcc
##
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    execute_process(
        COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    # If we're on OSX
    if (APPLE AND NOT (GCC_VERSION VERSION_GREATER 4.8.2 OR GCC_VERSION VERSION_EQUAL 4.8.2))
        message(FATAL_ERROR "When building under OSX, ${PROJECT_NAME} requires "
                            "either clang or g++ >= 4.8.2")
    elseif (NOT (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7))
        message(FATAL_ERROR "${PROJECT_NAME} requires g++ 4.7 or greater.")
    endif ()

    set (GCC TRUE)
    #Perhaps the logger now works in both worlds
    set (CMAKE_CXX_FLAGS "-DHAVE_LOGGER ${CMAKE_CXX_FLAGS}")
    set (HAVE_LOGGER TRUE)

    set(WARNING_IGNORE_FLAGS "${WARNING_IGNORE_FLAGS} -Wno-unused-local-typedefs")
    set (BOOST_CONFIGURE_TOOLSET "")
    set (BOOST_EXTRA_FLAGS "")
    # Tentatively, we support clang now
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(CLANG TRUE)
    # If we're on OSX using clang, try to use libc++ instead of libstdc++
    if (APPLE)
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
        set (BOOST_CONFIGURE_TOOLSET "--with-toolset=clang")
        set (BOOST_EXTRA_FLAGS toolset=clang cxxflags="-stdlib=libc++"  linkflags="-stdlib=libc++")
    endif()
    # There's currently a bug with clang-3.4 & Boost 1.55 -- this hack fixes it
    # but we should do something better (does this break things if CPU doesn't
    # have 128-bit support)?
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_HAS_INT128") 
    set(CMAKE_SHARK_CXX_FLAGS "-DBOOST_HAS_INT128") 
    # We should have the logger under clang
    set (CMAKE_CXX_FLAGS "-DHAVE_LOGGER ${CMAKE_CXX_FLAGS}")
    set (HAVE_LOGGER TRUE)
    if (APPLE)
        set (NON_APPLECLANG_LIBS "")
    endif()
else ()
    message(FATAL_ERROR "Your C++ compiler does not support C++11.")
endif ()

##
# Record this top-level path
##
set (GAT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# Have CMake tell us what it's doing
set (CMAKE_VERBOSE_MAKEFILE true)

##
# Super-secret override
##
if ( DEFINED CUSTOM_BOOST_PATH )
	set (CMAKE_INCLUDE_PATH ${CUSTOM_BOOST_PATH} ${CMAKE_INCLUDE_PATH})
    set (CMAKE_LIBRARY_PATH ${CUSTOM_BOOST_PATH}/lib ${CMAKE_LIBRARY_PATH})
endif ( DEFINED CUSTOM_BOOST_PATH )

##
# We want dyanmic, multithreaded boost libraries
##
set (Boost_USE_STATIC_LIBS OFF)
set (Boost_USE_MULTITHREADED ON)
set (Boost_USE_STATIC_RUNTIME OFF)

##
# Set the latest version and look for what we need
##
set(Boost_ADDITIONAL_VERSIONS "1.54" "1.54.0")
find_package(Boost 1.54.0 COMPONENTS iostreams filesystem system thread timer program_options serialization)

## TODO: Figure out how to detect this automatically
# If the "assembler" is too old, tell TBB not to compile
# with -mrtm
if (NO_RTM)
    set (TBB_CXXFLAGS "-mno-rtm")
endif()

# Recurse into Sailfish source directory
add_subdirectory ( src )

include(ExternalProject)

##
# Either inform the user of how to obtain Boost, or, if they passed in the FETCH_BOOST
# option, go and grab it for them.
##
if ((NOT Boost_FOUND) AND (NOT FETCH_BOOST))
	message(FATAL_ERROR 
		"Sailfish cannot be compiled without Boost.\n"
        "It is recommended to visit http://www.boost.org/ and install Boost according to those instructions.\n"
        "This build system can also download and install a local version of boost for you (this takes a lot of time).\n"
        "To fetch and build boost locally, call cmake with -DFETCH_BOOST=TRUE"
    )
elseif(FETCH_BOOST)
    ## Let the rest of the build process know we're going to be fetching boost
    set(BOOST_WILL_RECONFIGURE TRUE)
    message("Build system will fetch and build Boost")
    message("==================================================================")
    ExternalProject_Add(libboost
        DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
        URL http://en.sourceforge.jp/frs/g_redir.php?m=jaist&f=%2Fboost%2Fboost%2F1.55.0%2Fboost_1_55_0.tar.gz
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/boost_1_55_0
        INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
        CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/external/boost_1_55_0/bootstrap.sh ${BOOST_CONFIGURE_TOOLSET} --prefix=<INSTALL_DIR>
        BUILD_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/external/boost_1_55_0/b2 ${BOOST_EXTRA_FLAGS} install
        BUILD_IN_SOURCE 1
        INSTALL_COMMAND ""
    )

    ##
    # After we've installed boost, 
    ##
    ExternalProject_Add_Step(libboost reconfigure
        COMMAND ${CMAKE_COMMAND} ${CMAKE_CURRENT_SOURCE_DIR} -DBOOST_WILL_RECONFUGRE=FALSE -DBOOST_RECONFIGURE=TRUE 
        DEPENDEES install
    )
    add_dependencies(minceutils libboost)
endif()

##
# If we're fetching boost and we need to have dummy paths for these variables
# so that CMake won't complain
##
if (BOOST_WILL_RECONFIGURE)
    message("Setting Temporary Boost paths")
    unset(Boost_FOUND CACHE)
    unset(Boost_INCLUDE_DIR CACHE)
    unset(Boost_INCLUDE_DIRS CACHE)
    unset(Boost_LIBRARY_DIRS CACHE) 
    unset(Boost_LIBRARIES CACHE)
    unset(BOOST_ROOT CACHE)
    set(Boost_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install/include)
    set(Boost_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/external/install/include)
    set(Boost_LIBRARY_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/external/install/lib)        
    set(Boost_FOUND TRUE)
endif()

##
#  If we had to fetch Boost, the reconfigure step will re-run cmake.  The second configuration
#  pass is executed with the BOOST_RECONFIGURE flag set. This should allow our newly
#  installed Boost to be found by CMake.
##
if (BOOST_RECONFIGURE)
    message("In Boost Reconfiguration")
    unset(Boost_FOUND CACHE)
    unset(Boost_INCLUDE_DIR CACHE)
    unset(Boost_INCLUDE_DIRS CACHE)
    unset(Boost_LIBRARY_DIRS CACHE) 
    unset(Boost_LIBRARIES CACHE)
    unset(BOOST_ROOT CACHE)
    unset(CMAKE_PREFIX_PATH CACHE)
    set(CMAKE_PREFIX_PATH ${CMAKE_CURRENT_SOURCE_DIR}/external/install)
    set(Boost_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install/include)
    set(Boost_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/external/install/include)
    set(Boost_LIBRARY_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/external/install/lib)
    find_package(Boost 1.55.0 COMPONENTS iostreams filesystem system thread timer program_options)
endif()


message("BOOST INCLUDE DIR = ${Boost_INCLUDE_DIR}")
message("BOOST INCLUDE DIRS = ${Boost_INCLUDE_DIRS}")
message("BOOST LIB DIR = ${Boost_LIBRARY_DIRS}")
message("BOOST LIBRAREIS = ${Boost_LIBRARIES}")
find_package (ZLIB)

set(EXTERNAL_LIBRARY_PATH $CMAKE_CURRENT_SOURCE_DIR/lib)

message("Build system will fetch and build the Cereal serialization library")
message("==================================================================")
include(ExternalProject)
ExternalProject_Add(libcereal
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    URL http://www.cs.cmu.edu/~robp/files/cereal-v1.0.0.tgz
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/cereal-1.0.0
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    UPDATE_COMMAND sh -c "mkdir -p <SOURCE_DIR>/build"
    BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/cereal-1.0.0/build
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND sh -c "mkdir -p <INSTALL_DIR>/include && cp -r <SOURCE_DIR>/include/cereal <INSTALL_DIR>/include"
)
add_dependencies(minceutils libcereal)

message("Build system will fetch and build Jellyfish")
message("==================================================================")
ExternalProject_Add(libjellyfish
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    URL ftp://ftp.genome.umd.edu/pub/jellyfish/jellyfish-2.1.3.tar.gz
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/jellyfish-2.1.3
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/external/jellyfish-2.1.3/configure --prefix=<INSTALL_DIR>
    BUILD_COMMAND ${MAKE}
    BUILD_IN_SOURCE 1
    INSTALL_COMMAND make install && 
                    cp config.h <INSTALL_DIR>/include/jellyfish-2.1.3/jellyfish/ &&
                    cp config.h <INSTALL_DIR>/include/
)
add_dependencies(minceutils libjellyfish)

find_package(TBB)
##
#
# Fetch and build Intel's Threading Building Blocks library.
#
##
if(NOT TBB_FOUND)
    set(TBB_WILL_RECONFIGURE TRUE)
    # Set the appropriate compiler
    if (CLANG)
        set(TBB_COMPILER "clang")
    else()
        set(TBB_COMPILER "gcc")
    endif()
    message("Build system will fetch and build Intel Threading Building Blocks")
    message("==================================================================")
    # These are useful for the custom install step we'll do later
    set(TBB_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/tbb43_20140724oss)
    set(TBB_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install)
    if ("${TBB_COMPILER}" STREQUAL "gcc")
        ## Don't know why it's a problem yet, but if we're using
        ## GCC, get rid of the DO_ITT_NOTIFY flag
        set(TBB_CXXFLAGS "${TBB_CXXFLAGS} -UDO_ITT_NOTIFY")
    endif()
    ExternalProject_Add(libtbb
        DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
        URL http://www.threadingbuildingblocks.org/sites/default/files/software_releases/source/tbb43_20140724oss_src.tgz
        DOWNLOAD_COMMAND curl -k -L http://www.threadingbuildingblocks.org/sites/default/files/software_releases/source/tbb43_20140724oss_src.tgz -o tbb_20140724oss_src.tgz &&
        tar -xzvf tbb_20140724oss_src.tgz
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/tbb43_20140724oss
        INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
        PATCH_COMMAND "${TBB_PATCH_STEP}"
        CONFIGURE_COMMAND ""
        BUILD_COMMAND make CXXFLAGS=${TBB_CXXFLAGS} lambdas=1 compiler=${TBB_COMPILER} cfg=release tbb_build_prefix=LIBS
        INSTALL_COMMAND sh -c "cp ${TBB_SOURCE_DIR}/build/LIBS_release/*.${SHARED_LIB_EXTENSION}* ${TBB_INSTALL_DIR}/lib && cp -r ${TBB_SOURCE_DIR}/include/* ${TBB_INSTALL_DIR}/include"
        BUILD_IN_SOURCE 1
        )
    SET( RECONFIG_FLAGS ${RECONFIG_FLAGS} -DTBB_WILL_RECONFIGURE=FALSE -DTBB_RECONFIGURE=TRUE)
    ExternalProject_Add_Step(libtbb reconfigure
        COMMAND ${CMAKE_COMMAND} ${CMAKE_CURRENT_SOURCE_DIR} ${RECONFIG_FLAGS}
        DEPENDEES install
        )
    add_dependencies(minceutils libtbb)
endif()
##
# If we're fetching tbb, we need to have dummy paths for these variables
# so that CMake won't complain
##
if(TBB_WILL_RECONFIGURE)
    set(TBB_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/external/install/include)
    set(TBB_LIBRARY_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/external/install/lib)
    set(TBB_LIBRARIES tbb tbbmalloc)
    #set(TBB_LIBRARIES ${CMAKE_CURRENT_SOURCE_DIR}/external/install/lib/libtbb.${SHARED_LIB_EXTENSION}
    #    ${CMAKE_CURRENT_SOURCE_DIR}/external/install/lib/libtbbmalloc.${SHARED_LIB_EXTENSION}
    #    )
endif()
##
# Similar to the Boost trick above, the libtbb reconfigure should force this code
# to be run on the second configuration pass, where it should appropriately set the
# TBB_INSTALL_DIR variable.
##
if (TBB_RECONFIGURE)
    unset(TBB_FOUND CACHE)
    unset(TBB_INSTALL_DIR CACHE)
    unset(CMAKE_PREFIX_PATH CACHE)
    set(CMAKE_PREFIX_PATH ${CMAKE_CURRENT_SOURCE_DIR}/external/install)
    set(TBB_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install)
    message("TBB_INSTALL_DIR = ${TBB_INSTALL_DIR}")
    find_package(TBB)
endif()
message("TBB_LIBRARIES = ${TBB_LIBRARIES}")

###
#
# Done building external dependencies.
#
###

set (CPACK_SOURCE_IGNORE_FILES 
"/build/"
"/bin/"
"/lib/"
"/external/"
"/src/obsolete/"
"/include/obsolete/"
".git/")

message("CPACK_SOURCE_IGNORE_FILES = ${CPACK_SOURCE_IGNORE_FILES}")

# build a CPack driven installer package
include (CPack)

