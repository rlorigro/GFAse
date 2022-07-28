######################################
# ---------------------------------- #
# -------- COMPILER VERSION -------- #
# ---------------------------------- #
######################################

cmake_minimum_required(VERSION 3.11)
project(GFAse VERSION 0.0.0)

message(STATUS "CMAKE_CXX_COMPILER_ID: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMAKE_CXX_COMPILER_VERSION: ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_SYSTEM: ${CMAKE_SYSTEM}")

# Compiler options.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (dev)
    message(STATUS "---- Building with debug options ----")

    # Memory debug
    set(CMAKE_CXX_FLAGS "-fexceptions -fsanitize=address")
    set(ASAN_OPTIONS=check_initialization_order=1)
    set(ASAN_OPTIONS=detect_leaks=1)

    add_definitions(-ggdb3 -O0 -Wall)       # Debugging + No optimization

else()
    message(STATUS "---- Building with optimization ----")

    # Standard compilation
    add_definitions(-O3 -Wall)              # Much optimization
endif()



#########################################
# ------------------------------------- #
# -------- SOURCES AND HEADERS -------- #
# ------------------------------------- #
#########################################

# Include header files
include_directories(
        "inc"
)

# Define our shared library sources. NOT test/executables.
set(SOURCES
        src/Bam.cpp
        src/BinaryIO.cpp
        src/BinarySequence.cpp
        src/Bipartition.cpp
        src/BubbleGraph.cpp
        src/ContactGraph.cpp
        src/Color.cpp
        src/FixedBinarySequence.cpp
        src/GfaReader.cpp
        src/gfa_to_handle.cpp
        src/graph_utility.cpp
        src/HaplotypePathKmer.cpp
        src/Hasher.cpp
        src/Hasher2.cpp
        src/handle_to_gfa.cpp
        src/IncrementalIdMap.cpp
        src/KmerSets.cpp
        src/misc.cpp
        src/MurmurHash2.cpp
        src/MurmurHash3.cpp
##        src/OverlapMap.cpp
        src/Phase.cpp
        src/PhaseAssign.cpp
        src/Sequence.cpp
        src/Sam.cpp
        src/Timer.cpp
)


project(GFAse)
add_library(GFAse STATIC ${SOURCES})

# To make sure the library is named Name.so,
# get rid of the "lib" prefix.
set_target_properties(GFAse PROPERTIES PREFIX "")

# Eliminate an extraneous -D during compilation.
set_target_properties(GFAse PROPERTIES DEFINE_SYMBOL "")


##############################################
# ------------------------------------------ #
# --------- LINKING NATIVE LIBRARY --------- #
# ------------------ OMP ------------------- #
# ------------------------------------------ #
##############################################


if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

    # assumes clang build
    # we can't reliably detect when we're using clang, so for the time being we assume
    # TODO: can't we though?

    # adapted from https://stackoverflow.com/questions/46414660/macos-cmake-and-openmp
    # find_package(OpenMP) does not work reliably on macOS, so we do its work ourselves
    if(EXISTS /usr/local/Cellar/libomp/)
        set(OMPDIR /usr/local/Cellar/libomp/)
    elseif(EXISTS /opt/local/include/libomp)
        set(OMPDIR /opt/local/include/libomp)
    endif()

    set (OpenMP_C "${CMAKE_C_COMPILER}")
    set (OpenMP_C_FLAGS " -Xpreprocessor -fopenmp -I${OMPDIR} -I/usr/local/include -L${OMPDIR} -L/usr/local/lib")
    set (OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
    set (OpenMP_CXX "${CMAKE_CXX_COMPILER}")
    set (OpenMP_CXX_FLAGS " -Xpreprocessor -fopenmp -I${OMPDIR} -I/usr/local/include -L${OMPDIR} -L/usr/local/lib")
    set (OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
    set (OpenMP_libomp_LIBRARY "omp")
    set (OpenMP_libgomp_LIBRARY "gomp")
    set (OpenMP_libiomp5_LIBRARY "iomp5")

    # and now add the OpenMP parameters to the compile flags
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS} -lomp")

    # Mac needs libdl and libomp when linking the library
    set(PLATFORM_EXTRA_LIB_FLAGS -ldl -lomp)


elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    find_package(OpenMP REQUIRED)

    # To link statically, it is necessary to call lower level commands, as opposed to just
    # target_link_library(target OpenMP::OpenMP_CXX), or else make tries to -lgomp which is apparently
    # bad because it brings in dynamic libraries dlsym, dlopen, dlwtf
    if (OPENMP_FOUND)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()

endif()


##############################################
# ------------------------------------------ #
# -------- LINKING EXTERNAL LIBRARY -------- #
# ------------------ bdsg ------------------ #
# ------------------------------------------ #
##############################################

# Need to explicitly enable ExternalProject functionality
include(ExternalProject)

message(STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")
message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}")
message(STATUS "CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")

if (dev)
    # Download or update library as an external project
    ExternalProject_Add(project_bdsg
            GIT_REPOSITORY https://github.com/vgteam/libbdsg.git
            DOWNLOAD_COMMAND ""
            UPDATE_COMMAND ""
            PREFIX ${CMAKE_SOURCE_DIR}/external/bdsg/
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external/bdsg/ -DRUN_DOXYGEN=OFF -DBUILD_PYTHON_BINDINGS=OFF
            BUILD_IN_SOURCE True
            INSTALL_DIR ${CMAKE_SOURCE_DIR}/external/bdsg/
            INSTALL_COMMAND make install
            )
else()
    # Download or update library as an external project
    ExternalProject_Add(project_bdsg
            GIT_REPOSITORY https://github.com/vgteam/libbdsg.git
            PREFIX ${CMAKE_SOURCE_DIR}/external/bdsg/
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external/bdsg/ -DRUN_DOXYGEN=OFF -DBUILD_PYTHON_BINDINGS=OFF
            BUILD_IN_SOURCE True
            INSTALL_DIR ${CMAKE_SOURCE_DIR}/external/bdsg/
            INSTALL_COMMAND make install
            )
endif()

# Define INSTALL_DIR as the install directory for external library
ExternalProject_Get_Property(project_bdsg INSTALL_DIR)

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(LIB_SUFFIX "dylib")
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(LIB_SUFFIX "so")
endif()


# Create new library for external project (so it can be linked with main library)
add_library(bdsg SHARED IMPORTED)
set_property(TARGET bdsg
        PROPERTY IMPORTED_LOCATION ${INSTALL_DIR}/lib/libbdsg.${LIB_SUFFIX})

add_library(divsufsort SHARED IMPORTED)
set_property(TARGET divsufsort
        PROPERTY IMPORTED_LOCATION ${INSTALL_DIR}/lib/libdivsufsort64.${LIB_SUFFIX})

add_library(libhandlegraph SHARED IMPORTED)
set_property(TARGET libhandlegraph
        PROPERTY IMPORTED_LOCATION ${INSTALL_DIR}/lib/libhandlegraph.${LIB_SUFFIX})

add_library(libsdsl SHARED IMPORTED)
set_property(TARGET libsdsl
        PROPERTY IMPORTED_LOCATION ${INSTALL_DIR}/lib/libsdsl.${LIB_SUFFIX})

# Define library as dependent on the downloaded project
add_dependencies(bdsg
        project_bdsg
        libsdsl
        libhandlegraph
        divsufsort)

# Define main library as dependent on the downloaded project (transitively)
add_dependencies(GFAse bdsg)

# Ensure that main library has access to primary dependencies' and secondary dependencies' headers
include_directories(external/bdsg/include/
        external/bdsg/include/handlegraph/
        external/bdsg/include/bdsg/)

# Add runtime path for main library so it can access external library
#set_property(TARGET GFAse PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/bin/bdsg")

message(STATUS "INSTALL_DIR: ${INSTALL_DIR}")


###############################################
## ------------------------------------------ #
## -------- LINKING EXTERNAL LIBRARY -------- #
## --------------- bamtools ----------------- #
## ------------------------------------------ #
###############################################
#
## Need to explicitly enable ExternalProject functionality
#include(ExternalProject)
#
#message(STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")
#message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
#message(STATUS "CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}")
#message(STATUS "CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")
#
#if (dev)
#    # Download or update library as an external project
#    ExternalProject_Add(project_bamtools
#            GIT_REPOSITORY https://github.com/pezmaster31/bamtools.git
#            DOWNLOAD_COMMAND ""
#            UPDATE_COMMAND ""
#            PREFIX ${CMAKE_SOURCE_DIR}/external/bamtools/
#            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external/bamtools/
#            BUILD_IN_SOURCE True
#            INSTALL_DIR ${CMAKE_SOURCE_DIR}/external/bamtools/
#            INSTALL_COMMAND make install
#            )
#else()
#    # Download or update library as an external project
#    ExternalProject_Add(project_bamtools
#            GIT_REPOSITORY https://github.com/pezmaster31/bamtools.git
#            PREFIX ${CMAKE_SOURCE_DIR}/external/bamtools/
#            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external/bamtools/
#            BUILD_IN_SOURCE True
#            INSTALL_DIR ${CMAKE_SOURCE_DIR}/external/bamtools/
#            INSTALL_COMMAND make install
#            )
#endif()
#
## Define INSTALL_DIR as the install directory for external library
#ExternalProject_Get_Property(project_bamtools INSTALL_DIR)
#
#if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
#    error("NO MACs")
#elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
#    set(LIB_SUFFIX "a")
#endif()
#
#
## Create new library for external project (so it can be linked with main library)
#add_library(bamtools STATIC IMPORTED)
#set_property(TARGET bamtools
#        PROPERTY IMPORTED_LOCATION ${INSTALL_DIR}/lib/libbamtools.${LIB_SUFFIX})
#
## Define library as dependent on the downloaded project
#add_dependencies(bamtools
#        project_bamtools
#        )
#
## Define main library as dependent on the downloaded project (transitively)
#add_dependencies(GFAse bamtools)
#
## Ensure that main library has access to primary dependencies' and secondary dependencies' headers
#include_directories(external/bamtools/include/
#        external/bamtools/include/bamtools/
#        external/bamtools/include/bamtools/api
#        external/bamtools/include/bamtools/algorithms
#        external/bamtools/include/bamtools/shared
#        )
#
## Add runtime path for main library so it can access external library
##set_property(TARGET GFAse PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/bin/bamtools")
#
#message(STATUS "INSTALL_DIR: ${INSTALL_DIR}")


##############################################
# ------------------------------------------ #
# -------- LINKING EXTERNAL LIBRARY -------- #
# ---------------- htslib ------------------ #
# ------------------------------------------ #
##############################################

## Need to explicitly enable ExternalProject functionality
#include(ExternalProject)
#
#message(STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")
#message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
#message(STATUS "CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}")
#message(STATUS "CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")
#
## Dependencies
#find_package(ZLIB REQUIRED)
#
#find_library(BZ2LIB REQUIRED)
#
#find_library(CURLLIB REQUIRED)
#
##if (dev)
#    # Download or update library as an external project
#    ExternalProject_Add(project_htslib
#            GIT_REPOSITORY https://github.com/samtools/htslib.git
##            DOWNLOAD_COMMAND ""
##            UPDATE_COMMAND ""
#            PREFIX ${CMAKE_SOURCE_DIR}/external/htslib/
#            CMAKE_ARGS ""
#            CONFIGURE_COMMAND autoreconf -i COMMAND ./configure --disable-lzma --disable-libcurl --disable-bz2 --disable-s3 --disable-gcs --without-libdeflate --disable-plugins
#            BUILD_COMMAND $(MAKE)
#            INSTALL_COMMAND $(MAKE) prefix=${CMAKE_SOURCE_DIR}/external/htslib/ install
#            BUILD_IN_SOURCE True
#            INSTALL_DIR ${CMAKE_SOURCE_DIR}/external/htslib/
#            )
##else()
##    # Download or update library as an external project
##    ExternalProject_Add(project_htslib
##            GIT_REPOSITORY https://github.com/samtools/htslib.git
##            PREFIX ${CMAKE_SOURCE_DIR}/external/htslib/
##            CMAKE_ARGS ""
##            BUILD_COMMAND "autoreconf -i"
##            INSTALL_COMMAND "make install"
##            BUILD_IN_SOURCE True
##            INSTALL_DIR ${CMAKE_SOURCE_DIR}/external/htslib/
##            INSTALL_COMMAND make install
##            )
##endif()
#
## Define INSTALL_DIR as the install directory for external library
#ExternalProject_Get_Property(project_htslib INSTALL_DIR)
#
#if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
#    error("NO MACs")
#elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
#    set(LIB_SUFFIX "a")
#endif()
#
#
## Create new library for external project (so it can be linked with main library)
#add_library(hts STATIC IMPORTED)
#set_property(TARGET hts
#        PROPERTY IMPORTED_LOCATION ${INSTALL_DIR}/lib/libhts.${LIB_SUFFIX})
#
## Define library as dependent on the downloaded project
#add_dependencies(hts
#        project_htslib
#        )
#
## Define main library as dependent on the downloaded project (transitively)
#add_dependencies(GFAse hts)
#
## Ensure that main library has access to primary dependencies' and secondary dependencies' headers
#include_directories(
#        external/
#        external/htslib/
#        external/htslib/include/
#        external/htslib/include/htslib
#        )

# ----------------------------


#set(HTS_SOURCES
#        ${INSTALL_DIR}/bcf_sr_sort.c
#        ${INSTALL_DIR}/bgzf.c
#        ${INSTALL_DIR}/bgzip.c
#        ${INSTALL_DIR}/errmod.c
#        ${INSTALL_DIR}/faidx.c
#        ${INSTALL_DIR}/header.c
#        ${INSTALL_DIR}/hfile.c
#        ${INSTALL_DIR}/hfile_gcs.c
#        ${INSTALL_DIR}/hfile_libcurl.c
#        ${INSTALL_DIR}/hfile_s3.c
#        ${INSTALL_DIR}/hfile_s3_write.c
#        ${INSTALL_DIR}/hts.c
#        ${INSTALL_DIR}/hts_expr.c
#        ${INSTALL_DIR}/htsfile.c
#        ${INSTALL_DIR}/hts_os.c
#        ${INSTALL_DIR}/kfunc.c
#        ${INSTALL_DIR}/kstring.c
#        ${INSTALL_DIR}/md5.c
#        ${INSTALL_DIR}/multipart.c
#        ${INSTALL_DIR}/plugin.c
#        ${INSTALL_DIR}/probaln.c
#        ${INSTALL_DIR}/realn.c
#        ${INSTALL_DIR}/regidx.c
#        ${INSTALL_DIR}/region.c
#        ${INSTALL_DIR}/sam.c
#        ${INSTALL_DIR}/synced_bcf_reader.c
#        ${INSTALL_DIR}/tabix.c
#        ${INSTALL_DIR}/tbx.c
#        ${INSTALL_DIR}/textutils.c
#        ${INSTALL_DIR}/thread_pool.c
#        ${INSTALL_DIR}/vcf.c
#        ${INSTALL_DIR}/vcf_sweep.c
#        ${INSTALL_DIR}/vcfutils.c
#        )
#
#add_library(htslib ${HTSLIB_SOURCES})
#add_dependencies(htslib HTSLIB_CONFIGURED)
#target_link_libraries(htslib pthread bz2 z curl)
#
#message(STATUS "INSTALL_DIR: ${INSTALL_DIR}")


# --------------------

message(STATUS "Fetching htslib")


# https://github.com/samtools/htslib/releases/download/1.15.1/htslib-1.15.1.tar.bz2

find_package(ZLIB REQUIRED)

# Dependencies
find_library(ZLIB NAMES z)
if(${ZLIB} STREQUAL "ZLIB-NOTFOUND")
    message(WARNING "Couldn't find the 'z' library")
endif()

find_library(BZ2LIB bz2)
if(${BZ2LIB} STREQUAL "BZ2LIB-NOTFOUND")
    message(WARNING "Couldn't find the 'bz2' library")
endif()

find_library(CURLLIB curl)
if(${CURLLIB} STREQUAL "CURLLIB-NOTFOUND")
    message(WARNING "Couldn't find the 'curl' library")
endif()

#set(HTSLIB_SOURCES
#        ${project_htslib_SOURCE_DIR}/cram/cram_codecs.c
#        ${project_htslib_SOURCE_DIR}/cram/cram_decode.c
#        ${project_htslib_SOURCE_DIR}/cram/cram_encode.c
#        ${project_htslib_SOURCE_DIR}/cram/cram_external.c
#        ${project_htslib_SOURCE_DIR}/cram/cram_index.c
#        ${project_htslib_SOURCE_DIR}/cram/cram_io.c
#        ${project_htslib_SOURCE_DIR}/cram/cram_stats.c
#        ${project_htslib_SOURCE_DIR}/cram/mFILE.c
#        ${project_htslib_SOURCE_DIR}/cram/open_trace_file.c
#        ${project_htslib_SOURCE_DIR}/cram/pooled_alloc.c
#        ${project_htslib_SOURCE_DIR}/cram/string_alloc.c
#        ${project_htslib_SOURCE_DIR}/bcf_sr_sort.c
#        ${project_htslib_SOURCE_DIR}/bgzf.c
#        ${project_htslib_SOURCE_DIR}/bgzip.c
#        ${project_htslib_SOURCE_DIR}/errmod.c
#        ${project_htslib_SOURCE_DIR}/faidx.c
#        ${project_htslib_SOURCE_DIR}/header.c
#        ${project_htslib_SOURCE_DIR}/hfile.c
#        ${project_htslib_SOURCE_DIR}/hfile_gcs.c
#        ${project_htslib_SOURCE_DIR}/hfile_libcurl.c
#        ${project_htslib_SOURCE_DIR}/hfile_s3.c
#        ${project_htslib_SOURCE_DIR}/hfile_s3_write.c
#        ${project_htslib_SOURCE_DIR}/hts.c
#        ${project_htslib_SOURCE_DIR}/hts_expr.c
#        ${project_htslib_SOURCE_DIR}/htsfile.c
#        ${project_htslib_SOURCE_DIR}/hts_os.c
#        ${project_htslib_SOURCE_DIR}/kfunc.c
#        ${project_htslib_SOURCE_DIR}/kstring.c
#        ${project_htslib_SOURCE_DIR}/md5.c
#        ${project_htslib_SOURCE_DIR}/multipart.c
#        ${project_htslib_SOURCE_DIR}/plugin.c
#        ${project_htslib_SOURCE_DIR}/probaln.c
#        ${project_htslib_SOURCE_DIR}/realn.c
#        ${project_htslib_SOURCE_DIR}/regidx.c
#        ${project_htslib_SOURCE_DIR}/region.c
#        ${project_htslib_SOURCE_DIR}/sam.c
#        ${project_htslib_SOURCE_DIR}/synced_bcf_reader.c
#        ${project_htslib_SOURCE_DIR}/tabix.c
#        ${project_htslib_SOURCE_DIR}/tbx.c
#        ${project_htslib_SOURCE_DIR}/textutils.c
#        ${project_htslib_SOURCE_DIR}/thread_pool.c
#        ${project_htslib_SOURCE_DIR}/vcf.c
#        ${project_htslib_SOURCE_DIR}/vcf_sweep.c
#        ${project_htslib_SOURCE_DIR}/vcfutils.c
#        )

include(FetchContent)

FetchContent_Declare(
        project_htslib
        URL      https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
#        PREFIX ${CMAKE_SOURCE_DIR}/external/htslib/
#        CMAKE_ARGS ""
#        CONFIGURE_COMMAND autoreconf -i COMMAND ./configure --disable-lzma --disable-libcurl --disable-bz2 --disable-s3 --disable-gcs --without-libdeflate --disable-plugins
#        INSTALL_COMMAND $(MAKE) prefix=${CMAKE_SOURCE_DIR}/external/htslib/ install
#        INSTALL_DIR ${CMAKE_SOURCE_DIR}/external/htslib/
#        BINARY_DIR ${CMAKE_SOURCE_DIR}/external/htslib/
)

FetchContent_MakeAvailable(project_htslib)
#FetchContent_populate(project_htslib)

set(HTSLIB_SOURCES
        ${project_htslib_SOURCE_DIR}/cram/cram_codecs.c
        ${project_htslib_SOURCE_DIR}/cram/cram_decode.c
        ${project_htslib_SOURCE_DIR}/cram/cram_encode.c
        ${project_htslib_SOURCE_DIR}/cram/cram_external.c
        ${project_htslib_SOURCE_DIR}/cram/cram_index.c
        ${project_htslib_SOURCE_DIR}/cram/cram_io.c
        ${project_htslib_SOURCE_DIR}/cram/cram_samtools.c
        ${project_htslib_SOURCE_DIR}/cram/cram_stats.c
        ${project_htslib_SOURCE_DIR}/cram/files.c
        ${project_htslib_SOURCE_DIR}/cram/mFILE.c
        ${project_htslib_SOURCE_DIR}/cram/open_trace_file.c
        ${project_htslib_SOURCE_DIR}/cram/pooled_alloc.c
        ${project_htslib_SOURCE_DIR}/cram/rANS_static.c
        ${project_htslib_SOURCE_DIR}/cram/sam_header.c
        ${project_htslib_SOURCE_DIR}/cram/string_alloc.c
        ${project_htslib_SOURCE_DIR}/os/rand.c
#        ${project_htslib_SOURCE_DIR}/test/fieldarith.c
#        ${project_htslib_SOURCE_DIR}/test/hfile.c
#        ${project_htslib_SOURCE_DIR}/test/hts_endian.c
#        ${project_htslib_SOURCE_DIR}/test/sam.c
#        ${project_htslib_SOURCE_DIR}/test/test-bcf-sr.c
#        ${project_htslib_SOURCE_DIR}/test/test-regidx.c
#        ${project_htslib_SOURCE_DIR}/test/test-vcf-api.c
#        ${project_htslib_SOURCE_DIR}/test/test-vcf-sweep.c
#        ${project_htslib_SOURCE_DIR}/test/test_bgzf.c
#        ${project_htslib_SOURCE_DIR}/test/test_view.c
#        ${project_htslib_SOURCE_DIR}/test/thrash_threads1.c
#        ${project_htslib_SOURCE_DIR}/test/thrash_threads2.c
#        ${project_htslib_SOURCE_DIR}/test/thrash_threads3.c
#        ${project_htslib_SOURCE_DIR}/test/thrash_threads4.c
#        ${project_htslib_SOURCE_DIR}/test/thrash_threads5.c
#        ${project_htslib_SOURCE_DIR}/test/thrash_threads6.c
        ${project_htslib_SOURCE_DIR}/bcf_sr_sort.c
        ${project_htslib_SOURCE_DIR}/bgzf.c
        ${project_htslib_SOURCE_DIR}/bgzip.c
        ${project_htslib_SOURCE_DIR}/errmod.c
        ${project_htslib_SOURCE_DIR}/faidx.c
        ${project_htslib_SOURCE_DIR}/hfile.c
        ${project_htslib_SOURCE_DIR}/hfile_gcs.c
        ${project_htslib_SOURCE_DIR}/hfile_libcurl.c
        ${project_htslib_SOURCE_DIR}/hfile_net.c
#        ${project_htslib_SOURCE_DIR}/hfile_s3.c
        ${project_htslib_SOURCE_DIR}/hts.c
        ${project_htslib_SOURCE_DIR}/htsfile.c
        ${project_htslib_SOURCE_DIR}/kfunc.c
        ${project_htslib_SOURCE_DIR}/knetfile.c
        ${project_htslib_SOURCE_DIR}/kstring.c
        ${project_htslib_SOURCE_DIR}/md5.c
        ${project_htslib_SOURCE_DIR}/multipart.c
#        ${project_htslib_SOURCE_DIR}/plugin.c
        ${project_htslib_SOURCE_DIR}/probaln.c
        ${project_htslib_SOURCE_DIR}/realn.c
        ${project_htslib_SOURCE_DIR}/regidx.c
        ${project_htslib_SOURCE_DIR}/sam.c
        ${project_htslib_SOURCE_DIR}/synced_bcf_reader.c
        ${project_htslib_SOURCE_DIR}/tabix.c
        ${project_htslib_SOURCE_DIR}/tbx.c
        ${project_htslib_SOURCE_DIR}/textutils.c
        ${project_htslib_SOURCE_DIR}/thread_pool.c
        ${project_htslib_SOURCE_DIR}/vcf.c
        ${project_htslib_SOURCE_DIR}/vcf_sweep.c
        ${project_htslib_SOURCE_DIR}/vcfutils.c
        )

message(STATUS "htslib_SOURCE_DIR: ${project_htslib_SOURCE_DIR}")

## HTSLIB configure
#add_custom_command(
#        OUTPUT "${htslib_SOURCE_DIR}/external/htslib/config.h"
#        COMMAND autoreconf -i
#        COMMAND ./configure --disable-lzma --disable-libcurl --disable-bz2 --disable-s3 --disable-gcs --without-libdeflate --disable-plugins
#        COMMAND make prefix=${CMAKE_SOURCE_DIR}/external/htslib/ install
#        WORKING_DIRECTORY ${htslib_SOURCE_DIR}
#)


set(HTS_DIR _deps/project_htslib-src)

# HTSLIB configure
add_custom_target(
        BUILD_HTS
        ALL
#        OUTPUT ${PROJECT_SOURCE_DIR}/external/htslib/config.h # //TODO: at build time ${htslib_SOURCE_DIR} is not actually populated!! is an empty string
        WORKING_DIRECTORY ${HTS_DIR}
#        COMMENT "AAAAAAA ${HTS_DIR}"
        COMMAND pwd
        COMMAND autoconf
        COMMAND autoheader
        COMMAND ./configure --disable-lzma --disable-libcurl --disable-bz2 --disable-s3 --disable-gcs --without-libdeflate --disable-plugins
        COMMAND $(MAKE) prefix=${CMAKE_SOURCE_DIR}/external/htslib/ install
)

# Ensure that main library has access to primary dependencies' and secondary dependencies' headers
include_directories(
        ${CMAKE_BINARY_DIR}/_deps/project_htslib-src
#        ${CMAKE_BINARY_DIR}/_deps/project_htslib-src/htslib/
#        ${CMAKE_BINARY_DIR}/_deps/project_htslib-src/cram/
#        ${CMAKE_BINARY_DIR}/_deps/project_htslib-src/os/
        ${CMAKE_SOURCE_DIR}/external/
)

#add_custom_target(HTSLIB_CONFIGURED DEPENDS "${CMAKE_BINARY_DIR}/_deps/project_htslib-src/config.h")

add_library(htslib ${HTSLIB_SOURCES})
add_dependencies(htslib BUILD_HTS)
add_dependencies(GFAse htslib)

target_link_libraries(htslib pthread bz2 z curl)


##############################################
# ------------------------------------------ #
# -------- LINKING EXTERNAL LIBRARY -------- #
# ---------------- sparsepp ---------------- #
# ------------------------------------------ #
##############################################


ExternalProject_Add(project_sparsepp
        GIT_REPOSITORY https://github.com/greg7mdp/sparsepp.git
        TIMEOUT 5
        PREFIX "${CMAKE_CURRENT_BINARY_DIR}"
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND "" # Disable install step, is a header only lib!
        )

add_dependencies(GFAse project_sparsepp)

# Specify include dir
ExternalProject_Get_Property(project_sparsepp SOURCE_DIR)
include_directories(${SOURCE_DIR}/sparsepp)


##############################################
# ------------------------------------------ #
# -------- LINKING EXTERNAL LIBRARY -------- #
# ----------------- snarls ----------------- #
# ------------------------------------------ #
##############################################

#
## Need to explicitly enable ExternalProject functionality
#include(ExternalProject)
#
#message(STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")
#message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
#message(STATUS "CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}")
#message(STATUS "CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")
#
#if (dev)
#    # Download or update library as an external project
#    ExternalProject_Add(project_snarls
#            GIT_REPOSITORY https://github.com/vgteam/libsnarls.git
#            GIT_TAG main
#            DOWNLOAD_COMMAND ""
#            UPDATE_COMMAND ""
#            PREFIX ${CMAKE_SOURCE_DIR}/external/snarls/
#            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external/snarls/ -DRUN_DOXYGEN=OFF -DBUILD_PYTHON_BINDINGS=OFF
#            BUILD_IN_SOURCE True
#            INSTALL_DIR ${CMAKE_SOURCE_DIR}/external/snarls/
#            INSTALL_COMMAND make install
#            )
#else()
#    # Download or update library as an external project
#    ExternalProject_Add(project_snarls
#            GIT_REPOSITORY https://github.com/vgteam/libsnarls.git
#            GIT_TAG main
#            PREFIX ${CMAKE_SOURCE_DIR}/external/snarls/
#            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external/snarls/ -DRUN_DOXYGEN=OFF -DBUILD_PYTHON_BINDINGS=OFF
#            BUILD_IN_SOURCE True
#            INSTALL_DIR ${CMAKE_SOURCE_DIR}/external/snarls/
#            INSTALL_COMMAND make install
#            )
#endif()
#
## Define INSTALL_DIR as the install directory for external library
#ExternalProject_Get_Property(project_snarls INSTALL_DIR)
#
## Create new library for external project (so it can be linked with main library)
#add_library(snarls SHARED IMPORTED)
#set_property(TARGET snarls
#        PROPERTY IMPORTED_LOCATION ${INSTALL_DIR}/lib/libsnarls.${LIB_SUFFIX})
#
##add_library(libsdsl SHARED IMPORTED)
##set_property(TARGET libsdsl
##        PROPERTY IMPORTED_LOCATION ${INSTALL_DIR}/lib/libsdsl.${LIB_SUFFIX})
#
## Define library as dependent on the downloaded project
#add_dependencies(
#        snarls
#        project_snarls
#)
#
## Define main library as dependent on the downloaded project (transitively)
#add_dependencies(GFAse snarls)
#
## Ensure that main library has access to primary dependencies' and secondary dependencies' headers
#include_directories(
#        external/snarls/include/
#        external/snarls/include/libvgio/
#        external/snarls/include/structures/
#        external/snarls/include/dynamic/
#        external/snarls/include/gtest/
#        external/snarls/include/mio/
#        external/snarls/include/tsl/
#        external/snarls/include/vg/
#        )
#
#message(STATUS "INSTALL_DIR: ${INSTALL_DIR}")
#
#
############################################
# ---------------------------------------- #
# -------- Generating executables -------- #
# ---------------------------------------- #
############################################

#
##set(Boost_DEBUG 1)
#find_package( Boost 1.60 COMPONENTS program_options REQUIRED )
#include_directories( ${Boost_INCLUDE_DIR} )
#
#

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# -------- TESTS --------

set(TESTS
        test_assign_phase
        test_bfs
        test_bubble_graph
        test_binary_sequence
        test_binary_sequence_performance
        test_connected_component_finder
        test_contact_graph
        test_fixed_binary_sequence
        test_fixed_binary_sequence_performance_2
        test_fixed_binary_sequence_sparsepp_performance
        test_gfareader
        test_htslib
        test_htslib_bam_reader
        test_kmer_unordered_set
        test_phase_haplotype_paths
        test_nonbinary_sequence_performance
        test_nonbinary_sequence_sparsepp_performance
        test_rgb_to_hex
        test_set_intersection
        test_timer
)

foreach(FILENAME_PREFIX ${TESTS})
    add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
    target_link_libraries(${FILENAME_PREFIX}
            GFAse
            Threads::Threads
#            bamtools
            ZLIB::ZLIB
            bdsg
#            ${Boost_LIBRARIES}
            htslib
            divsufsort
            libhandlegraph
            libsdsl)

endforeach()


# -------- EXECUTABLES --------

set(EXECUTABLES
        assign_phases_via_diploid_alignment
        count_kmers
        create_bandage_path_color_table
        compute_minhash
        compute_minhash2
        extract_haplotype_kmers
        evaluate_phasing
        find_bubbles
        get_mapq_distribution_from_sam
        get_degree_stats
        get_path_lengths
        generate_contact_map_from_bam
        locate_kmer_matches
        phase_haplotype_paths
        phase_contacts
        phase_contacts_as_tripartition
        remove_empty_nodes
        split_connected_components
        unzip
)

foreach(FILENAME_PREFIX ${EXECUTABLES})
    add_executable(${FILENAME_PREFIX} src/executable/${FILENAME_PREFIX}.cpp)
    target_link_libraries(${FILENAME_PREFIX}
            GFAse
            Threads::Threads
#            bamtools
            ZLIB::ZLIB
            bdsg
#            ${Boost_LIBRARIES}
            divsufsort
            htslib
            libhandlegraph
            libsdsl)

endforeach()

#set_target_properties(gfase PROPERTIES LINK_FLAGS "-static" )
#SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -static-libstdc++ -static-libgcc")


# -------- final steps --------

# Where to install
set(BINARY_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)

## The shared library goes to the bin directory.
#install(TARGETS
#        GFAse
#        gfase
#        DESTINATION ${BINARY_INSTALL_DIR})

# Some potentially helpful messages
message(STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")
message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}")
message(STATUS "CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")
