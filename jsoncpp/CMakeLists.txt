# vim: et ts=4 sts=4 sw=4 tw=0

# ==== Define cmake build policies that affect compilation and linkage default behaviors
#
# Set the JSONCPP_NEWEST_VALIDATED_POLICIES_VERSION string to the newest cmake version
# policies that provide successful builds. By setting JSONCPP_NEWEST_VALIDATED_POLICIES_VERSION
# to a value greater than the oldest policies, all policies between
# JSONCPP_OLDEST_VALIDATED_POLICIES_VERSION and CMAKE_VERSION (used for this build)
# are set to their NEW behaivor, thereby suppressing policy warnings related to policies
# between the JSONCPP_OLDEST_VALIDATED_POLICIES_VERSION and CMAKE_VERSION.
#
# CMake versions greater than the JSONCPP_NEWEST_VALIDATED_POLICIES_VERSION policies will
# continue to generate policy warnings "CMake Warning (dev)...Policy CMP0XXX is not set:"
#
set(JSONCPP_OLDEST_VALIDATED_POLICIES_VERSION "3.8.0")
set(JSONCPP_NEWEST_VALIDATED_POLICIES_VERSION "3.13.2")
cmake_minimum_required(VERSION ${JSONCPP_OLDEST_VALIDATED_POLICIES_VERSION})
if("${CMAKE_VERSION}" VERSION_LESS "${JSONCPP_NEWEST_VALIDATED_POLICIES_VERSION}")
    #Set and use the newest available cmake policies that are validated to work
    set(JSONCPP_CMAKE_POLICY_VERSION "${CMAKE_VERSION}")
else()
    set(JSONCPP_CMAKE_POLICY_VERSION "${JSONCPP_NEWEST_VALIDATED_POLICIES_VERSION}")
endif()
cmake_policy(VERSION ${JSONCPP_CMAKE_POLICY_VERSION})
#
# Now enumerate specific policies newer than JSONCPP_NEWEST_VALIDATED_POLICIES_VERSION
# that may need to be individually set to NEW/OLD
#
foreach(pnew "") # Currently Empty
    if(POLICY ${pnew})
        cmake_policy(SET ${pnew} NEW)
    endif()
endforeach()
foreach(pold "") # Currently Empty
    if(POLICY ${pold})
        cmake_policy(SET ${pold} OLD)
    endif()
endforeach()

# ==== Define language standard configurations requiring at least c++11 standard
if(CMAKE_CXX_STANDARD EQUAL "98" )
    message(FATAL_ERROR "CMAKE_CXX_STANDARD:STRING=98 is not supported.")
endif()

#####
##  Set the default target properties
if(NOT CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 11) # Supported values are ``11``, ``14``, and ``17``.
endif()
if(NOT CMAKE_CXX_STANDARD_REQUIRED)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()
if(NOT CMAKE_CXX_EXTENSIONS)
    set(CMAKE_CXX_EXTENSIONS OFF)
endif()

# ====

# Ensures that CMAKE_BUILD_TYPE has a default value
if(NOT DEFINED CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel Coverage.")
endif()

project(JSONCPP
        # Note: version must be updated in three places when doing a release. This
        # annoying process ensures that amalgamate, CMake, and meson all report the
        # correct version.
        # 1. /meson.build
        # 2. /include/json/version.h
        # 3. /CMakeLists.txt
        # IMPORTANT: also update the SOVERSION!!
        VERSION 1.9.3 # <major>[.<minor>[.<patch>[.<tweak>]]]
        LANGUAGES CXX)

message(STATUS "JsonCpp Version: ${JSONCPP_VERSION_MAJOR}.${JSONCPP_VERSION_MINOR}.${JSONCPP_VERSION_PATCH}")
set( JSONCPP_SOVERSION 22 )

option(JSONCPP_WITH_TESTS "Compile and (for jsoncpp_check) run JsonCpp test executables" ON)
option(JSONCPP_WITH_POST_BUILD_UNITTEST "Automatically run unit-tests as a post build step" ON)
option(JSONCPP_WITH_WARNING_AS_ERROR "Force compilation to fail if a warning occurs" OFF)
option(JSONCPP_WITH_STRICT_ISO "Issue all the warnings demanded by strict ISO C and ISO C++" ON)
option(JSONCPP_WITH_PKGCONFIG_SUPPORT "Generate and install .pc files" ON)
option(JSONCPP_WITH_CMAKE_PACKAGE "Generate and install cmake package files" ON)
option(JSONCPP_WITH_EXAMPLE "Compile JsonCpp example" OFF)
option(BUILD_SHARED_LIBS "Build jsoncpp_lib as a shared library." OFF)

# Enable runtime search path support for dynamic libraries on OSX
if(APPLE)
    set(CMAKE_MACOSX_RPATH 1)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mmacosx-version-min=10.9")
endif()

# Adhere to GNU filesystem layout conventions
include(GNUInstallDirs)

set(DEBUG_LIBNAME_SUFFIX "" CACHE STRING "Optional suffix to append to the library name for a debug build")

set(JSONCPP_USE_SECURE_MEMORY "0" CACHE STRING "-D...=1 to use memory-wiping allocator for STL" )

configure_file( "${PROJECT_SOURCE_DIR}/version.in"
                "${PROJECT_BINARY_DIR}/version"
                NEWLINE_STYLE UNIX )

macro(UseCompilationWarningAsError)
    if(MSVC)
        # Only enabled in debug because some old versions of VS STL generate
        # warnings when compiled in release configuration.
        if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.12.0)
          add_compile_options($<$<CONFIG:Debug>:/WX>)
        else()
          set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /WX ")
        endif()
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.12.0)
          add_compile_options(-Werror)
        else()
          set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
        endif()
        if(JSONCPP_WITH_STRICT_ISO)
          if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.12.0)
            add_compile_options(-pedantic-errors)
          else()
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic-errors")
          endif()
        endif()
    endif()
endmacro()

# Include our configuration header
include_directories( ${jsoncpp_SOURCE_DIR}/include )

if(MSVC)
    # Only enabled in debug because some old versions of VS STL generate
    # unreachable code warning when compiled in release configuration.
    if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.12.0)
      add_compile_options($<$<CONFIG:Debug>:/W4>)
    else()
      set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /W4 ")
    endif()
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    # using regular Clang or AppleClang
    if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.12.0)
      add_compile_options(-Wall -Wconversion -Wshadow -Werror=conversion -Werror=sign-compare)
    else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wconversion -Wshadow -Werror=conversion -Werror=sign-compare")
    endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # using GCC
    if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.12.0)
      add_compile_options(-Wall -Wconversion -Wshadow -Wextra)
    else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wconversion -Wshadow -Wextra")
    endif()
    # not yet ready for -Wsign-conversion

    if(JSONCPP_WITH_STRICT_ISO)
      if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.12.0)
        add_compile_options(-pedantic)
      else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic")
      endif()
    endif()
    if(JSONCPP_WITH_WARNING_AS_ERROR)
      if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.12.0)
        add_compile_options(-Werror=conversion)
      else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror=conversion")
      endif()
    endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    # using Intel compiler
    if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.12.0)
      add_compile_options(-Wall -Wconversion -Wshadow -Wextra -Werror=conversion)
    else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wconversion -Wshadow -Wextra -Werror=conversion")
    endif()

    if(JSONCPP_WITH_STRICT_ISO AND NOT JSONCPP_WITH_WARNING_AS_ERROR)
      if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.12.0)
        add_compile_options(-pedantic)
      else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic")
      endif()
    endif()
endif()

find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif(CCACHE_FOUND)

if(JSONCPP_WITH_WARNING_AS_ERROR)
    UseCompilationWarningAsError()
endif()

if(JSONCPP_WITH_PKGCONFIG_SUPPORT)
    configure_file(
        "pkg-config/jsoncpp.pc.in"
        "pkg-config/jsoncpp.pc"
        @ONLY)
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/pkg-config/jsoncpp.pc"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")
endif()

if(JSONCPP_WITH_CMAKE_PACKAGE)
        include (CMakePackageConfigHelpers)
        install(EXPORT jsoncpp
                DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/jsoncpp
                FILE        jsoncppConfig.cmake)
        write_basic_package_version_file ("${CMAKE_CURRENT_BINARY_DIR}/jsoncppConfigVersion.cmake"
                VERSION ${PROJECT_VERSION}
                COMPATIBILITY SameMajorVersion)
        install(FILES ${CMAKE_CURRENT_BINARY_DIR}/jsoncppConfigVersion.cmake
                DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/jsoncpp)
endif()

if(JSONCPP_WITH_TESTS)
  enable_testing()
  include(CTest)
endif()

# Build the different applications
add_subdirectory( src )

#install the includes
add_subdirectory( include )

#install the example
if(JSONCPP_WITH_EXAMPLE)
  add_subdirectory( example )
endif()
