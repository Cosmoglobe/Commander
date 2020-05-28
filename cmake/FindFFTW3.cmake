
#--------------------------------------------------------------------------------
# Copyright (c) 2012-2013, Lars Baehren <lbaehren@gmail.com>
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
#
#  * Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#--------------------------------------------------------------------------------

# - Check for the presence of FFTW3
#
# The following variables are set when FFTW3 is found:
#  FFTW3_FOUND      = Set to true, if all components of FFTW3 have been found.
#  FFTW3_INCLUDES   = Include path for the header files of FFTW3
#  FFTW3_LIBRARIES  = Link these to use FFTW3
#  FFTW3_LFLAGS     = Linker flags (optional)

if (NOT FFTW3_FOUND)

  if (NOT FFTW3_ROOT_DIR)
    set (FFTW3_ROOT_DIR ${CMAKE_INSTALL_PREFIX})
  endif (NOT FFTW3_ROOT_DIR)

  ##_____________________________________________________________________________
  ## Check for the header files

  find_path (FFTW3_INCLUDES fftw3.h
    HINTS ${FFTW3_ROOT_DIR} ${CMAKE_INSTALL_PREFIX}
    PATH_SUFFIXES include include/fftw3
    )

  ##_____________________________________________________________________________
  ## Check for the library

  find_library (FFTW3_LIBRARIES fftw3
    HINTS ${FFTW3_ROOT_DIR} ${CMAKE_INSTALL_PREFIX}
    PATH_SUFFIXES lib
    )

  ##_____________________________________________________________________________
  ## Actions taken when all components have been found

  find_package_handle_standard_args (FFTW3 DEFAULT_MSG FFTW3_LIBRARIES FFTW3_INCLUDES)

  if (FFTW3_FOUND)
    if (NOT FFTW3_FIND_QUIETLY)
      message (STATUS "Found components for FFTW3")
      message (STATUS "FFTW3_ROOT_DIR  = ${FFTW3_ROOT_DIR}")
      message (STATUS "FFTW3_INCLUDES  = ${FFTW3_INCLUDES}")
      message (STATUS "FFTW3_LIBRARIES = ${FFTW3_LIBRARIES}")
    endif (NOT FFTW3_FIND_QUIETLY)
  else (FFTW3_FOUND)
    if (FFTW3_FIND_REQUIRED)
      message (FATAL_ERROR "Could not find FFTW3!")
    endif (FFTW3_FIND_REQUIRED)
  endif (FFTW3_FOUND)

  ##_____________________________________________________________________________
  ## Mark advanced variables

  mark_as_advanced (
    FFTW3_ROOT_DIR
    FFTW3_INCLUDES
    FFTW3_LIBRARIES
    )

endif (NOT FFTW3_FOUND)
