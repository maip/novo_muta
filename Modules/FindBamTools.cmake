# Copyright 2012 Tobias Marschall
#
# This file is part of CLEVER.
#
# CLEVER is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CLEVER is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with CLEVER. If not, see <http://www.gnu.org/licenses/>.


# - Try to find BamTools (https://github.com/pezmaster31/bamtools)
# Once done this will define
# BamTools_FOUND
# BamTools_INCLUDE_DIRS
# BamTools_LIBRARIES
# BamTools_DEFINITIONS

set(BamTools_PREFIX "/usr/local/")

find_path(BamTools_INCLUDE_DIR bamtools/api/api_global.h HINTS ${BamTools_PREFIX}/include )
find_library(BamTools_LIBRARY NAMES bamtools HINTS ${BamTools_PREFIX}/lib/bamtools )

set(BamTools_LIBRARIES ${BamTools_LIBRARY})
set(BamTools_INCLUDE_DIRS ${BamTools_INCLUDE_DIR} ${BamTools_INCLUDE_DIR}/bamtools )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set BamTools_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(BamTools "BamTools library (https://github.com/pezmaster31/bamtools) not found. If it is in a non-standard place, you have to set the variable BamTools_PREFIX, for example by adding -DBamTools_PREFIX=<path> to your cmake call" BamTools_LIBRARY BamTools_INCLUDE_DIR)

IF (BAMTOOLS_FOUND)
set(BamTools_FOUND TRUE)
ENDIF()

mark_as_advanced(BamTools_INCLUDE_DIR BamTools_LIBRARY)