# 
#  Copyright (C) 2024 Gilles Villard
#
#   This file is part of mapml. mapml is free software: you
#   can redistribute it and/or modify it under the terms of the GNU Lesser
#   General Public License as published by the Free Software Foundation,
#   either version 2.1 of the License, or (at your option) any later version.
#
#   mapml is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU Lesser General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public License
#   along with mapml. If not, see <http://www.gnu.org/licenses/>. 


if WITH_MAPLE

AUTOMAKE_OPTIONS = gnu

AM_CPPFLAGS = $(MAPLE_CPPFLAGS) 


AM_LDFLAGS = $(MAPLE_LDFLAGS)

LIBS = $(MAPLE_LIBS) 

include_pmldir= $(includedir)/pml

include_pml_HEADERS =  mapml.h conversion.h matpoly_export.h

lib_LTLIBRARIES = libmapml.la

libmapml_la_SOURCES = mapml.h conversion.h conversion.c dummy_export.c \
				     matpoly_export.h  matpoly_export.c 

endif
