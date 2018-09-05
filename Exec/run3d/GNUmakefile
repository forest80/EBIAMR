#AMREX_HOME defines the directory in which we will find the BoxLib directory
AMREX_HOME ?= ../../../amrex

#TOP defines the directory in which we will find Source, Exec, etc.
TOP = ../..

#
# Variables for the user to set ...
#

DIM        = 3
COMP	   = gnu
DEBUG	   = FALSE
USE_MPI    = FALSE
USE_OPENMP = FALSE

PRECISION = DOUBLE

USE_FLCTS = TRUE
USE_FLCTS = FALSE

USE_VELOCITY = TRUE
USE_VELOCITY = FALSE

USE_F90_SOLVERS = FALSE
USE_F90_SOLVERS = TRUE

USE_XBLAS = TRUE
USE_XBLAS = FALSE

EBASE     = amr

ifeq (${USE_FLCTS}, TRUE)
  DEFINES += -DBL_DO_FLCT
endif

# This is for the time-dependent turbulence forcing term
DEFINES += -DDO_IAMR_FORCE

Bpack   := ./Make.package
Blocs   := .

include $(TOP)/Exec/Make.IAMR
