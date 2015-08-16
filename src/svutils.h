/*
 *
 *                This source code is KINDA part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.5.3
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2010, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon

 * Written by Ahnaf Siddiqui and Dr. Sameer Varma.
 * Copyright (c) 2015, University of South Florida.
 */

#ifndef _svutils_h
#define _svutils_h

#include <stdarg.h>
#include <stdio.h>
#include "gmx_fatal.h"
#include "macros.h"
#include "smalloc.h"
#include "statutil.h"
#include "typedefs.h"

/* Enumeration of input and output files */
enum {TRAJ1, TRAJ2, NDX1, NDX2, COORD_PDB, COORD_DAT, RES_DAT, MAXFILES};

/*
 * Gets file names from command line
 * int files[] are the desired files based on files enum defined above
 * Requested file names will be stored in fnames
 */
void get_file_args(int argc, char *argv[], const char *desc[], int desc_size, int files[], const char *fnames[], int num_files);

/* Opens the logfile and logs initial time/date */
void init_log(char *program);

/* Closes the logfile */
void close_log(void);

/* Prints to both stdout and the logfile */
void print_log(char const *fmt, ...);

/* 
 * Logs fatal error to logfile and also calls gmx_fatal
 * Hint: Use FARGS for the first 3 arguments.
 */
void log_fatal(int fatal_errno, const char *file, int line, char const *fmt, ...);

#endif