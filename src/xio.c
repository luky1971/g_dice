/*
 * Written by Ahnaf Siddiqui.
 *
 * Link with -lgmx.$(CPU) when compiling.
 */

#include <stdio.h>
#include "/usr/local/gromacs/include/gromacs/xtcio.h"

int main(int argc, char *argv[]) {
	t_fileio *fio = NULL;

	fio = open_xtc(argv[1], "r");

	printf("%s\n", fio->fn);
}