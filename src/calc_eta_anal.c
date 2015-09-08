#include "svmanalys.h"

#define NUM_FILES 2
#define NUM_MU 20
#define NUM_SIG 14

void calc_eta_anal(int natoms, const char *fname);

static output_env_t oenv = NULL;

int main(int argc, char *argv[]) {
	const char *desc[] = {
		"calc_eta_anal calculates analytical eta values.",
		"It produces two output files, for analytical mu changed etas and analytical sigma changed etas."
	};
	const char *mu_fname, *sig_fname;
	int num_mu = NUM_MU, num_sig = NUM_SIG;

	t_filenm fnm[] = {
		{efDAT, "-fmu", "eta_anal_mu.dat", ffWRITE},
		{efDAT, "-fsig", "eta_anal_sig.dat", ffWRITE}
	};

	t_pargs pa[] = {
		{"-nmu", FALSE, etINT, {&num_mu}, "number of mu changes to calculate"},
		{"-nsig", FALSE, etINT, {&num_sig}, "number of sigma changes to calculate"}
	};

	parse_common_args(&argc, argv, 0, NUM_FILES, fnm, 
		asize(pa), pa, asize(desc), desc, 0, NULL, &oenv);

	mu_fname = opt2fn("-fmu", NUM_FILES, fnm);
	sig_fname = opt2fn("-fsig", NUM_FILES, fnm);

	calc_eta_anal(num_mu, num_sig, mu_fname, sig_fname);

}

void calc_eta_anal(int num_mu, int num_sig, const char *mu_fname, const char *sig_fname) {
	FILE *f_mu, *f_sig;
	int i;

	f_mu = fopen(mu_fname, "w");
	f_sig = fopen(sig_fname, "w");

	for(i = 2; i <= num_mu + 1; i++) {
		fprintf(f_mu, "%f\n", 
			erf(((i - 1.0) * 0.5) / sqrt(2)));
	}

	for(i = 2; i <= num_sig + 1; i++) {
		fprintf(f_sig, "%f\n", 
			exp((-2.0 * log(i) / (i * i - 1))) - exp((-2.0 * i * i * log(i)) / (i * i - 1.0)));
	}

	fclose(f_mu);
	fclose(f_sig);
}