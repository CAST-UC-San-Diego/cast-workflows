"""
associaTR class for running GWAS
"""

import numpy as np
import os
import pandas as pd
import tempfile
import trtools.associaTR

class AssociaTRRunner:
    def __init__(self, ptcovar, trvcf, region=None, covars=[]):
        self.ptcovar = ptcovar
        self.pt_npy = None
        self.covar_npy = None
        self.trvcf = trvcf
        self.region = region
        self.covars = covars
        self.gwas = None
        self.method = "associaTR"
        self.setup()

    def setup(self):
        # Set up phenotypes and covariates
        self.pt_npy = tempfile.NamedTemporaryFile(suffix='.npy')
        self.covar_npy = tempfile.NamedTemporaryFile(suffix='.npy')
        np.save(self.pt_npy.name, \
            self.ptcovar[["person_id", "phenotype"]].to_numpy(dtype=float))
        np.save(self.covar_npy.name, \
            self.ptcovar[["person_id"]+self.covars].to_numpy(dtype=float))

        # Fetch genotypes file if remote
        if self.trvcf.startswith("gs://"):
            trvcf = self.trvcf.split("/")[-1]
            if not os.path.isfile(trvcf):
                os.system("gsutil -u ${GOOGLE_PROJECT} cp %s ."%self.trvcf)
            self.trvcf = trvcf

    def RunGWAS(self):
        # Set associaTR options
        vcftype = None
        same_samples = False
        sample_list = None
        non_major_cutoff = 20
        beagle_dosages = True
        plotting_phenotype = None
        paired_genotype_plot = False
        plot_genotype_residuals = False
        plotting_ci_alphas = []
        imputed_ukb_strs_paper_period_check = False
        outfile = tempfile.NamedTemporaryFile(suffix='.tsv')
        trtools.associaTR.perform_gwas(
            outfile.name,
            self.trvcf,
            "phenotype",
            [self.pt_npy, self.covar_npy],
            vcftype,
            same_samples,
            sample_list,
            self.region,
            non_major_cutoff,
            beagle_dosages,
            plotting_phenotype,
            paired_genotype_plot,
            plot_genotype_residuals,
            plotting_ci_alphas,
            imputed_ukb_strs_paper_period_check
        )
        gwas = pd.read_csv(outfile.name, sep="\t")
        gwas.rename({"p_phenotype": "p_value", "coeff_phenotype": "beta", \
        	"se_phenotype": "standard_error"}, inplace=True, axis=1)
        gwas["-log10pvalue"] = -np.log10(gwas["p_value"])
        self.gwas = gwas
