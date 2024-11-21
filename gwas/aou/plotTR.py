#!/usr/bin/env python3

"""
Make a plot summarizing TR/trait association and allele freqs

Example:
./plotTR.py --phenotype platelet_count --tr-vcf CBL_test.filtered.sorted.vcf.gz --region chr11:119206290-119206323 --outprefix CBL_platelet_count
"""

# Allow making plots even with no x-forward
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import argparse
import numpy as np
import os
import pandas as pd
import scipy.stats
import trtools.utils.tr_harmonizer as trh
import trtools.utils.utils as utils

from scipy.stats import sem
import scipy.stats as st 

from utils import MSG, ERROR

SAMPLEFILE = os.path.join(os.environ["WORKSPACE_BUCKET"], "samples", \
    "passing_samples_v7.csv")

def GetPTCovarPath(phenotype):
    return os.path.join(os.getenv('WORKSPACE_BUCKET'), \
        "phenotypes", "%s_phenocovar.csv"%phenotype)


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--phenotype", help="Phenotypes file path, or phenotype name", type=str, required=True)
    parser.add_argument("--samples", help="List of sample IDs, sex to keep", type=str, default=SAMPLEFILE)
    parser.add_argument("--tr-vcf", help="VCF file with TR genotypes.", type=str, required=True)
    parser.add_argument("--region", help="chr:start-end of target TR", type=str, required=True)
    parser.add_argument("--outprefix", help="Output prefix", type=str, required=True)
    parser.add_argument("--min-samples-per-dosage", help="Filter dosages with fewer samples", type=int, default=100)
    args = parser.parse_args()

    # Load phenotype data
    if args.phenotype.endswith(".csv"):
        ptcovar_path = args.phenotype
    else:
        ptcovar_path = GetPTCovarPath(args.phenotype)
    data = pd.read_csv(ptcovar_path)
    data["person_id"] = data["person_id"].apply(str)

    # Load TR genotypes
    invcf = utils.LoadSingleReader(args.tr_vcf, checkgz=True)
    samples = invcf.samples
    region = invcf(args.region)
    nrecords = 0
    for record in region:
        trrecord = trh.HarmonizeRecord(trh.VcfTypes["hipstr"], record)
        afreqs = trrecord.GetAlleleFreqs()
        genotypes = trrecord.GetLengthGenotypes()
        dosages = [sum(item)/2 for item in genotypes]
        trdf = pd.DataFrame({"person_id": samples, "tr_dosage": dosages})
        nrecords += 1
    if nrecords == 0:
        ERROR("No matching TR records found")
    if nrecords > 1:
        ERROR("Multiple matching TR records found")

    # Merge phenotype and TR dosages
    df = pd.merge(data, trdf, on=["person_id"])
    df.to_csv('plot.csv',index=False,sep=',')
    pltdata = df.groupby("tr_dosage", as_index=False).agg(phenotype_mean=("phenotype", np.mean), n=("phenotype", len),sem=("phenotype", sem))
    pltdata = pltdata[pltdata["n"]>args.min_samples_per_dosage]
    pltdata.to_csv('CI.csv',index=False,sep=',')
    # Compute CI
    pltdata[['CI_lower','CI_upper']] = pltdata.apply(lambda x:st.norm.interval(alpha=0.95, loc=x['phenotype_mean'], scale=x['sem']),axis=1).tolist()
    pltdata = pltdata[pltdata["n"]>args.min_samples_per_dosage]
    #print(pltdata)
    #pltdata.to_csv('CI.csv',index=False,sep=',')





	# Plot - TODO
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(pltdata["tr_dosage"], pltdata["phenotype_mean"])
    fig.savefig(args.outprefix + ".png")
if __name__ == "__main__":
    main()