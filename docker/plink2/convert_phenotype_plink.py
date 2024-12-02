#!/usr/bin/env python3



"""
Convert AoU phenotypes to plink friendly format.
Output cleaned phenotype and comined covariates files
for the specified phenotype

Usage:
./convert_phenotype_plink.py --phenotype <phenotype>

"""

import os
import pandas as pd
import subprocess
import sys
import csv
import argparse


# Get token and set up project
token_fetch_command = subprocess.run(['gcloud', 'auth', 'application-default', 'print-access-token'], \
    capture_output=True, check=True, encoding='utf-8')
token = str.strip(token_fetch_command.stdout)
project = os.getenv("GCS_REQUESTER_PAYS_PROJECT")
bucket = os.getenv("WORKSPACE_BUCKET")
#print(f"Workspace Bucket: {bucket}")
#print all the environmental variables 
#print(dict(os.environ))

def GetPTCovarPath(phenotype):
    return os.path.join(bucket, \
        "phenotypes", "%s_phenocovar.csv"%phenotype)

def DownloadPT(filename):
    """
	Download phenotype_covar.csv locally

	Arguments
	---------
	filename : str
	   GCP path
    """
    cmd = "gsutil cp {filename} .".format(filename=filename)
    output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read()
    print(output.decode("utf-8"))

    return filename.split("/")[-1]  # Return local filename


def GetFloatFromPC(x):
    x = x.replace("[","").replace("]","")

    return float(x)

def LoadAncestry(ancestry_pred_path,project):
    """
	Download ancestry_pred.tsv locally

	Arguments
	---------
	ancestry_pred_path : str
	   GCP path
    project : str
    google project for downloading
    """
    if ancestry_pred_path.startswith("gs://"):
        os.system(f"gsutil -u {project} cp {ancestry_pred_path} .")
        ancestry_pred_path = "ancestry_preds.tsv"
    ancestry =  pd.read_csv(ancestry_pred_path, sep="\t")
    ancestry.rename({"research_id": "IID"}, axis=1, inplace=True)
    ancestry['IID'] = ancestry['IID'].astype(str)
    num_pcs = len(ancestry["pca_features"].values[0].split(","))
    pcols = ["PC_%s"%i for i in range(num_pcs)]
    ancestry[pcols] = ancestry["pca_features"].str.split(",", expand=True)
    for p in pcols:
        ancestry[p] = ancestry[p].apply(lambda x: GetFloatFromPC(x), 1)
    ancestry.insert(0, 'FID', 0)

    return ancestry

def convert_csv_to_plink (ptfile):
    df = pd.read_csv(ptfile)
    df.insert(0, 'FID', 0)
    df.rename(columns={"person_id": "IID"}, inplace=True)
    df['IID'] = df['IID'].astype(str)

    return df

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--phenotype", help="Comma-separated list of phenotypes to process. If not provided, process all available phenotypes", type=str, required=True)
    parser.add_argument("--num-pcs", help="Number of PCs to use as covariates", type=int, default=10)
    parser.add_argument("--ancestry-pred-path", help="Path to ancestry predictions",type=str, default="gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv")
    parser.add_argument("--ptcovars", help="Comma-separated list of phenotype-specific covariates. Default: age", type=str, default="age")
    parser.add_argument("--sharedcovars", help="Comma-separated list of shared covariates (besides PCs). Default: sex_at_birth_Male", type=str, default="sex_at_birth_Male")
    parser.add_argument("--dryrun", help="Don't actually run the workflow. Just set up", action="store_true")
    args = parser.parse_args()


    # Set up paths
    if args.phenotype.endswith(".csv"):
        ptcovar_path = args.phenotype
    else:
        ptcovar_path = GetPTCovarPath(args.phenotype)


    # Get covarlist
    pcols = ["PC_%s"%i for i in range(1, args.num_pcs+1)]
    shared_covars = [item for item in args.sharedcovars.split(",") if item != ""]
    pt_covars = [item for item in args.ptcovars.split(",") if item != ""]
    covars = pt_covars + shared_covars

    # Set up data frame with phenotype and covars
    ancestry = LoadAncestry(args.ancestry_pred_path,project)
    plink = convert_csv_to_plink(DownloadPT(ptcovar_path))

    # Extract phenotype and covars only
    data = pd.merge(plink[["FID","IID"]+covars], ancestry[["IID"]+pcols],on=["IID"],how="inner")
    plink_pheno = plink[["FID","IID","phenotype"]]

    # Output files
    pheno_file_path = f"{args.phenotype}_pheno_plink.txt"
    covar_file_path = f"{args.phenotype}_covar_combined.txt"

    plink_pheno.to_csv(pheno_file_path, sep="\t", index=False)
    data.to_csv(covar_file_path, sep="\t", index=False)
    
    sys.stderr.write(f"Done converting {args.phenotype} to plink format.")

if __name__ == "__main__":
    plink_pheno,data = main()