#!/usr/bin/env python3
"""
Run phenotype preprocessing on a list of phenotypes.
Upload results to ${WORKSPACE_BUCKET}/phenotypes/
Output manifest file

Usage:
./process_phenotypes.py <phenotype_list.csv>

Phenotype list has columns:
* `phenotype`
* `concept_id`
* `units`
* `outlier_sd`
* `min`
* `max`
* `drugcovars`
"""

import hashlib
import os
import pandas as pd
import subprocess
import sys

try:
    ptfile = sys.argv[1]
except:
    sys.stderr.write(__doc__)
    sys.exit(1)

def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def RunCmd(cmd):
    res = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    res.wait()
    retcode = res.returncode
    if retcode != 0:
    	sys.stderr.write("Command failed: %s"%cmd)
    	sys.exit(1)

def UploadGCS(outfile, gcsfile):
	cmd = "gsutil cp {outfile} {gcsfile}".format(outfile=outfile, gcsfile=gcsfile)
	RunCmd(cmd)

manifest = open(ptfile.replace(".csv",".manifest.csv"), "w")
ptdata = pd.read_csv(ptfile)
for index, row in ptdata.iterrows():
    if row["phenotype"].startswith("#"): continue
    sys.stderr.write("Processing {phenotype}\n".format(phenotype=row["phenotype"]))
    cmd = """./aou_phenotype_preprocessing.py \
           --phenotype {phenotype} \
           --concept-id {concept} \
           --units \"{units}\" \
           --range {minval},{maxval} \
           --outlier-sd {outlier_sd}""".format(phenotype=row["phenotype"], \
           	concept=row["concept_id"], units=row["units"], \
            minval=row["min"], maxval=row["max"], \
           	outlier_sd=row["outlier_sd"])
    if str(row["drugcovars"]) != "nan":
        cmd += " --drugexposure-covariate-concept-ids {drug}".format(drug=row["drugcovars"])
    if str(row["ppi"]) == "ppi":
        cmd += " --ppi"

    outfile = row["phenotype"]+"_phenocovar.csv"
    gcsfile = os.path.join(os.environ["WORKSPACE_BUCKET"], "phenotypes", outfile)
    RunCmd(cmd)
    UploadGCS(outfile, gcsfile)
    numsamples = len(open(outfile, "r").readlines())-1
    outitems = [row["phenotype"], row["concept_id"], \
    	row["units"], row["outlier_sd"], row["min"], row["max"], row["drugcovars"], \
    	outfile, numsamples, md5(outfile)]
    manifest.write(",".join([str(item) for item in outitems])+"\n")

manifest.close()