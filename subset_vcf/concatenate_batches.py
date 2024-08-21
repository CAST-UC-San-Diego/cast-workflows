#!/usr/bin/env python3
"""
Usage:

./concanate_batches.py <cromshell_job_output.json> <outprefix> <bcftools_path>

"""

DEBUG = False

import json
import os
import subprocess
import sys

jobdata = json.load(open(sys.argv[1], "r"))
outprefix = sys.argv[2]
bcftools_path = sys.argv[3]

def SortByCoordinate(vcf_files):
	""" Sort files by coordinates in filename """
	# Filename syntax: bucket/batch99-chr11-50000000-60000000.vcf.gz
	files = [] # (startcoord, file)
	for f in vcf_files:
		coord = int(os.path.basename(f).split("-")[2])
		files.append((coord, f))
	files = sorted(files, key = lambda x: x[0])
	return [item[1] for item in files]

# Gather files for each batch
batch_files = {} # batchname -> {"vcf": [], "index": []}
for region in jobdata["subset_vcf.vcf_output_array"]:
	for f in region:
		batchname = os.path.basename(f).split("-")[0]
		if batchname not in batch_files.keys():
			batch_files[batchname] = {}
			batch_files[batchname]["vcf"] = []
			batch_files[batchname]["index"] = []
		batch_files[batchname]["vcf"].append(f)
for region in jobdata["subset_vcf.vcf_index_array"]:
	for f in region:
		batchname = os.path.basename(f).split("-")[0]
		batch_files[batchname]["index"].append(f)

os.environ["GCS_REQUESTER_PAYS_PROJECT"] = os.environ["GOOGLE_PROJECT"]
# Process one batch at a time
for batch in batch_files.keys():
	vcf_files = batch_files[batch]["vcf"]
	index_files = batch_files[batch]["index"]
	output_fname = "%s-%s.vcf.gz"%(outprefix, batch)
	print("##### Processing %s ######"%batch)

	cmds = []

	# Run the bcftools concat command
	cmd = "%s concat %s -Oz -o %s"%(bcftools_path, " ".join(SortByCoordinate(vcf_files)), output_fname)
	cmd += " && tabix -p vcf %s"%(output_fname)
	cmds.append(cmd)

	# Upload to our workspace bucket
	cmd = "gsutil cp %s ${WORKSPACE_BUCKET}/acaf_batches/%s/%s"%(output_fname, outprefix, output_fname)
	cmds.append(cmd)

	cmd = "gsutil cp %s.tbi ${WORKSPACE_BUCKET}/acaf_batches/%s/%s.tbi"%(output_fname, outprefix, output_fname)
	cmds.append(cmd)

	# Remove this batch before we move on
	cmd = "rm %s; rm %s.tbi"%(output_fname, output_fname)
	cmds.append(cmd)

	if DEBUG:
		print(cmds)
	else:
		# Refresh credentials
		token_fetch_command = subprocess.run(['gcloud', 'auth', 'application-default', 'print-access-token'], \
			capture_output=True, check=True, encoding='utf-8')
		token = str.strip(token_fetch_command.stdout)
		os.environ["GCS_OAUTH_TOKEN"] = token

		for cmd in cmds:
			output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read()
			print(output.decode("utf-8"))
