# Running GWAS in All of Us

**under construction**

**Note: all of this needs to be run on the AoU workbench.**

## Quickstart

Assuming the phenotype you are analyzing has already been preprocessed (see below), you can run GWAS (by default using Hail):

1. Example test on a small region (can be run on default Hail environment)
```
./aou_gwas.py --phenotype ALT --num-pcs 10 --region chr11:119206339-119308149
```
will output `ALT_hail_ALL_chr11_119206339_119308149.gwas.tab` 

2. Example full GWAS. (TODO: which Hail compute resources)
```
./aou_gwas.py --phenotype ALT --num-pcs 10
```
will output `ALT_hail_ALL.gwas.tab`

Finding the appropriate spark parameters can be challenging/expensive. Instead, we recommend using the batched version instead of the full GWAS (see instructions below).

3. Small example associaTR run (requires installing TRTools - see setup note below)

```
./aou_gwas.py --phenotype platelet_count --num-pcs 10 --method associaTR --tr-vcf ${WORKSPACE_BUCKET}/cromwell-execution/targetTR/4519d903-dde6-4c8e-b1ee-bcc2d7cd6dd7/call-sort_index/CBL_test.filtered.sorted.vcf.gz
```
will output `platelet_count_associaTR_ALL.gwas.tab`

Note: currently need to install TRTools from branch `associatr-updates` due to conflict with newer numpy version.

## Detailed usage

Required arguments:

* `--phenotype <STR>`: phenotype ID to use to name the output file

Optional arguments:

* `--method <STR>`: Which GWAS method to use. Default: hail
* `--samples <FILE>`: csv file with list of samples to restrict to. Needs columns "person_id" and "sex_at_birth_Male". Typically this would a list of samples that passed upstream sample-level QC info. Defaults to `${WORKSPACE_BUCKET}/samples/passing_samples_v7.1.csv`.
* `--ancestry-pred-path <FILE>`: Path to file with PC info. Defaults to the AoU v7 ancestry predictions.
* `--region <STR>`: Region to restrict to (chrom:start-endd)
* `--num-pcs <INT>`: Number of population PCs to include as covariates. Default: 10
* `--ptcovars <STR>`: Comma-separated list of phenotype-specific covariates. Default: age
* `--sharedcovars <STR>`: Comma-separated list of shared covariates (besides PCs). Default: sex_at_birth_Male
* `--plot`: Output Manhattan and QQ plots
* `--norm`: Quantile or zscore normalize phenotype, type=str
* `--norm-by-sex`: Apply normalization by sex
* `--sample-call-rate`: Apply min sample call rate QC. Default: 0.9
* `--variant-call-rate`: Apply min sample call rate QC. Default: 0.9 
* `--MAF`: Apply minor allele frequency QC. Default=0.01
* `--HWE`: Apply HWE p-value cutoff QC. Default=1e-15
* `--GQ`: Apply minimun genotype score QC. Default=20
* `--logistic`: Apply logistic regression, default is linear regression
* `--test`: Select test for logistic regression: wald,lrt,score,firth



## Important file locations:

* Samples passing QC: `${WORKSPACE_BUCKET}/samples/passing_samples_v7.1.csv`
* Preprocessed phenotypes+covariates: `${WORKSPACE_BUCKET}/phenotypes/${phenotype}_phenocovar.csv`
* GWAS results: `${WORKSPACE_BUCKET}/gwas/${phenotype}_hail_${cohort}.gwas.tab` where `cohort` is parsed from the `--samples` argument. If using the default for `--samples`, `cohort` will be `ALL`.
* Manifest of existing preprocessed phenotypes: `phenotypes_manifest.csv` (this repo)

Available cohort options (edit 4/24/24: updated for v7.1):
* `${WORKSPACE_BUCKET}/samples/EUR_WHITE.csv` (119184 samples) (self-reported race matches `ANCESTRY_PRED`)
* `${WORKSPACE_BUCKET}/samples/AFR_BLACK.csv` (44386 samples) (self-reported race matches `ANCESTRY_PRED`)
* `${WORKSPACE_BUCKET}/samples/NOT_AFR_BLACK.csv` (174828 samples) (all samples except where race is African American/Black or `ANCESTRY_PRED` is AFR)

## Setup

The GWAS instructions below assume you have:

1. Run sample preprocessing to get a set of samples passing all filters. See the `sample_preprocessing` folder for more details.

We will store the list of samples passing QC at:

```
${WORKSPACE_BUCKET}/samples/passing_samples_v7.1.csv
```

2. Run phenotype preprocessing on your phenotype of interest. See the `phenotype_preprocessing` folder for more details.

We will store phenotypes+shared covariates at:
```
${WORKSPACE_BUCKET}/phenotypes/${phenotype}_phenocovar.csv
```

3. Code setup

To get this repository:
```
git clone https://github.com/cast-genomics/cast-workflows
cd cast-workflows
git checkout gwas-aou
```

For running associaTR, you need to install a special TRTools branch for now:

```
cd
git clone https://github.com/gymrek-lab/TRTools
cd TRTools
git checkout associatr-updates
pip install -e .
```

## Batching the full GWAS

Follow the instructions above to determine the command you would run for a full gwas. Write this command to a file (e.g. `full_gwas_cmd.txt`). Then run the following:

`python batch_gwas.py <phenotype_name> <full_gwas_cmd _file_name>`

OR

`nohup python batch_gwas.py <phenotype_name> <full_gwas_cmd _file_name> > nohup.log 2>&1 &`

This will run a gwas per chromosome. The per chromosome gwas results will get stored at `${WORKSPACE_BUCKET}/gwas/<phenotype_name>`. You will need a spark cluster to run the batched gwas. We have found that the following spark parameters are sufficient (but the full batched gwas will still take on the order of 12 hours to complete):

<img width="687" alt="Screen Shot 2024-03-08 at 4 10 44 PM" src="https://github.com/CAST-genomics/cast-workflows/assets/16807372/4b22fa94-efbc-4c7b-974b-dab86d43db48">


The next step is to combine per chromosome results. This does not require a spark cluster, a regular compute environment is sufficient. First, determine the name of the gwas files that got written. You can do this by running the following:

`gsutil ls ${WORKSPACE_BUCKET}/gwas/<phenotype_name>`

Look for any file ending in `.gwas.tab`. Copy the file name up until the region (i.e. from the start up to but excluding 'chr'). For example, if the file is `ldl_cholesterol_hail_EUR_WHITE_chr1_1_248956422.gwas.tab`, the prefix you should copy is `ldl_cholesterol_hail_EUR_WHITE`.

Then run the following:
`python combine_gwas_batches.py <phenotype_name> <gwas_file_prefix>`

This will write out the gwas results across all chromosomes to a single file and copy it to `${WORKSPACE_BUCKET}/gwas/<phenotype_name>` as if you had run the full gwas described above (quickstart point 2). Note that the gwas file is likely to be on the order of 1GB and will take 5-10 min to copy over to the bucket.
