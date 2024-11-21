#!/usr/bin/env python3

"""
Preprocess AoU phenotypes.
Output cleaned phenotype and covariates files
for the specified phenotype

Requires these environment variables to be set:
WORKSPACE_CDR
WORKSPACE_BUCKET
"""

import argparse
import aou_queries
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

SAMPLEFILE = os.path.join(os.environ["WORKSPACE_BUCKET"], "samples", \
    "passing_samples_v7.csv")
WORKSPACE_BUCKET = os.getenv('WORKSPACE_BUCKET')
WORKSPACE_CDR = os.getenv('WORKSPACE_CDR')

def MSG(msg_str):
    """
    Write a helpful progress message to stderr

    Arguments
    ---------
    msg_str : str
       Message to write
    """
    sys.stderr.write("[aou_phenotype_preprocessing]: %s\n"%msg_str.strip())

def ERROR(msg_str):
    """
    Write an error message to stderr then quit
    with exit code 1

    Arguments
    ---------
    msg_str : str
       Message to write
    """
    MSG("ERROR: " + msg_str)
    sys.exit(1)

def SQLToDF(sql):
    """
    Extract Google bigquery results to pandas df

    Arguments
    ---------
    sql : str
       Query string

    Returns
    -------
    df : pandas.DataFrame
       Dataframe with the results
    """
    df = pd.read_gbq(
        sql,
        dialect="standard",
        use_bqstorage_api=("BIGQUERY_STORAGE_API_ENABLED" in os.environ),
        progress_bar_type="tqdm_notebook")
    return df

def OverlapDrugMeasurement(measurement_datetime, start, end):
    """
    Determine whether a phenotype measurement overlaps
    the start/end time of a drug

    Currently: counts as no if phenotype measurement
    is before drug start

    Arguments
    ---------
    measurement_datetime : datetime
       Time phenotype measurement taken
    start : datetime
       Time drug started
    end : datetime
       Time drug ended

    Returns
    -------
    res : int
       0 (not on drug) or 1 (on drug)
    """
    if measurement_datetime < start:
        res = 0
    else:
        res = 1
    return res

def GetDrugData(concept_id):
    """
    Extract drug data for a concept id

    Arguments
    ---------
    concept_id : int
       AoU concept ID

    Returns
    -------
    drugdata : pandas.DataFrame
       Includes columns: person_id, start (drug start), end (drug end)
    """
    drug_sql = aou_queries.ConstructDrugExposureSQL(concept_id)
    drugdata = SQLToDF(drug_sql)
    drugdata = drugdata.groupby(["person_id"]).agg(start=('drug_exposure_start_datetime', np.min), \
        end=('drug_exposure_end_datetime', np.max)).reset_index()
    return drugdata


def my_median(series):
    """
    compute median of median 
    odd number use np.median
    even number pick either one

    Arguments
    ---------
    series : pandas.DataFrame column


    Returns
    -------
    median: pandas.DataFrame
       
    """
    if len(series)%2==1:
        return np.median(series)
    else:
        my_median = sorted(series)[int(len(series)/2)]
        return my_median

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--phenotype", help="Phenotype ID", type=str, required=True)
    parser.add_argument("--samples", help="List of sample IDs,sex to keep", type=str, default=SAMPLEFILE)
    parser.add_argument("--concept-id", help="Concept ID for phenotype", type=str, required=True)
    parser.add_argument("--drugexposure-covariate-concept-ids", help="Comma-separated list of conceptid:conceptname to use as drug exposure covariates", type=str)
    parser.add_argument("--units", help="Comma-separated list of acceptable units. Accepted shorthands: blood", type=str, required=True)
    parser.add_argument("--range", help="min, max acceptable phenotype values", type=str)
    parser.add_argument("--ppi", help="Whether or not the phenotype is from the PPI measurements " + \
                                      "as opposed to LOINC. Only physical measurements (e.g. height) " + \
                                      "are in PPI.",
                                      action="store_true", default=False)
    parser.add_argument("--outlier-sd", help="filter samples with phenotype values exceeding this number of SDs",
                                        type=int, required=False)

    args = parser.parse_args()
    MSG("Processing %s"%args.phenotype)

    # Set up dataframes
    demog = SQLToDF(aou_queries.demographics_sql)
    ptdata = SQLToDF(aou_queries.ConstructTraitSQL(args.concept_id, args.ppi))
    data = pd.merge(ptdata, demog, on="person_id", how="inner")
    MSG("After merge, have %s data points"%data.shape[0])

    # Restrict to samples we want to keep
    sampfile = args.samples
    if sampfile.startswith("gs://"):
        sampfile = sampfile.split("/")[-1]
        if not os.path.isfile(sampfile):
            os.system("gsutil -u ${GOOGLE_PROJECT} cp %s ."%(args.samples))
    samples = pd.read_csv(sampfile)
    data = pd.merge(data, samples)
    MSG("After filter samples, have %s data points"%data.shape[0])

    # Filtering
    data.dropna(axis=0, subset=['value_as_number'],inplace=True)
    MSG("After filter NA, have %s data points"%data.shape[0])
    MSG("  Allowable units: %s"%str(aou_queries.GetUnits(args.units)))
    MSG("  Unique units observed: %s"%(str(set(data["unit_concept_name"]))))
    data = data[data["unit_concept_name"].isin(aou_queries.GetUnits(args.units))]
    MSG("After filter units, have %s data points"%data.shape[0])
    if args.range is not None:
        minval, maxval = aou_queries.GetPhenotypeRange(args.range)
        data = data[(data["value_as_number"]>minval) & (data["value_as_number"]<maxval)]
        MSG("After filter range, have %s data points"%data.shape[0])

    # Filter outlier values based on number of SDs
    num_sds = getattr(args, "outlier_sd", None)
    if num_sds is not  None:
        avg = data["value_as_number"].mean()
        sd = data["value_as_number"].std()
        minval = avg - num_sds * sd
        maxval = avg + num_sds * sd
        data = data[(data["value_as_number"] >= minval) & (data["value_as_number"] <= maxval)]
        MSG("After outlier filtering, have %s data points"%data.shape[0])

    # Determine a single representative value per person
    data['Year'] = data['measurement_datetime'].dt.year
    median_per_year = data.groupby(['person_id','Year']).agg(median_year=('value_as_number', my_median)).reset_index()
    median_of_medians = median_per_year.groupby(['person_id']).agg(median_median=('median_year', my_median)).reset_index()
    median_of_medians.rename({"median_median": "value_as_number"}, inplace=True, axis=1)

    # Merge back to whole dataframe to only keep median of median value per person
    filtered = pd.merge(data, median_of_medians, on=["person_id", "value_as_number"])
    MSG("After merge median medians, have %s data points"%filtered.shape[0])

    # De-duplicate to keep one entry per person
    filtered = filtered.sort_values("measurement_datetime").drop_duplicates(subset=["person_id"], keep="last")
    MSG("After dedup, have %s data points"%filtered.shape[0])

    # Output histogram of phenotype values
    plt.hist(filtered["value_as_number"])
    plt.savefig(args.phenotype+"_histogram.png", dpi=300)

    # Record age info
    filtered["age"] = filtered['measurement_datetime'].dt.year - \
        filtered["date_of_birth"].dt.year

    # Optionally add any trait specific drugexposure covariates
    covar_cols = []
    if args.drugexposure_covariate_concept_ids is not None:
        covar_concepts = args.drugexposure_covariate_concept_ids.strip().split(",")
        for ci in covar_concepts:
            concept_id, concept_name = ci.split(":")
            drugdata = GetDrugData(concept_id)
            drugdata = pd.merge(drugdata, filtered[["person_id", "measurement_datetime"]], on="person_id")
            drugdata[concept_name] = drugdata.apply(lambda x: \
                OverlapDrugMeasurement(x["measurement_datetime"], x["start"], x["end"]), 1)
            filtered = pd.merge(filtered, drugdata[["person_id", concept_name]], \
                on="person_id", how="left").fillna(value={concept_name: 0})
            MSG("After add %s, filtered has %s data points"%(concept_name, filtered.shape[0]))
            covar_cols.append(concept_name)

    # Output final phenotype value
    MSG("Final file has %s data points"%filtered.shape[0])
    filtered.rename({"value_as_number": "phenotype"}, inplace=True, axis=1)
    filtered[["person_id", "phenotype", "age", "sex_at_birth_Male"]+covar_cols].to_csv(args.phenotype+"_phenocovar.csv", index=False)

if __name__ == "__main__":
    main()

