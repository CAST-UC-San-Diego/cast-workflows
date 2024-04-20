chr=$1
from=$2
to=$3
sp=$4
fp=${sp#*.}
s1_samples_file=$5
s1_samples="${s1_samples_file%.*}"
s2_samples_file=$6
s2_samples="${s2_samples_file%.*}"
phen=$7 #e.g. ldl_cholesterol
#cojo_numsignals=$4
#cr_numsignals=$4

fsl_s1=0.0001
fsl_s2=0.0001
maf=0.01
min_num_snps=25
plink_file_prefix=gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/exome/plink_bed/acaf_threshold.chr${chr}
plink_file_prefix=gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/acaf_threshold_v7.1/plink_bed/acaf_threshold.chr${chr}

logfile=pipsort.log
results=results_${phen}_${fp}.txt

#0. copy over samples files if needed
if [ -e "$s1_samples_file" ]; then
    echo "sample s1 file already exists"
else
    gsutil cp "${WORKSPACE_BUCKET}/samples/${s1_samples_file}" ./
fi
if [ -e "$s2_samples_file" ]; then
    echo "sample s2 file already exists"
else
    gsutil cp "${WORKSPACE_BUCKET}/samples/${s2_samples_file}" ./
fi

#1. get plink files
gsutil -q stat ${WORKSPACE_BUCKET}/pipsort/plink/${chr}_${from}_${to}_${phen}_plink.bed
status=$?
if [[ $status == 0 ]]; then
    echo "plink file exists"
else
    echo "plink file does not exist $chr $to $from $phen" >> $logfile
    exit 1
fi
gsutil cp "${WORKSPACE_BUCKET}/pipsort/plink/${chr}_${from}_${to}_${phen}_plink.*" ./


#2. get phenotype file for gwas
if [ -e "${phen}_phenocovar.csv" ]; then
    echo "no need to copy phen data"
else
    gsutil cp "${WORKSPACE_BUCKET}/pipsort/phenotypes/${phen}_phenocovar.csv" ./
fi
python convert_phen_to_plink_format.py ${phen}_phenocovar.csv ${phen}_plink_format.tab

#3. get gwas data and add rsid to gwas file
gwas_file_s1_pre=${phen}_hail_${s1_samples}
gwas_file_s2_pre=${phen}_hail_${s2_samples}

if [ -e "$gwas_file_s1_pre}.gwas.tab" ]; then
    echo "no need to copy gwas1"
else
    gsutil cp "${WORKSPACE_BUCKET}/pipsort/gwas/${gwas_file_s1_pre}.gwas.tab" ./
fi
if [ -e "$gwas_file_s2_pre}.gwas.tab" ]; then
    echo "no need to copy gwas2"
else
    gsutil cp "${WORKSPACE_BUCKET}/pipsort/gwas/${gwas_file_s2_pre}.gwas.tab" ./
fi


python subset_gwas_to_loci.py $chr $from $to ${gwas_file_s1_pre}.gwas.tab ${gwas_file_s1_pre}_${chr}_${from}_${to}.gwas.tab
python subset_gwas_to_loci.py $chr $from $to ${gwas_file_s2_pre}.gwas.tab ${gwas_file_s2_pre}_${chr}_${from}_${to}.gwas.tab
python add_rsid_col.py ${gwas_file_s1_pre}_${chr}_${from}_${to}.gwas.tab
python add_rsid_col.py ${gwas_file_s2_pre}_${chr}_${from}_${to}.gwas.tab


#3. subset to snps and samples I need
python convert_samples_to_plink_format.py $s1_samples_file plink_${s1_samples_file}
python convert_samples_to_plink_format.py $s2_samples_file plink_${s2_samples_file}


python extract_all_snps.py ${gwas_file_s1_pre}_${chr}_${from}_${to}.gwas.tab ${gwas_file_s2_pre}_${chr}_${from}_${to}.gwas.tab union_snps.txt
plink2 --bfile ${chr}_${from}_${to}_${phen}_plink --chr $chr --keep plink_$s1_samples_file --make-bed --out s1_data --pheno ${phen}_plink_format.tab --prune --extract union_snps.txt
plink2 --bfile ${chr}_${from}_${to}_${phen}_plink --chr $chr --keep plink_$s2_samples_file --make-bed --out s2_data --pheno ${phen}_plink_format.tab --prune --extract union_snps.txt

python sort_gwas_results.py --infile ${gwas_file_s1_pre}_${chr}_${from}_${to}.gwas.tab --pos_col pos --rsid_col rsid
python sort_gwas_results.py --infile ${gwas_file_s2_pre}_${chr}_${from}_${to}.gwas.tab --pos_col pos --rsid_col rsid


#exclude high p vals, outputs s1_snps.txt and s2_snps.txt and s1_gwas_temp.txt and s2_gwas_temp.txt
python remove_high_pvals_not_common.py --gwas1 ${gwas_file_s1_pre}_${chr}_${from}_${to}.gwas.tab --gwas2 ${gwas_file_s2_pre}_${chr}_${from}_${to}.gwas.tab --fsl1 $fsl_s1 --fsl2 $fsl_s2 --pval_col p_value --rsid_col rsid

python get_gwas_data_in_mscaviar_style.py s1_gwas_temp.txt s1_ldl_gwas.txt
python get_gwas_data_in_mscaviar_style.py s2_gwas_temp.txt s2_ldl_gwas.txt
