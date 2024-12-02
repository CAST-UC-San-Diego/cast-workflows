version 1.0

workflow tr_gwas {
    input {
        Array[File] pgens = [] 
        Array[File] psams = [] 
        Array[File] pvars = [] 
        Array[File] phenotypes = []
        Array[File] cohorts = []
        String GOOGLE_PROJECT = ""
        String GCS_OAUTH_TOKEN = ""
        String WORKSPACE_BUCKET = ""
    }

    ### Separate workflow for each phenotype ###

    scatter (pheno_file in phenotypes) {      
       
        call convert_phenotype {
            input:
                pheno=pheno_file,
                GOOGLE_PROJECT=GOOGLE_PROJECT,
                GCS_OAUTH_TOKEN=GCS_OAUTH_TOKEN,
                WORKSPACE_BUCKET=WORKSPACE_BUCKET
        }

            scatter (cohort in cohorts) {
                
                call run_tr_gwas {
                    input:
                        pgens=pgens,
                        psams=psams,
                        pvars=pvars,
                        pheno=convert_phenotype.outfile_pheno,
                        covar=convert_phenotype.outfile_covar,
                        samples=cohort,
                        out_prefix="${convert_phenotype.outpheno_name}"
                       
                }
            }
    }   
    

    output {
        Array[Array[File]] gwas_outputs = run_tr_gwas.outfile
    }

    meta {
        description: "Use plink2 to run genome-wide TR GWAS. It runs GWAS on each chromosome for each phenotype and merges summary statstics across all chromosomes."
    }
}

task convert_phenotype {
    input {
        File pheno
        String GOOGLE_PROJECT = ""
        String GCS_OAUTH_TOKEN = ""
        String WORKSPACE_BUCKET = ""
    }
    
    String pheno_name = basename(pheno,"_phenocovar.csv")

    command <<<
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        export GCS_REQUESTER_PAYS_PROJECT="~{GOOGLE_PROJECT}"
        export WORKSPACE_BUCKET="~{WORKSPACE_BUCKET}"
        python3 /usr/bin/convert_phenotype_plink.py --phenotype ~{pheno_name} 
    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/plink2:latest"
        preemptible: 1
    }

    output {
        File outfile_pheno = "${pheno_name}_pheno_plink.txt"
        File outfile_covar = "${pheno_name}_covar_combined.txt"
        String outpheno_name = "${pheno_name}"
    }
}

task run_tr_gwas {
    input {
        Array[File] pgens
        Array[File] psams
        Array[File] pvars
        File pheno
        File covar 
        File samples
        String out_prefix
        Int total = length(pgens)
    }

    String sample_name = basename(samples,"_plink.txt")

    command <<<

        #check if cohort is EUR or non_AFR, exclude PC6
        
        echo ~{sample_name} | grep -q -E "EUR|NOT_AFR"
        if [ $? -eq 0 ]; then
            cut --complement -f 10 ~{covar} > "~{sample_name}_~{out_prefix}_covar.txt"
            covar_file="~{sample_name}_~{out_prefix}_covar.txt"

        else 
            covar_file=~{covar}           
        fi

        # Run GWAS on each chrom
        PFILEARRAY=(~{sep=" " pgens})
        gwas_outfiles=""
        # bash array are 0-indexed 
        for (( c = 0; c < ~{total}; c++ ))
        do
            pfile=${PFILEARRAY[$c]}
            pfile_outprefix="${pfile%.pgen}"
            chrom_outprefix=$(basename $pfile .pgen)
            plink2 --pfile ${pfile_outprefix} \
               --pheno ~{pheno} \
               --linear \
               --covar ${covar_file} \
               --keep ~{samples} \
               --covar-variance-standardize \
               --out "~{out_prefix}_${chrom_outprefix}_~{sample_name}"
               
            gwas_outfiles+="~{out_prefix}_${chrom_outprefix}_~{sample_name}.phenotype.glm.linear "

        done
        
        # Concatenate all results
        head -n 1 ${gwas_outfiles} > "~{out_prefix}_~{sample_name}_gwas.tab"

        for file in ${gwas_outfiles}
        do
            tail -n +2 "$file" >> "~{out_prefix}_~{sample_name}_gwas.tab"
        done
    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/plink2:latest"
        preemptible: 1
        memory: "4G"
        disks: "local-disk 330 SSD"
    }

    output {
       File outfile = "${out_prefix}_~{sample_name}_gwas.tab"
    }
}