version 1.0

workflow imputation {
    input {
        String vcf 
        String vcf_index 
        File ref_panel
        String out_prefix
        String GOOGLE_PROJECT = ""
        String GCS_OAUTH_TOKEN = ""
        Int? mem 
        Int? window_size 
        File samples_file 
	    File? regions_file
        Boolean subset_region = false
        Boolean beagle_region = false
        String? chrom
    }

    call subset_vcf {
    input:
        samples_file=samples_file,
        regions_file=regions_file,
        vcf=vcf,
        vcf_index=vcf_index,
        GOOGLE_PROJECT=GOOGLE_PROJECT,
        GCS_OAUTH_TOKEN=GCS_OAUTH_TOKEN,
        subset_region=subset_region,
        out_prefix=out_prefix
    }
    
    call index_vcf {
        input:
            vcf=subset_vcf.outfile
    }
    call beagle {
        input : 
          vcf=index_vcf.outvcf, 
          vcf_index=index_vcf.outvcf_index,
          ref_panel=ref_panel, 
          out_prefix=out_prefix,
          GOOGLE_PROJECT=GOOGLE_PROJECT,
          GCS_OAUTH_TOKEN=GCS_OAUTH_TOKEN,
          mem=mem,
          window_size=window_size,
          beagle_region=beagle_region,
          chrom=chrom
    }
    call sort_index_beagle {
        input :
            vcf=beagle.outfile
    }
    output {
        File outfile = sort_index_beagle.outvcf 
        File outfile_index = sort_index_beagle.outvcf_index
    }
    meta {
      description: "Run Beagle on a subset of samples on a single chromesome with default parameters"
    }
}


task subset_vcf {
    input {
        String vcf
        String vcf_index
        File samples_file
	    File? regions_file
        String GOOGLE_PROJECT = ""
        String GCS_OAUTH_TOKEN = ""
        String out_prefix=out_prefix
        Boolean subset_region = false
    }

    command <<<
        export GCS_REQUESTER_PAYS_PROJECT=~{GOOGLE_PROJECT}
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        # The "bcftools head" command was to check the header for the labeling if contigs e.g. chr21 vs 21.
        # bcftools head ~{vcf} > header.txt
        # Subsetting region for each chromesome
        if [[ "~{subset_region}" == true ]] ; then
        bcftools view -R ~{regions_file} -S ~{samples_file} ~{vcf} > ~{out_prefix}.vcf
        else:
            bcftools view -S ~{samples_file} ~{vcf} > ~{out_prefix}.vcf
        fi

    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
    }

    output {
        File outfile = "${out_prefix}.vcf"
    }    
}
task index_vcf {
    input {
      File vcf
    }

    String basename = basename(vcf, ".vcf")

    command <<<
        bgzip -c ~{vcf}> ~{basename}.vcf.gz && tabix -p vcf ~{basename}.vcf.gz
    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/vcfutils:latest"
    }

    output {
        File outvcf = "${basename}.vcf.gz"
        File outvcf_index = "${basename}.vcf.gz.tbi"
    }
}

task beagle {
    input {
        File vcf 
        File vcf_index 
        File ref_panel
        String out_prefix
        String GOOGLE_PROJECT = ""
        String GCS_OAUTH_TOKEN = ""
        Int? mem 
        Int? window_size
        Boolean beagle_region = false
        String? chrom
    } 

    command <<<

        export GCS_REQUESTER_PAYS_PROJECT=~{GOOGLE_PROJECT}
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        if [[ "~{beagle_region}" == false ]] ; then
        java -Xmx~{mem}g -jar /beagle.jar \
            gt=~{vcf} \
            ref=~{ref_panel} \
            window=~{window_size} \
            out=~{out_prefix}_output
        else
            java -Xmx~{mem}g -jar /beagle.jar \
            gt=~{vcf} \
            ref=~{ref_panel} \
            window=~{window_size} \
            chrom=~{chrom} \
            out=~{out_prefix}_output
        fi
    >>>
    
    #file upto 300mb use mem=25
    runtime {
        docker:"gcr.io/ucsd-medicine-cast/beagle:latest"
	    memory: mem + "GB"
    }

    output {
       File outfile = "${out_prefix}_output.vcf.gz"
    }
}

task sort_index_beagle {
    input {
      File vcf
    }

    String basename = basename(vcf, ".vcf.gz")

    command <<<
        zcat ~{vcf} | vcf-sort | bgzip -c > ~{basename}.sorted.vcf.gz && tabix -p vcf ~{basename}.sorted.vcf.gz
        echo "Number of TRs in the genotyped file"
        bcftools view -i 'ID="."' ~{basename}.sorted.vcf.gz | grep -v "^#" | wc -l
    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/vcfutils:latest"
    }

    output {
    File outvcf = "${basename}.sorted.vcf.gz"
    File outvcf_index = "${basename}.sorted.vcf.gz.tbi"
  }
}
