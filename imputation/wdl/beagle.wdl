version 1.0

workflow beagle {
    input {
        File vcf
        File vcf_index
        File ref_panel
        File ref_panel_index
        String out_prefix
    }

    call beagle {
        input : 
          vcf=vcf, 
          vcf_index=vcf_index,
          ref_panel=ref_panel, 
          ref_panel_index=ref_panel_index,
          out_prefix=out_prefix
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
      description: "Run Beagle on a single chromesome with default parameters"
    }
}

task beagle {
    input {
        File vcf
        File vcf_index
        File ref_panel
        File ref_panel_index
        String out_prefix
    } 

    command <<<
      #set tokens for AoU
      export GCS_REQUESTER_PAYS_PROJECT=~{GOOGLE_PROJECT}
      export GCS_OAUTH_TOKEN=~{GCS_OAUTH_TOKEN}

      java -Xmx4g -jar beagle.version.jar \
            gt=~{vcf} \
            ref=~{ref_panel} \
            out=~{out_prefix}_imputed_TR_SNPs
    >>>
    
    #need to create new docker
    runtime {
        docker:"gcr.io/ucsd-medicine-cast/hipstr-gymreklab:latest"
    }

    output {
       File outfile = "${out_prefix}.vcf.gz"
    }
}

task sort_index_beagle {
    input {
      File vcf
    }

    String basename = basename(vcf, ".vcf.gz")

    command <<<
        zcat ~{vcf} | vcf-sort | bgzip -c > ~{basename}.sorted.vcf.gz && tabix -p vcf ~{basename}.sorted.vcf.gz
    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/vcfutils:latest"
    }

    output {
    File outvcf = "${basename}.sorted.vcf.gz"
    File outvcf_index = "${basename}.sorted.vcf.gz.tbi"
  }
}