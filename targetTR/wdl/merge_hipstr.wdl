version 1.0

workflow merge_hipstr {
    input {
        Array[File] vcfs
        Array[File] vcf_indexes
        String out_prefix
	      Int? merge_mem = 4
    }

    call mergestr {
        input : 
          vcfs=vcfs,
          vcf_indexes=vcf_indexes,
          out_prefix=out_prefix+".merged",
	        merge_mem=merge_mem
    }

    output {
       File outfile = mergestr.outvcf
    }

    meta {
      description: "Merge VCFs from multiple HipSTR runs"
    }
}

task mergestr {
  input {
    Array[File] vcfs
    Array[File] vcf_indexes
    String out_prefix
    Int total = length(vcfs)
    Int? merge_mem = 4
  }

  command <<<
      touch vcf.list
      FILEARRAY=(~{sep=' ' vcfs}) # Load array into bash variable
      for (( c = 0; c < ~{total}; c++ )) # bash array are 0-indexed ;)
      do
           vcf-validator && echo ${FILEARRAY[$c]} >> vcf.list
           vcf-validator || echo "Failed: " ${FILEARRAY[$c]}
      done
      mergeSTR --vcfs-list vcf.list --out ~{out_prefix}
  >>>
    
  runtime {
      docker: "gcr.io/ucsd-medicine-cast/trtools-mergestr-files:latest"
      memory: merge_mem +"GB"
  }

  output {
      File outvcf = "${out_prefix}.vcf"
  }
}
