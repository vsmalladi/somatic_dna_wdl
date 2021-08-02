version 1.0

import "../wdl_structs.wdl"


task Deconstructsig {
    input {
        String pairId
        String outputPrefix = "~{pairId}.cosmic.v3.2.deconstructSigs.v1.9.0.signatures.highconfidence"
        File mainVcf
        String vepGenomeBuild 
        File runDeconstructSigs = "gs://nygc-comp-s-fd4e-input/run_deconstructSigs.R"
        File deconstructSigs
        Int memoryGb = 8
        Int diskSize = 1
    }

    command {
        Rscript \
        ~{runDeconstructSigs} \
        --highconf \
        --file ~{mainVcf} \
        --ref ~{vepGenomeBuild} \
        --cosmic ~{deconstructSigs} \
        --output ~{outputPrefix} \
        --samplename ~{pairId}
    }

    output {
       File sigs = "~{outputPrefix}.txt"
       File counts = "~{outputPrefix}.counts.txt"
       File sigInput = "~{outputPrefix}.input.txt"
       File reconstructed = "~{outputPrefix}.reconstructed.txt"
       File diff = "~{outputPrefix}.diff.txt"
    }

    runtime {
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/deconstructsigs:1.9.0"
    }
}
