version 1.0

import "../wdl_structs.wdl"


task Deconstructsig {
    input {
        String pairId
        String outputPrefix = "~{pairId}.cosmic.v3.2.deconstructSigs.signatures.highconfidence"
        File mainVcf
        String vepGenomeBuild
        String highConf = "TRUE"
        File cosmicSigs
        Int memoryGb = 8
        Int diskSize = 1
    }

    command {
        Rscript \
        /run_deconstructSigs.R \
        --highconf ~{highConf} \
        --file ~{mainVcf} \
        --ref ~{vepGenomeBuild} \
        --cosmic ~{cosmicSigs} \
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
        docker : "gcr.io/nygc-public/deconstructsigs:1.9.0.nygc.1"
    }
}
