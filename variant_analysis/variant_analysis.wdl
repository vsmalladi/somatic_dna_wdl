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
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/deconstructsigs@sha256:009ddb6ed3ec2a0290a88b1e7027dd3caac1a2f5f3df3e8f68f410481d9323a3"
    }
}
