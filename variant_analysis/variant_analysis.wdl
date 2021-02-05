version 1.0

import "../wdl_structs.wdl"

task DeconstructsigPrep38 {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String pairName
        String highconfidencePath = "~{pairName}.annotated.highconfidence.txt"
        File mainVcf
        File header
    }

    command {
        bcftools query \
        -i "TYPE=\"snp\" && HighConfidence=1" \
        -f "~{pairName}\\t%CHROM\\t%POS\\t%REF\\t%ALT\\n" \
        -f "~{pairName}\\t%CHROM\\t%POS\\t%REF\\t%ALT\\n" \
        ~{mainVcf} \
        | grep -P "\\tchrX\\t" -v \
        | grep -P "\\tchrY\\t" -v \
        | cat ~{header} - \
        > ~{highconfidencePath} 
    }

    output {
        File highconfidence = "~{highconfidencePath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task Deconstructsig {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String pairName
        String suffix = "~{pairName}.deconstructSigs.v1.8.0.signatures.highconfidence"
        File deconstructsigsBs
        File deconstructsigsFasta
        File highconfidence
    }

    command {
        Rscript \
        run_deconstructSigs.v2.R \
        ~{highconfidence} \
        ~{suffix} \
        ~{deconstructsigsBs} \
        ~{deconstructsigsFasta} \
        plotting.R \
    }

    output {
        File diff = "~{pairName}.deconstructSigs.v1.8.0.signatures.highconfidence.diff.txt"
        File trinuc = "~{pairName}.deconstructSigs.v1.8.0.signatures.highconfidence.trinucleotidecounts.txt"
        File input_file = "~{pairName}.deconstructSigs.v1.8.0.signatures.highconfidence.input.txt"
        File highconfidencePng = "~{pairName}.deconstructSigs.v1.8.0.signatures.highconfidence.png"
        File highconfidenceTxt = "~{pairName}.deconstructSigs.v1.8.0.signatures.highconfidence.txt"
        File reconstructed = "~{pairName}.deconstructSigs.v1.8.0.signatures.highconfidence.reconstructed.txt"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}