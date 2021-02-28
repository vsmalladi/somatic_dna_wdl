version 1.0

import "../wdl_structs.wdl"

task DeconstructsigPrep38 {
    input {
        String pairId
        String highconfidencePath = "~{pairId}.annotated.highconfidence.txt"
        File mainVcf
        Int memoryGb = 8
        Int diskSize = 1
    }

    command {
        set -e -o pipefail
         
        echo -e "Sample\tchr\tpos\tref\talt\n" > header.txt
        
        bcftools query \
        -i "TYPE=\"snp\" && HighConfidence=1" \
        -f "~{pairId}\\t%CHROM\\t%POS\\t%REF\\t%ALT\\n" \
        -f "~{pairId}\\t%CHROM\\t%POS\\t%REF\\t%ALT\\n" \
        ~{mainVcf} \
        | grep -P "\\tchrX\\t" -v \
        | grep -P "\\tchrY\\t" -v \
        | cat header.txt - \
        > ~{highconfidencePath} 
    }

    output {
        File highconfidence = "~{highconfidencePath}"
    }

    runtime {
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/bcftools:1.5"
    }
}

task Deconstructsig {
    input {
        Int memoryGb = 8
        Int diskSize = 1
        String pairId
        String suffix = "~{pairId}.deconstructSigs.v1.8.0.signatures.highconfidence"
        File deconstructsigsBs
        File deconstructsigsFasta
        File highconfidence
        File plotting_r = "gs://nygc-comp-s-fd4e-resources/plotting.R"
    }

    command {
        Rscript \
        /run_deconstructSigs.v2.R \
        ~{highconfidence} \
        ~{suffix} \
        ~{deconstructsigsBs} \
        ~{deconstructsigsFasta} \
        ~{plotting_r}
    }

    output {
        File diff = "~{pairId}.deconstructSigs.v1.8.0.signatures.highconfidence.diff.txt"
        File trinuc = "~{pairId}.deconstructSigs.v1.8.0.signatures.highconfidence.trinucleotidecounts.txt"
        File input_file = "~{pairId}.deconstructSigs.v1.8.0.signatures.highconfidence.input.txt"
        File highconfidencePng = "~{pairId}.deconstructSigs.v1.8.0.signatures.highconfidence.png"
        File highconfidenceTxt = "~{pairId}.deconstructSigs.v1.8.0.signatures.highconfidence.txt"
        File reconstructed = "~{pairId}.deconstructSigs.v1.8.0.signatures.highconfidence.reconstructed.txt"
    }

    runtime {
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/somatic_tools:0.9.2"
    }
}