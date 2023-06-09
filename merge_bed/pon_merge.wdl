version 1.0

import "../wdl_structs.wdl"

task MakeSinglePon {
    input {
        Int threads = 1
        Int memoryGb = 40
        Int diskSize
        String sampleId
        File vcfAnnotatedVep
        String simpleTable = "~{sampleId}.pon.table.bed"
        String singlePonPath = "~{sampleId}.pon.bed"
        File countPonSites = "gs://nygc-comp-s-fd4e-input/scripts/count_pon_sites.py"
        String openBraket = "{"
        String closeBraket = "}"
    }
    
    command {
        echo "Start header..."
        cat ~{vcfAnnotatedVep} \
        | awk '~{openBraket}if(/^#/)print;else exit~{closeBraket}' \
        | grep "^#CHROM" \
        | cut -f 1-7 \
        > ~{simpleTable}
        
        set -e -o pipefail
        
        echo "Finish table..."
        grep -v "^#" ~{vcfAnnotatedVep} \
        | cut -f 1-7 \
        >> ~{simpleTable}
        
        echo "Refine table..."
        python \
        ~{countPonSites} \
        ~{simpleTable} \
        ~{singlePonPath}
        
    }

    output {
        File singlePon = "~{singlePonPath}"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:20a48e2c422a43ce35e197243bda8dbf06c9a7b3175094524f74f8835cce85b6"
    }
}

task MergePons {
    input {
        Int threads = 1
        Int memoryGb = 40
        Int diskSize = 300
        String name = "WGS_1000g_GRCh38"
        Array[File] singlePons
        String ponPath = "~{name}.pon.bed"
        File mergePonSites = "gs://nygc-comp-s-fd4e-input/scripts/merge_pon_sites.py"
    }
    
    command {        
        python \
        ~{mergePonSites} \
        ~{ponPath} \
        ~{sep=" " singlePons}
        
    }

    output {
        File pon = "~{ponPath}"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:20a48e2c422a43ce35e197243bda8dbf06c9a7b3175094524f74f8835cce85b6"
    }
}

