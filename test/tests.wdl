version 1.0

import "../wdl_structs.wdl"

task SummarizeVcf {
    input {
        String pairId
        String name
        File oldOnlyVcf
        File newOnlyVcf
        File concordantVcf
        String detailedOutputTablePath = "~{pairId}.~{name}.detailed.output.csv"
        String summaryOutputTablePath = "~{pairId}.~{name}.counts.csv"
        
        Int diskSize = 10
        Int memoryGb = 20
        File compareVcfs = "gs://nygc-comp-s-fd4e-input/compare_vcfs.py"
    }
    
    command {
       python \
        ~{compareVcfs} \
        --old-only-vcf ~{oldOnlyVcf} \
        --new-only-vcf ~{newOnlyVcf} \
        --concordant-vcf ~{concordantVcf} \
        --output ~{detailedOutputTablePath} \
        --summary ~{summaryOutputTablePath}
    }
    
    output {
        File detailedOutputTable = "~{detailedOutputTablePath}"
        File summaryOutputTable = "~{summaryOutputTablePath}"
    }
    
    runtime {
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-internal-tools/somatic_tools:v1.1.1"
    }
    
    meta {
        summaryTsvColumns : ["Callset", "Confidence", "Count"]
        Callsets : ["Old unique", "New unique", "Concordant"]
    }
}

task SomPy {
    input {
        String pairId
        File oldVcf
        File newVcf
        IndexedReference referenceFa
        String outTablePath = "~{pairId}.stats.csv"
        Int diskSize = 10
        Int memoryGb = 20

    }
    
    command {
        set -e -o pipefail
        
        export HGREF=~{referenceFa.fasta}
        
        som.py \
        --reference ~{referenceFa.fasta} \
        --normalize-all \
        ~{oldVcf} \
        ~{newVcf}  \
        --scratch-prefix comparison \
        -o "~{pairId}"
    }
    
    output {
        File outTable = "~{outTablePath}"
        File oldOnlyVcf = "comparison/tpfn/0000.vcf.gz"
        File newOnlyVcf = "comparison/tpfn/0001.vcf.gz"
        File concordantVcf = "comparison/tpfn/0002.vcf.gz"
        File concordantFromNewVcf = "comparison/tpfn/0003.vcf.gz" # use new records for concordant annotation
    }
    
    runtime {
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-public/hap:v0.3.9"
    }
    
    meta {
        ouputColumns : ["None", "type", "total.truth", "total.query", 
            "tp", "fp", "fn", "unk", "ambi", "recall", "recall_lower", 
            "recall_upper", "recall2", "precision", "precision_lower", 
            "precision_upper", "na", "ambiguous", "fp.region.size", 
            "fp.rate", "sompyversion", "sompycmd"]
    }
}

task CompareBedPe {
    input {
        String pairId
        File oldBedpe
        File newBedpe
        String outFileSummaryPath = "~{pairId}.bedpe.concordance.tsv"
        String outFileBedpePath = "~{pairId}.bedpe.concordance.bedpe"
        Int slop = 300
        Float sizeMargin = 0.8
        Int diskSize = 10
        Int memoryGb = 20
        File bedpeConcordance = "gs://nygc-comp-s-fd4e-input/bedpe-concordance.r"
    }
    
    command {
        Rscript \
        ~{bedpeConcordance} \
        --old_bedpe=~{oldBedpe} \
        --new_bedpe=~{newBedpe} \
        --slop=~{slop} \
        --size_margin=~{sizeMargin} \
        --out_file_bedpe=~{outFileBedpePath} \
        --out_file_summary=~{outFileSummaryPath}
    }
    
    output {
        File outFileSummary = "~{outFileSummaryPath}"
        File outFileBedpe = "~{outFileBedpePath}"
    }
    
    runtime {
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-internal-tools/sv_cnv:1.0.0"
    }
}

task CompareBed {
    input {
        String pairId
        File oldBed
        File newBed
        String outFileSummaryPath = "~{pairId}.bed.concordance.tsv"
        String outFileBedPath = "~{pairId}.bed.concordance.bed"
        Float overlapFraction = 0.8
        Int diskSize = 10
        Int memoryGb = 20
        File bedConcordance = "gs://nygc-comp-s-fd4e-input/bed-concordance.r"
    }
    
    command {
        Rscript \
        ~{bedConcordance} \
        --old_bed=~{oldBed} \
        --new_bed=~{newBed} \
        --overlap_fraction=~{overlapFraction} \
        --out_file_bed=~{outFileBedPath} \
        --out_file_summary=~{outFileSummaryPath}
    }
    
    output {
        File outFileSummary = "~{outFileSummaryPath}"
        File outFileBed = "~{outFileBedPath}"
    }
    
    runtime {
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-internal-tools/sv_cnv:1.0.0"
    }
    
    meta {
        summaryTsvColumns : ["callset", "count"]
    }
}