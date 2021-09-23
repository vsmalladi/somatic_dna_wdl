version 1.0

import "../wdl_structs.wdl"

task CompareCnvGenes {
    input {
        String pairId
        String name
        File concordanceBed
        String cnvGenesTablePath = "~{pairId}.~{name}.cnv.gene.concordance.csv"
        Int diskSize = 1
        Int memoryGb = 1
        File compareBed = "gs://nygc-comp-s-fd4e-input/compare_bed.py"
    }
    
    command {
        python \
        ~{compareBed} \
        --bed-file ~{concordanceBed} \
        --output ~{cnvGenesTablePath}
    }
    
    output {
        File cnvGenesTable = "~{cnvGenesTablePath}"
    }
    
    runtime {
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-internal-tools/somatic_tools@sha256:46ab81b8dc09d6f8cf90c81f7d5692f23d73c134df6dbcd5298abde7f414dce3"
    }
    
    meta {
        summaryTsvColumns : ["pair_id", "gene", "tier", "type", "coordinates", "size_kb", "log2,callset"]
        internalOnly : "True, produces empty output for external BED files"
    }
}

task CompareSvGenes {
    input {
        String pairId
        String name
        File concordanceBedPe
        String svGenesTablePath = "~{pairId}.~{name}.sv.gene.concordance.csv"
        Int diskSize = 1
        Int memoryGb = 1
        File compareBedPe = "gs://nygc-comp-s-fd4e-input/compare_bedpe.py"
    }
    
    command {
        python \
        ~{compareBedPe} \
        --bedpe-file ~{concordanceBedPe} \
        --output ~{svGenesTablePath}
    }
    
    output {
        File svGenesTable = "~{svGenesTablePath}"
    }
    
    runtime {
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-internal-tools/somatic_tools@sha256:46ab81b8dc09d6f8cf90c81f7d5692f23d73c134df6dbcd5298abde7f414dce3"
    }
    
    meta {
        summaryTsvColumns : ["pair_id", "gene", "tier", "type", "coordinates", "size_kb", "log2,callset"]
        internalOnly : "True, produces empty output for external BEDPE files"
    }
}

task ConcateTables {
    input {
        String outputTablePath
        Array[File] tables
        Int diskSize = 1
        Int memoryGb = 1
    }
    
    command {
       python \
        /concate_tables.py \
        --tables ~{sep=" " tables} \
        --output ~{outputTablePath}
    }
    
    output {
        File outputTable = "~{outputTablePath}"
    }
    
    runtime {
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-internal-tools/somatic_tools@sha256:46ab81b8dc09d6f8cf90c81f7d5692f23d73c134df6dbcd5298abde7f414dce3"
    }
    
    meta {
        internalOnly : "False, produces output for external or internal CSV files"
    }
}

task SummarizeFlagstat {
    input {
        String sampleId
        File oldflagStat
        File newflagStat
        String outputTablePath = "~{sampleId}.flagstat.summary.csv"
        Int diskSize = 1
        Int memoryGb = 1
        File compareBams = "gs://nygc-comp-s-fd4e-input/compare_bam.py"
    }
    
    command {
       python \
        ~{compareBams} \
        --old-flagstat ~{oldflagStat} \
        --new-flagstat ~{newflagStat} \
        --sample-id ~{sampleId} \
        --output ~{outputTablePath}
    }
    
    output {
        File outputTable = "~{outputTablePath}"
    }
    
    runtime {
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-internal-tools/somatic_tools@sha256:46ab81b8dc09d6f8cf90c81f7d5692f23d73c134df6dbcd5298abde7f414dce3"
    }
    
    meta {
        outputTableColumns : ["sample_id", "version", "duplicates", "mate_diff_chr", "not_paired", "singletons", "unmapped"]
        internalOnly : "False, produces output for external or internal Flagstat files"
    }
}

task SummarizeMsi {
    input {
        String sampleId
        File oldmantisStatusFinal
        File newmantisStatusFinal
        String outputTablePath = "~{sampleId}.msi.summary.csv"
        Int diskSize = 1
        Int memoryGb = 1
        File compareMsis = "gs://nygc-comp-s-fd4e-input/compare_msi.py"
    }
    
    command {
       python \
        ~{compareMsis} \
        --old-mantis-final ~{oldmantisStatusFinal} \
        --new-mantis-final ~{newmantisStatusFinal} \
        --sample-id ~{sampleId} \
        --output ~{outputTablePath}
    }
    
    output {
        File outputTable = "~{outputTablePath}"
    }
    
    runtime {
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-internal-tools/somatic_tools@sha256:46ab81b8dc09d6f8cf90c81f7d5692f23d73c134df6dbcd5298abde7f414dce3"
    }
    
    meta {
        outputTableColumns : ["sample_id", "version", "MSI_Status", "MSI_Score", "Threshold"]
        internalOnly : "True, produces output for internal threshold adjusted MSI status files"
    }
}

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
        --pair-id ~{pairId} \
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
        docker: "gcr.io/nygc-internal-tools/somatic_tools@sha256:46ab81b8dc09d6f8cf90c81f7d5692f23d73c134df6dbcd5298abde7f414dce3"
    }
    
    meta {
        summaryTsvColumns : ["Callset", "Confidence", "Count"]
        Callsets : ["Old unique", "New unique", "Concordant"]
        internalOnly : "False, produces output for external or internal VCF files"
    }
}

task DraftComparison {
    input {
        String name
        File chromLengths
        File cosmicCensus
        Array[String] listOfChroms
        Array[File] cnvBeds
        Array[File] svBedPes
        Array[File] vcfDetails
        Array[File] oldBafs
        Array[File] newBafs
        File bedCounts
        File bedPeCounts
        File bedPeGenes
        File bedGenes
        File vcfCounts
        File flagstatSummary
        File msiSummary
        
        String mdPath = "~{name}.v7.comparison.report.md"
        String headerPath = "~{name}.v7.comparison.report.header.txt"
        
        Int diskSize = 10
        Int memoryGb = 20
        File draftComparison = "gs://nygc-comp-s-fd4e-input/draft_comparison.py"
        File plotComparison = "gs://nygc-comp-s-fd4e-input/plot_comparison.py"
        File prep = "gs://nygc-comp-s-fd4e-input/prep.py"
        File compose = "gs://nygc-comp-s-fd4e-input/compose.py"
        File ploter = "gs://nygc-comp-s-fd4e-input/ploter.py"
        File colorer = "gs://nygc-comp-s-fd4e-input/Colorer.py" 
    }
    
    command {
       python \
        ~{draftComparison} \
        --chrom-lengths ~{chromLengths} \
        --cosmic-census ~{cosmicCensus} \
        --chroms ~{sep=" " listOfChroms} \
        --name ~{name} \
        --cnv-summary ~{bedCounts} \
        --cnv-beds ~{sep=" " cnvBeds}  \
        --sv-summary ~{bedPeCounts} \
        --sv-bedpes ~{sep=" " svBedPes} \
        --sv-table ~{bedPeGenes} \
        --cnv-table ~{bedGenes} \
        --vcf-gene-tables ~{sep=" " vcfDetails} \
        --snv-summary ~{vcfCounts} \
        --old-bafs ~{sep=" " oldBafs} \
        --new-bafs ~{sep=" " newBafs} \
        --flagstat-summary ~{flagstatSummary} \
        --msi-summary ~{msiSummary}
    }
    
    output {
        File md = "~{mdPath}"
        File header = "~{headerPath}"
    }
    
    runtime {
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-public/bokeh_plotly:1.0.1"
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
        internalOnly : "False, produces output for external or internal VCF files"
    }
}

task CompareBedPe {
    input {
        String pairId
        File oldBedpe
        File newBedpe
        String outFileSummaryPath = "~{pairId}.sv.concordance.tsv"
        String outFileBedpePath = "~{pairId}.sv.concordance.bedpe"
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
        --out_file_summary=~{outFileSummaryPath} \
        --pair_id=~{pairId}
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
    
    meta {
        summaryTsvColumns : ["callset", "sv_type", "count"]
        internalOnly : "False, produces output for external or internal BEDPE files"
    }
}

task CompareBed {
    input {
        String pairId
        File oldBed
        File newBed
        String outFileSummaryPath = "~{pairId}.cnv.concordance.tsv"
        String outFileBedPath = "~{pairId}.cnv.concordance.bed"
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
        --out_file_summary=~{outFileSummaryPath} \
        --pair_id=~{pairId}
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
        internalOnly : "False, produces output for external or internal BED files"
    }
}