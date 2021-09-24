version 1.0

import "../wdl_structs.wdl"

task FilterHighConfidence {
    input {            
        File vcf
        String filteredVcfPath
        Int diskSize = 1
        Int memoryGb = 1
    }
    
    command {
        bash \
        /filter_confidence.sh \
        ~{vcf} \
        ~{filteredVcfPath}
    }
    
    output {
        File filteredVcf = "~{filteredVcfPath}"
    }
    
    runtime {
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-internal-tools/somatic_reports@sha256:ba18bd0fd7ce55af26e4d255f85529a728951bab02596e0528391b2b52d05045"
    }
    
    meta {
        summaryTsvColumns : ["mean_cov", "sample_id"]
        internalOnly : "False, but depends on filename to add sample_id"
    }
}

task SummarizeMantis {
    input {            
        String pipeline = "v7"
        String name
        Array[File] mantisStatusFinal
        String mantisStatusTablePath = "~{name}.~{pipeline}.msi.csv"
        Int diskSize = 1
        Int memoryGb = 1
    }
    
    command {
        python \
        /prep_Mantis.py \
        --output ~{mantisStatusTablePath} \
        --mantis-files ~{sep=" " mantisStatusFinal}
    }
    
    output {
        File mantisStatusTable = "~{mantisStatusTablePath}"
    }
    
    runtime {
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-internal-tools/somatic_reports@sha256:ba18bd0fd7ce55af26e4d255f85529a728951bab02596e0528391b2b52d05045"
    }
    
    meta {
        summaryTsvColumns : ["mean_cov", "sample_id"]
        internalOnly : "False, but depends on filename to add sample_id"
    }
}

task SummarizeQualityByCycle {
    input {            
        String pipeline = "v7"
        String name
        Array[File] qualityByCycleMetrics
        String qualityByCycleTablePath = "~{name}.~{pipeline}.quality_by_cycle_metrics.csv"
        Int diskSize = 1
        Int memoryGb = 1
    }
    
    command {
        python \
        /prep_QualityByCycle.py \
        --output ~{qualityByCycleTablePath} \
        --quality-files ~{sep=" " qualityByCycleMetrics}
    }
    
    output {
        File qualityByCycleTable = "~{qualityByCycleTablePath}"
    }
    
    runtime {
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-internal-tools/somatic_reports@sha256:ba18bd0fd7ce55af26e4d255f85529a728951bab02596e0528391b2b52d05045"
    }
    
    meta {
        summaryTsvColumns : ["mean_cov", "sample_id"]
        internalOnly : "False, but depends on filename to add sample_id"
    }
}



task SummarizeInsertSizeMetrics {
    input {            
        String pipeline = "v7"
        String name
        Array[File] insertSizeMetrics
        String insertSizeTablePath = "~{name}.~{pipeline}.insert_size.csv"
        Int diskSize = 1
        Int memoryGb = 1
    }
    
    command {
        python \
        /prep_InsertSizeMetrics.py \
        --output ~{insertSizeTablePath} \
        --insert-sizes ~{sep=" " insertSizeMetrics}
    }
    
    output {
        File insertSizeTable = "~{insertSizeTablePath}"
    }
    
    runtime {
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-internal-tools/somatic_reports@sha256:ba18bd0fd7ce55af26e4d255f85529a728951bab02596e0528391b2b52d05045"
    }
    
    meta {
        summaryTsvColumns : ["mean_cov", "sample_id"]
        internalOnly : "False, but depends on filename to add sample_id"
    }
}


task SummarizeCollectWgsMetrics {
    input {            
        String pipeline = "v7"
        String name
        Array[File] collectWgsMetrics
        String coverageTablePath = "~{name}.~{pipeline}.cov.csv"
        Int diskSize = 1
        Int memoryGb = 1
    }
    
    command {
        python \
        /prep_CollectWgsMetrics.py \
        --output ~{coverageTablePath} \
        --coverage-files ~{sep=" " collectWgsMetrics}
    }
    
    output {
        File coverageTable = "~{coverageTablePath}"
    }
    
    runtime {
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-internal-tools/somatic_reports@sha256:ba18bd0fd7ce55af26e4d255f85529a728951bab02596e0528391b2b52d05045"
    }
    
    meta {
        summaryTsvColumns : ["mean_cov", "sample_id"]
        internalOnly : "False, but depends on filename to add sample_id"
    }
}

task SummarizeSvs {
    input {            
        String pipeline = "v7"
        String name
        Array[File] highConfidenceSvTables
        Array[File] allSomaticSvTables
        String svTablePath = "~{name}.~{pipeline}.summarySv.csv"
        Int diskSize = 1
        Int memoryGb = 1
    }
    
    command {
        python \
        /summarize_svs.py \
        --output ~{svTablePath} \
        --high-conf-sv-tables ~{sep=" " highConfidenceSvTables} \
        --all-somatic-sv-tables ~{sep=" " allSomaticSvTables}
    }
    
    output {
        File summarySvTable = "~{svTablePath}"
    }
    
    runtime {
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-internal-tools/somatic_reports@sha256:ba18bd0fd7ce55af26e4d255f85529a728951bab02596e0528391b2b52d05045"
    }
    
    meta {
        summaryTsvColumns : ["Confidence", "Count", "pair_id"]
        internalOnly : "False, but depends on filename to add sample_id"
    }
}

task SummarizeFlagStat {
    input {            
        String pipeline = "v7"
        String name
        Array[File] flagStats
        String flagStatTablePath = "~{name}.~{pipeline}.flagstat.csv"
        Int diskSize = 1
        Int memoryGb = 1
    }
    
    command {
        python \
        /prep_flagstat.py \
        --output ~{flagStatTablePath} \
        --flagstat-files ~{sep=" " flagStats}
    }
    
    output {
        File flagStatTable = "~{flagStatTablePath}"
    }
    
    runtime {
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-internal-tools/somatic_reports@sha256:ba18bd0fd7ce55af26e4d255f85529a728951bab02596e0528391b2b52d05045"
    }
    
    meta {
        summaryTsvColumns : ["flagstat", "value", "sample_id"]
        internalOnly : "False, but depends on filename to add sample_id"
    }
}

task DraftSampleReport {
    input {
        String pairId
        Array[String] listOfChroms
        File chromLengths
        Array[File] cnvTable
        Array[File] cnvGenesTable
        
        Array[File] svGenesTable
        Array[File] svTable
        
        Array[File] detailedVcfTable
        Array[File] summaryVcfTable
        
        String mdPath = "~{pairId}.v7.final.report.md"
        String headerPath = "~{pairId}.v7.final.report.header.txt"
        
        Int diskSize = 10
        Int memoryGb = 20
    }
    
    command {
       python \
        /draft_sample_report.py \
        --chrom-lengths ~{chromLengths} \
        --chroms ~{sep=" " listOfChroms} \
        --pair-id ~{pairId} \
        --cnv-gene-tables ~{sep=" " cnvTable} \
        --cnv-tables ~{sep=" " cnvGenesTable} \
        --sv-gene-tables ~{sep=" " svGenesTable} \
        --sv-tables ~{sep=" " svTable} \
        --snv-gene-tables ~{sep=" " detailedVcfTable} \
        --snv-summary-tables ~{sep=" " summaryVcfTable} \
    }
    
    output {
        File md = "~{mdPath}"
        File header = "~{headerPath}"
    }
    
    runtime {
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-internal-tools/somatic_reports@sha256:ba18bd0fd7ce55af26e4d255f85529a728951bab02596e0528391b2b52d05045"
    }

}

task DescribeBaf {
    input {            
        String pipeline = "v7"
        String pairId
        File bedpe
        File chromLengths
        Array[String] listOfChroms
        String svTablePath = "~{pairId}.~{pipeline}.baf.csv"
        Int diskSize = 1
        Int memoryGb = 1
    }
    
    command {
        python \
        /prep_fusions.py \
        --sv-bedpes ~{bedpe} \
        --output ~{svTablePath} \
        --chrom-lengths ~{chromLengths} \
        --chroms ~{sep=" " listOfChroms} \
        --pair-id ~{pairId}
    }
    
    output {
        File svTable = "~{svTablePath}"
    }
    
    runtime {
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-internal-tools/somatic_reports@sha256:ba18bd0fd7ce55af26e4d255f85529a728951bab02596e0528391b2b52d05045"
    }
    
    meta {
        summaryTsvColumns : ["pair_id", "type", "coordinates", "log2", "chrom", "start", "end", "cumulative_start_pos", "cumulative_end_pos", "xs", "kind"]
        internalOnly : "True, produces empty output for external BED files"
    }
}

task DescribeBedPe {
    input {            
        String pipeline = "v7"
        String name
        String pairId
        File bedpe
        File chromLengths
        Array[String] listOfChroms
        String svTablePath = "~{pairId}.~{name}.~{pipeline}.fusions.sv.csv"
        Int diskSize = 1
        Int memoryGb = 1
    }
    
    command {
        python \
        /prep_fusions.py \
        --sv-bedpes ~{bedpe} \
        --output ~{svTablePath} \
        --chrom-lengths ~{chromLengths} \
        --chroms ~{sep=" " listOfChroms} \
        --pair-id ~{pairId}
    }
    
    output {
        File svTable = "~{svTablePath}"
    }
    
    runtime {
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-internal-tools/somatic_reports@sha256:ba18bd0fd7ce55af26e4d255f85529a728951bab02596e0528391b2b52d05045"
    }
    
    meta {
        summaryTsvColumns : ["pair_id", "type", "coordinates", "log2", "chrom", "start", "end", "cumulative_start_pos", "cumulative_end_pos", "xs", "kind"]
        internalOnly : "True, produces empty output for external BED files"
    }
}

task DescribeBedPeGenes {
    input {            
        String pipeline = "v7"
        String pairId
        File bedpe
        String svGenesTablePath = "~{pairId}.~{pipeline}.genes.sv.csv"
        Int diskSize = 1
        Int memoryGb = 1
    }
    
    command {
        python \
        /describe_bedpe_genes.py \
        --bedpe-file ~{bedpe} \
        --output ~{svGenesTablePath} \
        --pair-id ~{pairId}
    }
    
    output {
        File svGenesTable = "~{svGenesTablePath}"
    }
    
    runtime {
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-internal-tools/somatic_reports@sha256:ba18bd0fd7ce55af26e4d255f85529a728951bab02596e0528391b2b52d05045"
    }
    
    meta {
        summaryTsvColumns : ["pair_id", "type", "coordinates", "log2", "chrom", "start", "end", "cumulative_start_pos", "cumulative_end_pos", "xs", "kind"]
        internalOnly : "True, produces empty output for external BED files"
    }
}

task SummarizeFinalVcf {
    input {
        String pairId
        String pipeline = "v7"
        File vcf
        String detailedOutputTablePath = "~{pairId}.~{pipeline}.detailed.output.csv"
        String summaryOutputTablePath = "~{pairId}.~{pipeline}.counts.csv"
        
        Int diskSize = 10
        Int memoryGb = 20
    }
    
    command {
       python \
        /summarize_vcfs.py \
        --pair-id ~{pairId} \
        --vcf ~{vcf} \
        --output ~{detailedOutputTablePath} \
        --summary ~{summaryOutputTablePath}
    }
    
    output {
        File detailedVcfTable = "~{detailedOutputTablePath}"
        File summaryVcfTable = "~{summaryOutputTablePath}"
    }
    
    runtime {
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-internal-tools/somatic_reports@sha256:ba18bd0fd7ce55af26e4d255f85529a728951bab02596e0528391b2b52d05045"
    }
    
    meta {
        summaryTsvColumns : ["Callset", "Confidence", "Count"]
        Callsets : ["Old unique", "New unique", "Concordant"]
        internalOnly : "False, produces output for external or internal VCF files"
    }
}

task DescribeBed {
    input {
        String pipeline = "v7"
        String pairId
        File bed
        File chromLengths
        Array[String] listOfChroms
        String cnvTablePath = "~{pairId}.~{pipeline}.cnv.csv"
        Int diskSize = 1
        Int memoryGb = 1
    }
    
    command {
        python \
        /describe_bed.py \
        --bed-file ~{bed} \
        --output ~{cnvTablePath} \
        --chrom-lengths ~{chromLengths} \
        --chroms ~{sep=" " listOfChroms} \
        --pair-id ~{pairId}
    }
    
    output {
        File cnvTable = "~{cnvTablePath}"
    }
    
    runtime {
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-internal-tools/somatic_reports@sha256:ba18bd0fd7ce55af26e4d255f85529a728951bab02596e0528391b2b52d05045"
    }
    
    meta {
        #summaryTsvColumns : ["pair_id", "type", "coordinates", "log2", "chrom", "start", "end", "cumulative_start_pos", "cumulative_end_pos", "xs", "kind"]
        internalOnly : "True, produces empty output for external BED files"
    }
}

task DescribeBedGenes {
    input {            
        String pipeline = "v7"
        String pairId
        File bed
        String cnvGenesTablePath = "~{pairId}.~{pipeline}.genes.cnv.csv"
        Int diskSize = 1
        Int memoryGb = 1
    }
    
    command {
        python \
        /describe_bed_genes.py \
        --bed-file ~{bed} \
        --output ~{cnvGenesTablePath} \
        --pair-id ~{pairId}
    }
    
    output {
        File cnvGenesTable = "~{cnvGenesTablePath}"
    }
    
    runtime {
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-internal-tools/somatic_reports@sha256:ba18bd0fd7ce55af26e4d255f85529a728951bab02596e0528391b2b52d05045"
    }
    
    meta {
        summaryTsvColumns : ["pair_id", "type", "coordinates", "log2", "chrom", "start", "end", "cumulative_start_pos", "cumulative_end_pos", "xs", "kind"]
        internalOnly : "True, produces empty output for external BED files"
    }
}

task SummarizeHla {
    input {
        String sampleId
        File oldkouramiResult
        File newkouramiResult
        String outputTablePath = "~{sampleId}.hla.summary.csv"
        String outputConcordancePath = "~{sampleId}.hla.concordance.csv"
        Int diskSize = 1
        Int memoryGb = 1
    }

    command {
       python \
        /compare_hla.py \
        --old-kourami ~{oldkouramiResult} \
        --new-kourami ~{newkouramiResult} \
        --sample-id ~{sampleId} \
        --output ~{outputTablePath} \
        --concordance ~{outputConcordancePath}
    }

    output {
        File outputTable = "~{outputTablePath}"
        File outputConcordance = "~{outputConcordancePath}"
    }

    runtime {
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-internal-tools/somatic_reports@sha256:ba18bd0fd7ce55af26e4d255f85529a728951bab02596e0528391b2b52d05045"
    }

    meta {
        outputTableColumns : ["sample_id", "version", "A", "B", "C", "DQA1", "DQB1", "DRB1"]
        outputConcordanceColumns : ["gene", "sample_id", "discordance_count", "concordance_count"]
    }
}

task CompareCnvGenes {
    input {
        String pairId
        String name
        File concordanceBed
        String cnvGenesTablePath = "~{pairId}.~{name}.cnv.gene.concordance.csv"
        Int diskSize = 1
        Int memoryGb = 1
    }
    
    command {
        python \
        /compare_bed.py \
        --bed-file ~{concordanceBed} \
        --output ~{cnvGenesTablePath}
    }
    
    output {
        File cnvGenesTable = "~{cnvGenesTablePath}"
    }
    
    runtime {
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-internal-tools/somatic_reports@sha256:ba18bd0fd7ce55af26e4d255f85529a728951bab02596e0528391b2b52d05045"
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
    }
    
    command {
        python \
        /compare_bedpe.py \
        --bedpe-file ~{concordanceBedPe} \
        --output ~{svGenesTablePath}
    }
    
    output {
        File svGenesTable = "~{svGenesTablePath}"
    }
    
    runtime {
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-internal-tools/somatic_reports@sha256:ba18bd0fd7ce55af26e4d255f85529a728951bab02596e0528391b2b52d05045"
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
    }
    
    command {
       python \
        /compare_bam.py \
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
        docker: "gcr.io/nygc-internal-tools/somatic_reports@sha256:ba18bd0fd7ce55af26e4d255f85529a728951bab02596e0528391b2b52d05045"
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
    }
    
    command {
       python \
        /compare_msi.py \
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
        docker: "gcr.io/nygc-internal-tools/somatic_reports@sha256:ba18bd0fd7ce55af26e4d255f85529a728951bab02596e0528391b2b52d05045"
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
    }
    
    command {
       python \
        /compare_vcfs.py \
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
        docker: "gcr.io/nygc-internal-tools/somatic_reports@sha256:ba18bd0fd7ce55af26e4d255f85529a728951bab02596e0528391b2b52d05045"
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
    }
    
    command {
       python \
        /draft_comparison.py \
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
        docker: "gcr.io/nygc-internal-tools/somatic_reports@sha256:ba18bd0fd7ce55af26e4d255f85529a728951bab02596e0528391b2b52d05045"
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
        docker: "gcr.io/nygc-public/hap@sha256:6fc9a90c944d4e8434977909d0a175bde9f0cbcdd8a8e6806a65d1aa24b4a683"
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
        docker: "gcr.io/nygc-internal-tools/somatic_reports@sha256:ba18bd0fd7ce55af26e4d255f85529a728951bab02596e0528391b2b52d05045"
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
        docker: "gcr.io/nygc-internal-tools/somatic_reports@sha256:ba18bd0fd7ce55af26e4d255f85529a728951bab02596e0528391b2b52d05045"
    }
    
    meta {
        summaryTsvColumns : ["callset", "count"]
        internalOnly : "False, produces output for external or internal BED files"
    }
}