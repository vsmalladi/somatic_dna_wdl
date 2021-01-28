version 1.0

import "../wdl_structs.wdl"

task MultipleMetrics {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String MultipleMetricsBase
        IndexedReference referenceFa
        Bam finalBam
        String sampleId
    }

    command {
        gatk \
        CollectMultipleMetrics \
        --java-options "-Xmx24576m -XX:ParallelGCThreads=1" \
        --PROGRAM CollectAlignmentSummaryMetrics \
        --PROGRAM CollectInsertSizeMetrics \
        --PROGRAM QualityScoreDistribution \
        --PROGRAM MeanQualityByCycle \
        --REFERENCE_SEQUENCE ~{referenceFa.fasta} \
        --INCLUDE_UNPAIRED true \
        --VALIDATION_STRINGENCY SILENT \
        -O ~{MultipleMetricsBase} \
        -I ~{finalBam.bam}
    }

    output {
        File alignmentSummaryMetrics = "~{sampleId}.MultipleMetrics.alignment_summary_metrics"
        File qualityByCyclePdf = "~{sampleId}.MultipleMetrics.quality_by_cycle.pdf"
        File baseDistributionByCycleMetrics = "~{sampleId}.MultipleMetrics.base_distribution_by_cycle_metrics"
        File qualityByCycleMetrics = "~{sampleId}.MultipleMetrics.quality_by_cycle_metrics"
        File baseDistributionByCyclePdf = "~{sampleId}.MultipleMetrics.base_distribution_by_cycle.pdf"
        File qualityDistributionPdf = "~{sampleId}.MultipleMetrics.quality_distribution.pdf"
        File qualityDistributionMetrics = "~{sampleId}.MultipleMetrics.quality_distribution_metrics"
        File insertSizeHistogramPdf = "~{sampleId}.MultipleMetrics.insert_size_histogram.pdf"
        File insertSizeMetrics = "~{sampleId}.MultipleMetrics.insert_size_metrics"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task MultipleMetricsPreBqsr {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String MultipleMetricsBasePreBqsrBasename
        IndexedReference referenceFa
        Bam mergedDedupBam
        String sampleId
    }

    command {
        gatk \
        CollectMultipleMetrics \
        --java-options "-Xmx24576m -XX:ParallelGCThreads=1" \
        --PROGRAM QualityScoreDistribution \
        --PROGRAM MeanQualityByCycle \
        --PROGRAM CollectGcBiasMetrics \
        --VALIDATION_STRINGENCY SILENT \
        --INCLUDE_UNPAIRED true \
        --REFERENCE_SEQUENCE ~{referenceFa.fasta} \
        -O ~{MultipleMetricsBasePreBqsrBasename} \
        -I ~{mergedDedupBam.bam}
    }

    output {
        File qualityDistributionPdfPreBqsr = "~{sampleId}.dedup.MultipleMetrics.quality_distribution.pdf"
        File qualityByCycleMetricsPreBqsr = "~{sampleId}.dedup.MultipleMetrics.quality_by_cycle_metrics"
        File qualityByCyclePdfPreBqsr = "~{sampleId}.dedup.MultipleMetrics.quality_by_cycle.pdf"
        File qualityDistributionMetricsPreBqsr = "~{sampleId}.dedup.MultipleMetrics.quality_distribution_metrics"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task CollectGcBiasMetrics {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String sampleId
        String gcBiasPdfPath = "~{sampleId}.GcBiasMetrics.gc_bias.pdf"
        String gcBiasMetricsPath = "~{sampleId}.GcBiasMetrics.gc_bias_metrics"
        String gcBiasSummaryPath = "~{sampleId}.GcBiasMetrics.gc_bias_summary"
        IndexedReference referenceFa
        Bam finalBam
    }

    command {
        gatk \
        CollectGcBiasMetrics \
        --java-options "-Xmx24576m -XX:ParallelGCThreads=1" \
        --CHART_OUTPUT ~{gcBiasPdfPath} \
        -O ~{gcBiasMetricsPath} \
        -I ~{finalBam.bam} \
        --SUMMARY_OUTPUT ~{gcBiasSummaryPath} \
        --REFERENCE_SEQUENCE ~{referenceFa.fasta} \
        --ASSUME_SORTED true \
        --VALIDATION_STRINGENCY SILENT
    }

    output {
        File gcBiasMetrics = "~{gcBiasMetricsPath}"
        File gcBiasSummary = "~{gcBiasSummaryPath}"
        File gcBiasPdf = "~{gcBiasPdfPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task Flagstat {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String sampleId
        String FlagStatPath = "~{sampleId}.FlagStat.txt"
        IndexedReference referenceFa
        Bam finalBam
    }

    command {
        gatk \
        FlagStat \
        --java-options "-Xmx24576m" \
        --verbosity INFO \
        --reference ~{referenceFa.fasta} \
        -I ~{finalBam.bam} \
        > ~{FlagStatPath}
    }

    output {
        File FlagStat = "~{FlagStatPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task HsMetrics {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String sampleId
        String HsMetricsPath = "~{sampleId}.HsMetrics.txt"
        String HsMetricsPerTargetCoveragePath = "~{sampleId}.HsMetrics.perTargetCoverage.txt"
        IndexedReference referenceFa
        Bam finalBam
        File hsMetricsIntervals
    }

    command {
        gatk \
        CollectHsMetrics \
        --java-options "-Xmx24576m -XX:ParallelGCThreads=1" \
        --BAIT_INTERVALS ~{hsMetricsIntervals} \
        --TARGET_INTERVALS ~{hsMetricsIntervals} \
        --BAIT_SET_NAME ~{sampleId} \
        --METRIC_ACCUMULATION_LEVEL ALL_READS \
        --MINIMUM_MAPPING_QUALITY 1 \
        --MINIMUM_BASE_QUALITY 0 \
        --CLIP_OVERLAPPING_READS false \
        --REFERENCE_SEQUENCE ~{referenceFa.fasta} \
        -O ~{HsMetricsPath} \
        -I ~{finalBam.bam} \
        --PER_TARGET_COVERAGE ~{HsMetricsPerTargetCoveragePath} \
        --VALIDATION_STRINGENCY SILENT
    }

    output {
        File HsMetrics = "~{HsMetricsPath}"
        File HsMetricsPerTargetCoverage = "~{HsMetricsPerTargetCoveragePath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task FormatHsMetrics {
    input {
        Int threads
        Int memory_gb
        String dockerImage
        String sampleId
        String HsMetricsPerTargetCoverageAutocorrPath = "~{sampleId}.HsMetrics.perTargetCoverage.txt.autocorr"
        File HsMetricsPerTargetCoverage
    }

    command {
        create_autocorrelation_input.v.0.1.pl \
        -input ~{HsMetricsPerTargetCoverage} \
        > ~{HsMetricsPerTargetCoverageAutocorrPath} \
    }

    output {
        File HsMetricsPerTargetCoverageAutocorr = "~{HsMetricsPerTargetCoverageAutocorrPath}"
    }

    runtime {
        cpu : threads
        memory : memory_gb + "GB"
        docker : dockerImage
    }
}

task Autocorrelations {
    input {
        Int threads
        Int memory_gb
        String dockerImage
        String sampleId
        File HsMetricsPerTargetCoverageAutocorr
    }

    command {
        R --no-save \
        --args \
        "./" \
        ~{HsMetricsPerTargetCoverageAutocorr} \
        ~{sampleId} \
        < ASP_modified_final.v.0.1.R \
    }

    output {
        File autocorroutput1100 = "~{sampleId}.autocorroutput.1.100.txt"
    }

    runtime {
        cpu : threads
        memory : memory_gb + "GB"
        docker : dockerImage
    }
}


task CollectOxoGMetricsWgs {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String sampleId
        String CollectOxoGMetricsPath = "~{sampleId}.CollectOxoGMetrics.txt"
        IndexedReference referenceFa
        Bam finalBam
    }

    command {
        gatk \
        CollectOxoGMetrics \
        --java-options "-Xmx24576m -XX:ParallelGCThreads=1" \
        --VALIDATION_STRINGENCY SILENT \
        -I ~{finalBam.bam} \
        -O ~{CollectOxoGMetricsPath} \
        --REFERENCE_SEQUENCE ~{referenceFa.fasta}
    }

    output {
        File CollectOxoGMetrics = "~{CollectOxoGMetricsPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task CollectWgsMetricsWgsDecoy {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String sampleId
        String CollectWgsMetricsPath = "~{sampleId}.CollectWgsMetrics.txt"
        IndexedReference referenceFa
        Bam finalBam
        File randomIntervals
    }

    command {
        gatk \
        CollectWgsMetrics \
        --java-options "-Xmx24576m -XX:ParallelGCThreads=1" \
        --VALIDATION_STRINGENCY SILENT \
        -I ~{finalBam.bam} \
        -O ~{CollectWgsMetricsPath} \
        --INTERVALS ~{randomIntervals} \
        -R ~{referenceFa.fasta} \
        --MINIMUM_MAPPING_QUALITY 0 \
        --COVERAGE_CAP 1000 \
        --MAX_RECORDS_IN_RAM 2000000 \
        --COUNT_UNPAIRED true \
        --MINIMUM_BASE_QUALITY 3
    }

    output {
        File CollectWgsMetrics = "~{CollectWgsMetricsPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task Binest {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String sampleId
        String binestCovPath = "~{sampleId}.binest.coverage.txt"
        File finalBam
    }

    command {
        binest \
        size \
        ~{finalBam.bamIndex} \
        > ~{binestCovPath}
    }

    output {
        File binestCov = "~{binestCovPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task PlotBinCov {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String genome
        File genomeTemplates
        String sampleId
        File binestCov
    }

    command {
        plot_bin_cov.R \
        "--binestOutput=~{binestCov}" \
        "--genome=~{genome}" \
        "--genomeTemplates=~{genomeTemplates}" \
        "--sample=~{sampleId}"
    }

    output {
        File normCoverageByChrPng = "~{sampleId}.binest.coverage.png"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task Pileup {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String sampleId
        String pileupsTablePath = "~{sampleId}_pileups_table.table"
        Bam finalBam
        File gnomadBiallelic
    }

    command {
        gatk \
        --java-options "-Xmx30g" \
        GetPileupSummaries \
        -I ~{finalBam.bam} \
        -V ~{gnomadBiallelic} \
        -O ~{pileupsTablePath}
    }

    output {
        File pileupsTable = "~{pileupsTablePath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task CalculateContamination {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String sampleId
        String contaminationTablePath = "~{sampleId}.contamination.table"
        File pileupsTable
    }

    command {
        gatk \
        --java-options "-Xmx30g" \
        CalculateContamination \
        -I ~{pileupsTable} \
        -O ~{contaminationTablePath}
    }

    output {
        File contaminationTable = "~{contaminationTablePath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task ConpairPileup {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        IndexedReference referenceFa
        String sampleId
        String pileupsConpairPath = "~{sampleId}_pileups_table.txt"
        Bam finalBam
        File markerFile
    }

    command {
        java \
        "-Xmx12g" \
        -jar GenomeAnalysisTK.jar \
        -T Pileup \
        -R ~{referenceFa.fasta} \
        -I ~{finalBam.bam} \
        -L ~{markerFile} \
        -o ~{pileupsConpairPath} \
        -verbose \
        -rf DuplicateRead \
        --filter_reads_with_N_cigar \
        --filter_mismatching_base_and_quals
    }

    output {
        File pileupsConpair = "~{pileupsConpairPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}









