version 1.0

import "../wdl_structs.wdl"

task MultipleMetrics {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        IndexedReference referenceFa
        Bam finalBam
        String sampleId
        String MultipleMetricsBase = "~{sampleId}.MultipleMetrics"
    }

    command {
        gatk \
        CollectMultipleMetrics \
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
        File alignmentSummaryMetrics = "~{MultipleMetricsBase}.alignment_summary_metrics"
        File qualityByCyclePdf = "~{MultipleMetricsBase}.quality_by_cycle.pdf"
        File baseDistributionByCycleMetrics = "~{MultipleMetricsBase}.base_distribution_by_cycle_metrics"
        File qualityByCycleMetrics = "~{MultipleMetricsBase}.quality_by_cycle_metrics"
        File baseDistributionByCyclePdf = "~{MultipleMetricsBase}.base_distribution_by_cycle.pdf"
        File qualityDistributionPdf = "~{MultipleMetricsBase}.quality_distribution.pdf"
        File qualityDistributionMetrics = "~{MultipleMetricsBase}.quality_distribution_metrics"
        File insertSizeHistogramPdf = "~{MultipleMetricsBase}.insert_size_histogram.pdf"
        File insertSizeMetrics = "~{MultipleMetricsBase}.insert_size_metrics"
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
        --java-options "-XX:ParallelGCThreads=1" \
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
        File qualityDistributionPdfPreBqsr = "~{MultipleMetricsBasePreBqsrBasename}.quality_distribution.pdf"
        File qualityByCycleMetricsPreBqsr = "~{MultipleMetricsBasePreBqsrBasename}.quality_by_cycle_metrics"
        File qualityByCyclePdfPreBqsr = "~{MultipleMetricsBasePreBqsrBasename}.quality_by_cycle.pdf"
        File qualityDistributionMetricsPreBqsr = "~{MultipleMetricsBasePreBqsrBasename}.quality_distribution_metrics"
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
        --java-options "-XX:ParallelGCThreads=1" \
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
        Int memoryGb
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
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task Autocorrelations {
    input {
        Int threads
        Int memoryGb
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
        memory : memoryGb + "GB"
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
        Bam finalBam
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

task CalculateContaminationPaired {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String pairName
        String contaminationTablePath = "~{pairName}.contamination.table"
        File pileupsNormalTable
        File pileupsTumorTable
    }

    command {
        gatk \
        CalculateContamination \
        -I ~{pileupsTumorTable} \
        -matched ~{pileupsNormalTable} \
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
        File markerBedFile
    }

    command {
        java \
        -jar GenomeAnalysisTK.jar \
        -T Pileup \
        -R ~{referenceFa.fasta} \
        -I ~{finalBam.bam} \
        -L ~{markerBedFile} \
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

task VerifyConcordanceAll {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        File pileupsTumorConpair
        File pileupsNormalConpair
        File markerTxtFile
        String pairName
        String concordanceAllPath = "~{pairName}.concordance.all.conpair-0.1.txt"
    }

    command {
        verify_concordance.py \
        -T ~{pileupsTumorConpair} \
        -N ~{pileupsNormalConpair} \
        -O ~{concordanceAllPath} \
        -D "./" \
        -M ~{markerTxtFile}
    }

    output {
        File concordanceAll = "~{concordanceAllPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task VerifyConcordanceHomoz {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        File pileupsTumorConpair
        File pileupsNormalConpair
        File markerTxtFile
        String pairName
        String concordanceHomozPath = "~{pairName}.concordance.homoz.conpair-0.1.txt"
    }

    command {
        verify_concordance.py \
        -T ~{pileupsTumorConpair} \
        -N ~{pileupsNormalConpair} \
        -O ~{concordanceHomozPath} \
        -D "./" \
        -M ~{markerTxtFile} \
        -H
    }

    output {
        File concordanceHomoz = "~{concordanceHomozPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task Contamination {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        File pileupsTumorConpair
        File pileupsNormalConpair
        File markerTxtFile
        String pairName
        String contaminationPath = "~{pairName}.contamination.conpair-0.1.txt"
    }

    command {
        verify_concordance.py \
        -T ~{pileupsTumorConpair} \
        -N ~{pileupsNormalConpair} \
        -O ~{contaminationPath} \
        -P 0.001 \
        -Q 10 \
        -D "./" \
        -M ~{markerTxtFile}
    }

    output {
        File contamination = "~{contaminationPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}
