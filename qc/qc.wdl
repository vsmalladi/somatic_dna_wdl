version 1.0

import "../wdl_structs.wdl"

task MultipleMetrics {
    input {
        Int threads = 2
        Int memoryGb = 40
        Int diskSize
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
        disks: "local-disk " + diskSize + " SSD"
        docker : "us.gcr.io/broad-gatk/gatk:4.1.0.0"
    }
}

task MultipleMetricsPreBqsr {
    input {
        Int threads = 2
        Int memoryGb = 40
        Int diskSize
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
        disks: "local-disk " + diskSize + " SSD"
        docker : "us.gcr.io/broad-gatk/gatk:4.1.0.0"
    }
}

task CollectGcBiasMetrics {
    input {
        Int threads = 2
        Int memoryGb = 32
        Int diskSize
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
        disks: "local-disk " + diskSize + " SSD"
        docker : "us.gcr.io/broad-gatk/gatk:4.1.0.0"
    }
}

task Flagstat {
    input {
        Int threads = 2
        Int memoryGb = 32
        Int diskSize
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
        disks: "local-disk " + diskSize + " SSD"
        docker : "us.gcr.io/broad-gatk/gatk:4.1.0.0"
    }
}

task HsMetrics {
    input {
        Int threads = 2
        Int memoryGb = 40
        Int diskSize
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
        disks: "local-disk " + diskSize + " SSD"
        docker : "us.gcr.io/broad-gatk/gatk:4.1.0.0"
    }
}

task FormatHsMetrics {
    input {
        Int threads = 1
        Int memoryGb = 4
        String sampleId
        String HsMetricsPerTargetCoverageAutocorrPath = "~{sampleId}.HsMetrics.perTargetCoverage.txt.autocorr"
        File HsMetricsPerTargetCoverage
    }

    command {
        perl /create_autocorrelation_input.v.0.1.pl \
        -input ~{HsMetricsPerTargetCoverage} \
        > ~{HsMetricsPerTargetCoverageAutocorrPath} \
    }

    output {
        File HsMetricsPerTargetCoverageAutocorr = "~{HsMetricsPerTargetCoverageAutocorrPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/somatic_tools"
    }
}

task Autocorrelations {
    input {
        Int threads = 1
        Int memoryGb = 4
        String sampleId
        File HsMetricsPerTargetCoverageAutocorr
    }

    command {
        R --no-save \
        --args \
        "./" \
        ~{HsMetricsPerTargetCoverageAutocorr} \
        ~{sampleId} \
        < somatic_tools/ASP_modified_final.v.0.1.R \
    }

    output {
        File autocorroutput1100 = "~{sampleId}.autocorroutput.1.100.txt"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/somatic_tools"
    }
}


task CollectOxoGMetricsWgs {
    input {
        Int threads = 2
        Int memoryGb = 8
        Int diskSize
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
        disks: "local-disk " + diskSize + " SSD"
        docker : "us.gcr.io/broad-gatk/gatk:4.1.0.0"
    }
}

task CollectWgsMetricsWgsDecoy {
    input {
        Int threads = 2
        Int memoryGb = 32
        Int diskSize
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
        disks: "local-disk " + diskSize + " SSD"
        docker : "us.gcr.io/broad-gatk/gatk:4.1.0.0"
    }
}

task Binest {
    input {
        Int threads = 1
        Int memoryGb = 4
        Int diskSize
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
        disks: "local-disk " + diskSize + " SSD"
        docker : "gcr.io/nygc-compbio/binest:0.8.4"
    }
}

task PlotBinCov {
    input {
        Int threads = 1
        Int memoryGb = 4
        File chromLengths
        String sampleId
        File binestCov
    }

    command {
        somatic_tools/plot_bin_cov.R \
        "--binestOutput=~{binestCov}" \
        "--chrom_lengths=~{chromLengths}" \
        "--sample=~{sampleId}"
    }

    output {
        File normCoverageByChrPng = "~{sampleId}.binest.coverage.png"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/somatic_tools"  # "gcr.io/nygc-public/r:v3.2.1"
    }
}

task Pileup {
    input {
        Int threads
        Int memoryGb
        Int diskSize
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
        disks: "local-disk " + diskSize + " SSD"
        docker : "us.gcr.io/broad-gatk/gatk:4.0.0.0"
    }
}

task CalculateContamination {
    input {
        Int threads
        Int memoryGb
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
        docker : "us.gcr.io/broad-gatk/gatk:4.0.0.0"
    }
}

task CalculateContaminationPaired {
    input {
        Int threads
        Int memoryGb
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
        docker : "us.gcr.io/broad-gatk/gatk:4.0.0.0"
    }
}

task ConpairPileup {
    input {
        Int threads
        Int memoryGb
        Int diskSize
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
        disks: "local-disk " + diskSize + " SSD"
        docker : "us.gcr.io/broad-gatk/gatk:3.4.0"
    }
}

task VerifyConcordanceAll {
    input {
        Int threads
        Int memoryGb
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
        docker : "us.gcr.io/broad-gatk/gatk:4.1.0.0"
    }
}

task VerifyConcordanceHomoz {
    input {
        Int threads
        Int memoryGb
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
        docker : "us.gcr.io/broad-gatk/gatk:4.1.0.0"
    }
}

task Contamination {
    input {
        Int threads
        Int memoryGb
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
        docker : "us.gcr.io/broad-gatk/gatk:4.1.0.0"
    }
}
