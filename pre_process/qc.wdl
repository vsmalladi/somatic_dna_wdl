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
        String outputDir = "."
        String MultipleMetricsBase = "~{outputDir}/~{sampleId}.MultipleMetrics"
    }

    command {
        mkdir -p $(dirname ~{MultipleMetricsBase})

        gatk CollectMultipleMetrics \
        --java-options "-Xmx32G -XX:ParallelGCThreads=1" \
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
        disks: "local-disk " + diskSize + " HDD"
        docker : "us.gcr.io/broad-gatk/gatk:4.1.0.0"
    }
}

task MultipleMetricsPreBqsr {
    input {
        Int threads = 2
        Int memoryGb = 40
        Int diskSize
        String outputDir = "."
        String MultipleMetricsBasePreBqsrBasename = "~{outputDir}/~{sampleId}.MultipleMetrics.dedup"
        IndexedReference referenceFa
        Bam mergedDedupBam
        String sampleId
    }

    command {
        mkdir -p $(dirname ~{MultipleMetricsBasePreBqsrBasename})

        gatk CollectMultipleMetrics \
        --java-options "-Xmx32G -XX:ParallelGCThreads=1" \
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
        disks: "local-disk " + diskSize + " HDD"
        docker : "us.gcr.io/broad-gatk/gatk:4.1.0.0"
    }
}

task CollectGcBiasMetrics {
    input {
        Int threads = 2
        Int memoryGb = 32
        Int diskSize
        String sampleId
        String outputDir = "."
        String gcBiasPdfPath = "~{outputDir}/~{sampleId}.GcBiasMetrics.gc_bias.pdf"
        String gcBiasMetricsPath = "~{outputDir}/~{sampleId}.GcBiasMetrics.gc_bias_metrics"
        String gcBiasSummaryPath = "~{outputDir}/~{sampleId}.GcBiasMetrics.gc_bias_summary"
        IndexedReference referenceFa
        Bam finalBam
    }

    command {
        mkdir -p $(dirname ~{gcBiasPdfPath})

        gatk CollectGcBiasMetrics \
        --java-options "-Xmx24G -XX:ParallelGCThreads=1" \
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
        disks: "local-disk " + diskSize + " HDD"
        docker : "us.gcr.io/broad-gatk/gatk:4.1.0.0"
    }
}

task Flagstat {
    input {
        Int threads = 2
        Int memoryGb = 32
        Int diskSize
        String sampleId
        String outputDir = "."
        String flagStatPath = "~{outputDir}/~{sampleId}.FlagStat.txt"
        IndexedReference referenceFa
        Bam finalBam
    }

    command {
        mkdir -p $(dirname ~{flagStatPath})

        gatk FlagStat \
        --java-options "-Xmx24G -XX:ParallelGCThreads=1" \
        --verbosity INFO \
        --reference ~{referenceFa.fasta} \
        -I ~{finalBam.bam} \
        > ~{flagStatPath}
    }

    output {
        File flagStat = "~{flagStatPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "us.gcr.io/broad-gatk/gatk:4.1.0.0"
    }
}

task HsMetrics {
    input {
        Int threads = 2
        Int memoryGb = 40
        Int diskSize
        String sampleId
        String outputDir = "."
        String hsMetricsPath = "~{outputDir}/~{sampleId}.HsMetrics.txt"
        String hsMetricsPerTargetCoveragePath = "~{outputDir}/~{sampleId}.HsMetrics.perTargetCoverage.txt"
        IndexedReference referenceFa
        Bam finalBam
        File hsMetricsIntervals
    }

    command {
        mkdir -p $(dirname ~{hsMetricsPath})

        gatk CollectHsMetrics \
        --java-options "-Xmx32G -XX:ParallelGCThreads=1" \
        --BAIT_INTERVALS ~{hsMetricsIntervals} \
        --TARGET_INTERVALS ~{hsMetricsIntervals} \
        --BAIT_SET_NAME ~{sampleId} \
        --METRIC_ACCUMULATION_LEVEL ALL_READS \
        --MINIMUM_MAPPING_QUALITY 1 \
        --MINIMUM_BASE_QUALITY 0 \
        --CLIP_OVERLAPPING_READS false \
        --REFERENCE_SEQUENCE ~{referenceFa.fasta} \
        -O ~{hsMetricsPath} \
        -I ~{finalBam.bam} \
        --PER_TARGET_COVERAGE ~{hsMetricsPerTargetCoveragePath} \
        --VALIDATION_STRINGENCY SILENT
    }

    output {
        File hsMetrics = "~{hsMetricsPath}"
        File hsMetricsPerTargetCoverage = "~{hsMetricsPerTargetCoveragePath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "us.gcr.io/broad-gatk/gatk:4.1.0.0"
    }
}

task FormatHsMetrics {
    input {
        Int threads = 1
        Int memoryGb = 4
        String sampleId
        String outputDir = '.'
        String hsMetricsPerTargetCoverageAutocorrPath = "~{outputDir}/~{sampleId}.HsMetrics.perTargetCoverage.txt.autocorr"
        File hsMetricsPerTargetCoverage
    }

    command {
        mkdir -p $(dirname ~{hsMetricsPerTargetCoverageAutocorrPath})

        perl /create_autocorrelation_input.v.0.1.pl \
        -input ~{hsMetricsPerTargetCoverage} \
        > ~{hsMetricsPerTargetCoverageAutocorrPath} \
    }

    output {
        File hsMetricsPerTargetCoverageAutocorr = "~{hsMetricsPerTargetCoverageAutocorrPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/somatic_tools:0.9.4"
    }
}

task Autocorrelations {
    input {
        Int threads = 1
        Int memoryGb = 4
        String sampleId
        String outputDir = '.'
        File hsMetricsPerTargetCoverageAutocorr
    }

    command {
        mkdir -p ~{outputDir}

        R --no-save \
        --args \
        ~{outputDir}/ \
        ~{hsMetricsPerTargetCoverageAutocorr} \
        ~{sampleId} \
        < /ASP_modified_final.v.0.1.R \
    }

    output {
        File autocorroutput1100 = "~{outputDir}/~{sampleId}.autocorroutput.1.100.txt"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/somatic_tools:0.9.4"
    }
}


task CollectOxoGMetricsWgs {
    input {
        Int threads = 2
        Int memoryGb = 8
        Int diskSize
        String sampleId
        String outputDir = "."
        String collectOxoGMetricsPath = "~{outputDir}/~{sampleId}.CollectOxoGMetrics.txt"
        IndexedReference referenceFa
        Bam finalBam
    }

    command {
        mkdir -p $(dirname ~{collectOxoGMetricsPath})

        gatk CollectOxoGMetrics \
        --java-options "-Xmx4G -XX:ParallelGCThreads=1" \
        --VALIDATION_STRINGENCY SILENT \
        -I ~{finalBam.bam} \
        -O ~{collectOxoGMetricsPath} \
        --REFERENCE_SEQUENCE ~{referenceFa.fasta}
    }

    output {
        File collectOxoGMetrics = "~{collectOxoGMetricsPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "us.gcr.io/broad-gatk/gatk:4.1.0.0"
    }
}

task CollectWgsMetrics {
    input {
        Int threads = 2
        Int memoryGb = 32
        Int diskSize
        String sampleId
        String outputDir = "."
        String collectWgsMetricsPath = "~{outputDir}/~{sampleId}.CollectWgsMetrics.txt"
        IndexedReference referenceFa
        Bam inputBam
        File randomIntervals
    }

    command {
        mkdir -p $(dirname ~{collectWgsMetricsPath})

        gatk CollectWgsMetrics \
        --java-options "-Xmx24G -XX:ParallelGCThreads=1" \
        --VALIDATION_STRINGENCY SILENT \
        -I ~{inputBam.bam} \
        -O ~{collectWgsMetricsPath} \
        --INTERVALS ~{randomIntervals} \
        -R ~{referenceFa.fasta} \
        --MINIMUM_MAPPING_QUALITY 0 \
        --COVERAGE_CAP 1000 \
        --MAX_RECORDS_IN_RAM 2000000 \
        --COUNT_UNPAIRED true \
        --MINIMUM_BASE_QUALITY 3
    }

    output {
        File collectWgsMetrics = "~{collectWgsMetricsPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "us.gcr.io/broad-gatk/gatk:4.1.0.0"
    }
}


task Binest {
    input {
        Int threads = 1
        Int memoryGb = 4
        Int diskSize
        String sampleId
        String outputDir = "."
        String binestCovPath = "~{outputDir}/~{sampleId}.binest.coverage.txt"
        Bam finalBam
    }

    command {
        mkdir -p $(dirname ~{binestCovPath})

        binest size ~{finalBam.bamIndex} > ~{binestCovPath}
    }

    output {
        File binestCov = "~{binestCovPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
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
        Rscript /plot_bin_cov.R \
        "--binest_output=~{binestCov}" \
        "--chrom_lengths=~{chromLengths}" \
        "--sample=~{sampleId}"
    }

    output {
        File normCoverageByChrPng = "~{sampleId}.binest.coverage.png"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/somatic_tools:0.9.4"  # "gcr.io/nygc-public/r:v3.2.1"
    }
}

task Pileup {
    input {
        Int memoryGb = 16
        Int diskSize
        String sampleId
        String outputDir = "."
        String pileupsTablePath = "~{outputDir}/~{sampleId}_pileups_table.table"
        Bam finalBam
        IndexedVcf gnomadBiallelic
    }

    command {
        mkdir -p $(dirname ~{pileupsTablePath})

        gatk GetPileupSummaries \
        --java-options "-Xmx12G -XX:ParallelGCThreads=1" \
        -I ~{finalBam.bam} \
        -V ~{gnomadBiallelic.vcf} \
        -O ~{pileupsTablePath}
    }

    output {
        File pileupsTable = "~{pileupsTablePath}"
    }

    runtime {
        cpu: 1
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "us.gcr.io/broad-gatk/gatk:4.0.0.0"
    }
}

task CalculateContamination {
    input {
        Int memoryGb = 16
        String sampleId
        String outputDir = "."
        String contaminationTablePath = "~{outputDir}/~{sampleId}.contamination.table"
        File pileupsTable
        Int diskSize
    }

    command {
        mkdir -p $(dirname ~{contaminationTablePath})

        gatk CalculateContamination \
        --java-options "-Xmx12G -XX:ParallelGCThreads=1" \
        -I ~{pileupsTable} \
        -O ~{contaminationTablePath}
    }

    output {
        File contaminationTable = "~{contaminationTablePath}"
    }

    runtime {
        cpu : 1
        memory : memoryGb + "GB"
        docker : "us.gcr.io/broad-gatk/gatk:4.0.0.0"
        disks: "local-disk " + diskSize + " HDD"
    }
}

task CalculateContaminationPaired {
    input {
        Int threads
        Int memoryGb = 8
        String pairName
        String outputDir = "."
        String contaminationTablePath = "~{outputDir}/~{pairName}.contamination.table"
        File pileupsNormalTable
        File pileupsTumorTable
    }

    command {
        mkdir -p $(dirname ~{contaminationTablePath})

        gatk CalculateContamination \
        --java-options "-Xmx4G -XX:ParallelGCThreads=1" \
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
        mkdir -p $(dirname ~{pileupsConpairPath})

        java \
        -jar GenomeAnalysisTK.jar \
        --java-options "-Xmx4G -XX:ParallelGCThreads=1" \
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
        disks: "local-disk " + diskSize + " HDD"
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
