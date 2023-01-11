version 1.0

import "../wdl_structs.wdl"

task MultipleMetrics {
    input {
        Int threads = 2
        Int memoryGb = 24
        Int diskSize
        IndexedReference referenceFa
        Bam finalBam
        String sampleId
        String outputDir = "."
        String MultipleMetricsBase = "~{outputDir}/~{sampleId}.MultipleMetrics"
    }

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)

    command {
        mkdir -p $(dirname ~{MultipleMetricsBase})

        gatk CollectMultipleMetrics \
        --java-options "-Xmx~{jvmHeap}m -XX:ParallelGCThreads=1" \
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
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-public/broadinstitute/gatk4@sha256:b3bde7bc74ab00ddce342bd511a9797007aaf3d22b9cfd7b52f416c893c3774c"
    }
}

task MultipleMetricsPreBqsr {
    input {
        Int threads = 2
        Int memoryGb = 16
        Int diskSize
        String outputDir = "."
        String MultipleMetricsBasePreBqsrBasename = "~{outputDir}/~{sampleId}.MultipleMetrics.dedup"
        IndexedReference referenceFa
        Bam mergedDedupBam
        String sampleId
    }

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)
    command {
        mkdir -p $(dirname ~{MultipleMetricsBasePreBqsrBasename})

        gatk CollectMultipleMetrics \
        --java-options "-Xmx~{jvmHeap}m -XX:ParallelGCThreads=1" \
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
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-public/broadinstitute/gatk4@sha256:b3bde7bc74ab00ddce342bd511a9797007aaf3d22b9cfd7b52f416c893c3774c"
    }
}

task CollectGcBiasMetrics {
    input {
        Int threads = 2
        Int memoryGb = 16
        Int diskSize
        String sampleId
        String outputDir = "."
        String gcBiasPdfPath = "~{outputDir}/~{sampleId}.GcBiasMetrics.gc_bias.pdf"
        String gcBiasMetricsPath = "~{outputDir}/~{sampleId}.GcBiasMetrics.gc_bias_metrics"
        String gcBiasSummaryPath = "~{outputDir}/~{sampleId}.GcBiasMetrics.gc_bias_summary"
        IndexedReference referenceFa
        Bam finalBam
    }

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)
    command {
        mkdir -p $(dirname ~{gcBiasPdfPath})

        gatk CollectGcBiasMetrics \
        --java-options "-Xmx~{jvmHeap}m -XX:ParallelGCThreads=1" \
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
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-public/broadinstitute/gatk4@sha256:b3bde7bc74ab00ddce342bd511a9797007aaf3d22b9cfd7b52f416c893c3774c"
    }
}

task Flagstat {
    input {
        Int threads = 2
        Int memoryGb = 8
        Int diskSize
        String sampleId
        String outputDir = "."
        String flagStatPath = "~{outputDir}/~{sampleId}.FlagStat.txt"
        IndexedReference referenceFa
        Bam finalBam
    }

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)
    command {
        mkdir -p $(dirname ~{flagStatPath})

        gatk FlagStat \
        --java-options "-Xmx~{jvmHeap}m -XX:ParallelGCThreads=1" \
        --verbosity INFO \
        --reference ~{referenceFa.fasta} \
        -I ~{finalBam.bam} \
        > ~{flagStatPath}
    }

    output {
        File flagStat = "~{flagStatPath}"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-public/broadinstitute/gatk4@sha256:b3bde7bc74ab00ddce342bd511a9797007aaf3d22b9cfd7b52f416c893c3774c"
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

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)
    command {
        mkdir -p $(dirname ~{hsMetricsPath})

        gatk CollectHsMetrics \
        --java-options "-Xmx~{jvmHeap}m -XX:ParallelGCThreads=1" \
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
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-public/broadinstitute/gatk4@sha256:b3bde7bc74ab00ddce342bd511a9797007aaf3d22b9cfd7b52f416c893c3774c"
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
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:2d61578c9bc8e5ef4a4b25c1e6a32883c67e0fcff6259e43c7da91dbc9c100d7"
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
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:2d61578c9bc8e5ef4a4b25c1e6a32883c67e0fcff6259e43c7da91dbc9c100d7"
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

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)
    command {
        mkdir -p $(dirname ~{collectOxoGMetricsPath})

        gatk CollectOxoGMetrics \
        --java-options "-Xmx~{jvmHeap}m -XX:ParallelGCThreads=1" \
        --VALIDATION_STRINGENCY SILENT \
        -I ~{finalBam.bam} \
        -O ~{collectOxoGMetricsPath} \
        --REFERENCE_SEQUENCE ~{referenceFa.fasta}
    }

    output {
        File collectOxoGMetrics = "~{collectOxoGMetricsPath}"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-public/broadinstitute/gatk4@sha256:b3bde7bc74ab00ddce342bd511a9797007aaf3d22b9cfd7b52f416c893c3774c"
    }
}

task CollectWgsMetrics {
    input {
        Int threads = 2
        Int memoryGb = 16
        Int diskSize
        String sampleId
        String outputDir = "."
        String collectWgsMetricsPath = "~{outputDir}/~{sampleId}.CollectWgsMetrics.txt"
        IndexedReference referenceFa
        Bam inputBam
        File randomIntervals
    }

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)
    command {
        mkdir -p $(dirname ~{collectWgsMetricsPath})

        gatk CollectWgsMetrics \
        --java-options "-Xmx~{jvmHeap}m -XX:ParallelGCThreads=1" \
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
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-public/broadinstitute/gatk4@sha256:b3bde7bc74ab00ddce342bd511a9797007aaf3d22b9cfd7b52f416c893c3774c"
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
        String binestSexPath = "~{outputDir}/~{sampleId}.binest.sex.txt"
        File finalBamIndex
        File referenceFai
    }

    command {
        mkdir -p $(dirname ~{binestCovPath})

        binest -f ~{referenceFai} size ~{finalBamIndex} > ~{binestCovPath}
        binest -f ~{referenceFai} sex ~{finalBamIndex} > ~{binestSexPath}
    }

    output {
        File binestCov = "~{binestCovPath}"
        File binestSex = "~{binestSexPath}"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-public/binest@sha256:86aa44cd68a70e811d4f6250dd1fafbd49c6b2a5e59462f131753cb9f0f5be1a"
    }
}

task PlotBinCov {
    input {
        Int threads = 1
        Int memoryGb = 4
        File chromLengths
        String sampleId
        File binestCov
        String outputDir = "."
    }

    command {
        mkdir -p ~{outputDir}

        cd ~{outputDir}
        Rscript /plot_bin_cov.R \
        "--binest_output=~{binestCov}" \
        "--chrom_lengths=~{chromLengths}" \
        "--sample=~{sampleId}"
    }

    output {
        File normCoverageByChrPng = "~{outputDir}/~{sampleId}.binest.coverage.png"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:2d61578c9bc8e5ef4a4b25c1e6a32883c67e0fcff6259e43c7da91dbc9c100d7"
    }
}

task Pileup {
    input {
        Int memoryGb = 16
        Int threads = 1
        Int diskSize
        String sampleId
        String outputDir = "."
        String pileupsTablePath = "~{outputDir}/~{sampleId}_pileups_table.table"
        Bam finalBam
        IndexedVcf gnomadBiallelic
    }

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)
    command {
        mkdir -p $(dirname ~{pileupsTablePath})

        gatk GetPileupSummaries \
        --java-options "-Xmx~{jvmHeap}m -XX:ParallelGCThreads=1" \
        -I ~{finalBam.bam} \
        -V ~{gnomadBiallelic.vcf} \
        -L ~{gnomadBiallelic.vcf} \
        -O ~{pileupsTablePath}
    }

    output {
        File pileupsTable = "~{pileupsTablePath}"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu: threads
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-public/broadinstitute/gatk4@sha256:b3bde7bc74ab00ddce342bd511a9797007aaf3d22b9cfd7b52f416c893c3774c"
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
        Int threads = 1
    }

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)
    command {
        mkdir -p $(dirname ~{contaminationTablePath})

        gatk CalculateContamination \
        --java-options "-Xmx~{jvmHeap}m -XX:ParallelGCThreads=1" \
        -I ~{pileupsTable} \
        -O ~{contaminationTablePath}
    }

    output {
        File contaminationTable = "~{contaminationTablePath}"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/broadinstitute/gatk4@sha256:b3bde7bc74ab00ddce342bd511a9797007aaf3d22b9cfd7b52f416c893c3774c"
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

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)
    command {
        mkdir -p $(dirname ~{contaminationTablePath})

        gatk CalculateContamination \
        --java-options "-Xmx~{jvmHeap}m -XX:ParallelGCThreads=1" \
        -I ~{pileupsTumorTable} \
        -matched ~{pileupsNormalTable} \
        -O ~{contaminationTablePath}
    }

    output {
        File contaminationTable = "~{contaminationTablePath}"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/broadinstitute/gatk4@sha256:b3bde7bc74ab00ddce342bd511a9797007aaf3d22b9cfd7b52f416c893c3774c"
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

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)
    command {
        mkdir -p $(dirname ~{pileupsConpairPath})

        java \
        -Xmx~{jvmHeap}m -XX:ParallelGCThreads=1 \
        -jar /usr/GenomeAnalysisTK.jar \
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
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-public/broadinstitute/gatk3@sha256:9f72be83047bf9774c6afb091d622c6e7e0c8e94111f4acc745a4e70b7a1b965"
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
        String concordanceAllPath = "~{pairName}.concordance.all.conpair.txt"
    }

    command {
        mkdir -p $(dirname ~{concordanceAllPath})
        export CONPAIR_DIR=/Conpair-0.2
        export PYTHONPATH=/Conpair-0.2/modules
        python /Conpair-0.2/scripts/verify_concordance.py \
        -T ~{pileupsTumorConpair} \
        -N ~{pileupsNormalConpair} \
        -O ~{concordanceAllPath} \
        -M ~{markerTxtFile}
    }

    output {
        File concordanceAll = "~{concordanceAllPath}"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/conpair@sha256:c427e731b280ba8ad83c79993970bab5e94dd42f2b972f580b817c4a1f88e8c9"
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
        String concordanceHomozPath = "~{pairName}.concordance.homoz.conpair.txt"
    }

    command {
        mkdir -p $(dirname ~{concordanceHomozPath})

        export CONPAIR_DIR=/Conpair-0.2
        export PYTHONPATH=/Conpair-0.2/modules
        python /Conpair-0.2/scripts/verify_concordance.py \
        -T ~{pileupsTumorConpair} \
        -N ~{pileupsNormalConpair} \
        -O ~{concordanceHomozPath} \
        -M ~{markerTxtFile} \
        -H
    }

    output {
        File concordanceHomoz = "~{concordanceHomozPath}"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/conpair@sha256:c427e731b280ba8ad83c79993970bab5e94dd42f2b972f580b817c4a1f88e8c9"
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
        String contaminationPath = "~{pairName}.contamination.conpair.txt"
    }

    command {
        mkdir -p $(dirname ~{contaminationPath})

        export CONPAIR_DIR=/Conpair-0.2
        export PYTHONPATH=/Conpair-0.2/modules
        python /Conpair-0.2/scripts/estimate_tumor_normal_contamination.py \
        -T ~{pileupsTumorConpair} \
        -N ~{pileupsNormalConpair} \
        -O ~{contaminationPath} \
        -M ~{markerTxtFile}
    }

    output {
        File contamination = "~{contaminationPath}"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/conpair@sha256:c427e731b280ba8ad83c79993970bab5e94dd42f2b972f580b817c4a1f88e8c9"
    }
}
