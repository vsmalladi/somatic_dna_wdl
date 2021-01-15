version 1.0

task GatkCollectHsMetrics {
    input {
        File bam
        File bamIndex
        File baitIntervals
        String sample
        File targetIntervals
        File reference
        File referenceIndex

        #String docker_tag
        Int memory_gb
    }
# Check which of these are needed long with the reference.
#    File referenceAnn
#    File referenceBwt
#    File referencePac
#    File referenceSa
#    File referenceDict

    command {
        gatk CollectHsMetrics \
            -I ~{bam} \
            -O ~{sample}.HsMetrics.txt \
            --BAIT_SET_NAME ~{sample} \
            --PER_TARGET_COVERAGE ~{sample}.HsMetrics.perTargetCoverage.txt \
            --METRIC_ACCUMULATION_LEVEL ALL_READS \
            --VALIDATION_STRINGENCY SILENT \
            --BAIT_INTERVALS ~{baitIntervals} \
            --TARGET_INTERVALS ~{targetIntervals} \
            --REFERENCE_SEQUENCE ~{reference}
    }
    runtime {
        memory: memory_gb + "GB"
        cpu: 1
        #docker: gatkImage
    }
    output {
        File outputHsMetrics = "~{sample}.HsMetrics.txt"
        File outputPerTargetCoverage = "~{sample}.HsMetrics.perTargetCoverage.txt"
    }
}

task GatkCollectMultipleMetrics {
    input {
        File bam
        File bamIndex
        String sample
        File reference
        File referenceIndex
    }
    command {
        gatk CollectMultipleMetrics \
            -I ~{bam} \
            -O ~{sample}.MultipleMetrics \
            --VALIDATION_STRINGENCY SILENT \
            --REFERENCE_SEQUENCE ~{reference} \
            --PROGRAM null \
            --PROGRAM CollectAlignmentSummaryMetrics \
            --PROGRAM MeanQualityByCycle \
            --PROGRAM CollectInsertSizeMetrics \
            --PROGRAM QualityScoreDistribution
   }
   runtime {
      memory: "16 GB"
      cpu: "1"
   }
   output {
       File multiple_metrics_alignment_summary = "~{sample}.MultipleMetrics.alignment_summary_metrics"
       File multiple_metrics_base_distribution_by_cycle_metrics = "~{sample}.MultipleMetrics.base_distribution_by_cycle_metrics"
       File multiple_metrics_base_distribution_by_cycle_pdf = "~{sample}.MultipleMetrics.base_distribution_by_cycle.pdf"
       File multiple_metrics_insert_size_metrics =  "~{sample}.MultipleMetrics.insert_size_metrics"
       File multiple_metrics_insert_size_histogram_pdf = "~{sample}.MultipleMetrics.insert_size_histogram.pdf"
       File multiple_metrics_quality_by_cycle_metrics = "~{sample}.MultipleMetrics.quality_by_cycle_metrics"
       File multiple_metrics_quality_by_cycle_pdf =  "~{sample}.MultipleMetrics.quality_by_cycle.pdf"
       File multiple_metrics_quality_distribution_metrics = "~{sample}.MultipleMetrics.quality_distribution_metrics"
       File multiple_metrics_quality_distribution_pdf = "~{sample}.MultipleMetrics.quality_distribution.pdf"
  }
}

task GatkFlagStat {
    input {
        File bam
        File bamIndex
        File reference
        File reference_ann
        File reference_bwt
        File reference_pac
        File reference_sa
        File referenceIndex
        File reference_dict
        String sample
    }
    command {
        gatk FlagStat \
            -I ~{bam} \
            -R ~{reference} \
            --verbosity INFO >> \
            ~{sample}.FlagStat.txt \
    }
    runtime {
        memory: "16 GB"
        cpus: "1"
    }
    output {
        File output_flagstat = "~{sample}.FlagStat.txt"
    }
}

task GatkCollectWgsMetrics {
    input {
        File bam
        File bamIndex
        String sample
        File reference
        File referenceIndex
        File intervals
    }

    command {
        gatk CollectWgsMetrics \
            -I ~{bam} \
            -O ~{sample}.CollectWgsMetrics.txt \
            -R ~{reference} \
            --MAX_RECORDS_IN_RAM 2000000
            --INTERVALS ~{intervals}
            --MINIMUM_MAPPING_QUALITY 0
            --COUNT_UNPAIRED true
            --MINIMUM_BASE_QUALITY 3
            --COVERAGE_CAP 1000
            --VALIDATION_STRINGENCY SILENT
    }
    runtime {
        memory: "16 GB"
        cpu: "1"
    }
    output {
        File output_wgs_metrics = "~{sample}.CollectWgsMetrics.txt"
    }
}

task GatkCollectGcBiasMetrics {
    input {
        File bam
        File bamIndex
        String sample
        File reference
        File referenceIndex
    }

    command {
        gatk CollectGcBiasMetrics \
            --INPUT ~{bam} \
            --OUTPUT ~{sample}.GcBiasMetrics.gc_bias_metrics \
            --CHART_OUTPUT ~{sample}.GcBiasMetrics.gc_bias.pdf \
            --SUMMARY_OUTPUT ~{sample}.GcBiasMetrics.gc_bias_summary \
            --REFERENCE_SEQUENCE ~{reference} \
            --ASSUME_SORTED true \
            --VALIDATION_STRINGENCY SILENT
    }
    runtime {
        memory: "16 GB"
        cpu: "1"
    }
    output {
        File gc_bias_metrics = "~{sample}.GcBiasMetrics.gc_bias_metrics"
        File chart_output = "~{sample}.GcBiasMetrics.gc_bias.pdf"
        File summary_output = "~{sample}.GcBiasMetrics.gc_bias_summary"
    }
}

task GatkCollectOxoGMetrics {
    input {
        File bam
        File bamIndex
        String sample
        File reference
        File referenceIndex
    }

    command {
        gatk CollectOxoGMetrics \
            --INPUT ~{bam} \
            --OUTPUT ~{sample}.CollectOxoGMetrics.txt \
            --REFERENCE ~{reference} \
            --VALIDATION_STRINGENCY=SILENT
    }
    runtime {

    }
    output {
        File output_oxog_metrics = "~{sample}.CollectOxoGMetrics.txt"
    }
}

task GatkGetPileupSummaries {
    input {
        File bam
        File bamIndex
        String sample
        File gnomadVcf
    }

    command {
        gatk GetPileupSummaries \
            --INPUT ~{bam} \
            --OUTPUT ~{sample}.pileups_table.table \
            --variant ~{gnomadVcf} \
            --VALIDATION_STRINGENCY SILENT
    }
    runtime {

    }
    output {
        File output_pileup_table  = "~{sample}.pileups_table.table"
    }
}

task GatkCalculateContamination {
    input {
        File bam
        File bamIndex
        String sample
    }

    command {
        gatk GetPileupSummaries \
            --INPUT ~{bam} \
            --OUTPUT ~{sample}.contamination.table \
    }
    runtime {

    }
    output {
        File output_contamination_table  = "~{sample}.contamination.table"
    }
}
