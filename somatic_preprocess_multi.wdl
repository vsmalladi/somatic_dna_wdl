version 1.0
import "./wdl_structs.wdl"
import "pre_process/pre_process.wdl" as preprocess
import "qc/qc_wkf.wdl" as qc

workflow SomaticWorkflow {
    input {
        BwaReference bwaReference
        IndexedReference indexedReference
        IndexedVcf MillsAnd1000G
        IndexedVcf Indels
        IndexedVcf DbSnp
        File chromFile
        File chromLengths  # Is this different from above?
        File hsMetricsIntervals
        File randomIntervals
        Array[Array[Fastqs]]+ tumorFastqs
        # normalFastqs can be missing for some.
        Array[Array[Fastqs]] normalFastqs
        Array[String] tumorIds
        Array[String] normalIds

    }

    scatter (i in range(length(tumorFastqs))) {
        call preprocess.Preprocess as tumorPrep {
            input:
                listOfFastqPairs = tumorFastqs[i],
                sampleId = tumorIds[i],
                bwaReference = bwaReference,
                indexedReference = indexedReference,
                MillsAnd1000G = MillsAnd1000G,
                Indels = Indels,
                DbSnp = DbSnp,
                chromFile = chromFile
        }

        call qc.QcMetrics as tumorQc {
            input:
                sampleId = tumorIds[i],
                mergedDedupBam = tumorPrep.mergedDedupBam,
                finalBam = tumorPrep.finalBam,
                referenceFa = indexedReference,
                hsMetricsIntervals = hsMetricsIntervals,
                randomIntervals = randomIntervals,
                chromLengths = chromLengths
        }
    }

    scatter (i in range(length(normalFastqs))) {
        call preprocess.Preprocess as normalPrep {
            input:
                listOfFastqPairs = normalFastqs[i],
                sampleId = normalIds[i],
                bwaReference = bwaReference,
                indexedReference = indexedReference,
                MillsAnd1000G = MillsAnd1000G,
                Indels = Indels,
                DbSnp = DbSnp,
                chromFile = chromFile
        }

        call qc.QcMetrics as normalQc {
            input:
                sampleId = normalIds[i],
                mergedDedupBam = normalPrep.mergedDedupBam,
                finalBam = normalPrep.finalBam,
                referenceFa = indexedReference,
                hsMetricsIntervals = hsMetricsIntervals,
                hsMetricsIntervals = hsMetricsIntervals,
                randomIntervals = randomIntervals,
                chromLengths = chromLengths
        }
    }

    output {
        Array[Bam] tumorBams = tumorPrep.finalBam
        Array[Array[File]] tumorQCFiles = tumorQc.QcFiles
        # QC files for samples.
        Array[Bam] normalBams = normalPrep.finalBam
        Array[Array[File]] normalQCFiles = normalQc.QcFiles
    }
}
