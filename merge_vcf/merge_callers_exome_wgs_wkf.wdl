version 1.0

import "merge_vcf.wdl" as mergeVcf
import "../calling/calling.wdl" as calling
import "../wdl_structs.wdl"

# note that we will still need to match dictionaries or provide matching reference files

workflow MergeCallers {
    input {
        Boolean external = false
        String tumorId
        String normalId
        String pairName
        Array[String]+ listOfChroms
        Array[IndexedVcf]+ allVcfCompressed

        IndexedReference referenceFa
    }

    scatter(vcfCompressed in allVcfCompressed) {
            File allVcfCompressedFile = vcfCompressed.vcf
        }
    Array[File]+ allVcfCompressedList = allVcfCompressedFile

    scatter(chrom in listOfChroms) {
        call mergeVcf.MergeExomeWgsCallers as allMergeExomeWgsCallers {
            input:
                chrom = chrom,
                pairName = pairName,
                allVcfCompressed = allVcfCompressed,
                allVcfCompressedList = allVcfCompressedList
        }

        #  =================================================================
        #                     Merge columns
        #  =================================================================

        call mergeVcf.MergeColumns {
            input:
                chrom = chrom,
                tumorId = tumorId,
                normalId = normalId,
                pairName = pairName,
                supportedChromVcf = allMergeExomeWgsCallers.mergedChromVcf
        }

        call mergeVcf.SnvstomnvsAnnotateExomeWgsCalled {
            input:
                chrom = chrom,
                pairName = pairName,
                filteredOutFile = MergeColumns.columnChromVcf
        }

    }

    output {
            Array[File] mergedChromVcf = MergeColumns.columnChromVcf
            Array[File] finalChromVcf = SnvstomnvsAnnotateExomeWgsCalled.finalChromVcf
        }
}
