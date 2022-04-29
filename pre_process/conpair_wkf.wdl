version 1.0

import "../wdl_structs.wdl"
import "qc.wdl"

workflow Conpair {
    # command
    input {
        Bam finalNormalBam
        Bam finalTumorBam
        IndexedReference referenceFa
        String tumor
        String normal
        String pairName

        File markerTxtFile
        File markerBedFile

        Int threads=1
        Int memoryGb=4
        String qcDir = "Sample_~{tumor}/qc"
    }

    Int additionalDiskSize = 10
    Int tumorSize = ceil(size(finalTumorBam.bam, "GB") + size(finalTumorBam.bamIndex, "GB"))
    Int normalSize = ceil(size(finalNormalBam.bam, "GB") + size(finalNormalBam.bamIndex, "GB"))

    call qc.ConpairPileup as tumorConpairPileup {
        input:
            markerBedFile = markerBedFile,
            referenceFa = referenceFa,
            finalBam = finalTumorBam,
            sampleId = tumor,
            memoryGb = memoryGb,
            threads = threads,
            diskSize = tumorSize + additionalDiskSize
    }

    call qc.ConpairPileup as normalConpairPileup {
        input:
            markerBedFile = markerBedFile,
            referenceFa = referenceFa,
            finalBam = finalNormalBam,
            sampleId = normal,
            memoryGb = memoryGb,
            threads = threads,
            diskSize = normalSize + additionalDiskSize
    }

    call qc.VerifyConcordanceAll {
        input:
            pileupsTumorConpair = tumorConpairPileup.pileupsConpair,
            pileupsNormalConpair = normalConpairPileup.pileupsConpair,
            markerTxtFile = markerTxtFile,
            pairName = pairName,
            memoryGb = memoryGb,
            threads = threads,
            concordanceAllPath = "~{qcDir}/~{pairName}.concordance.all.conpair.txt"
    }

    call qc.VerifyConcordanceHomoz {
        input:
            pileupsTumorConpair = tumorConpairPileup.pileupsConpair,
            pileupsNormalConpair = normalConpairPileup.pileupsConpair,
            markerTxtFile = markerTxtFile,
            pairName = pairName,
            memoryGb = memoryGb,
            threads = threads,
            concordanceHomozPath = "~{qcDir}/~{pairName}.concordance.homoz.conpair.txt"
    }

    call qc.Contamination {
        input:
            pileupsTumorConpair = tumorConpairPileup.pileupsConpair,
            pileupsNormalConpair = normalConpairPileup.pileupsConpair,
            markerTxtFile = markerTxtFile,
            pairName = pairName,
            memoryGb = memoryGb,
            threads = threads,
            contaminationPath = "~{qcDir}/~{pairName}.contamination.conpair.txt"
    }

    output {
        File tumorPileup = tumorConpairPileup.pileupsConpair
        File normalPileup = normalConpairPileup.pileupsConpair
        File concordanceAll = VerifyConcordanceAll.concordanceAll
        File concordanceHomoz = VerifyConcordanceHomoz.concordanceHomoz
        File contamination = Contamination.contamination
    }
}
