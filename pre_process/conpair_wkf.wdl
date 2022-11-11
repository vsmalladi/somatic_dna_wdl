version 1.0

import "../wdl_structs.wdl"
import "qc.wdl"

workflow Conpair {
    # command
    input {
        File tumorPileupsConpair
        File normalPileupsConpair
        String tumor
        String normal
        String pairName

        File markerTxtFile

        Int threads=1
        Int memoryGb=4
        String qcDir = "Sample_~{tumor}/qc"
    }

    call qc.VerifyConcordanceAll {
        input:
            pileupsTumorConpair = tumorPileupsConpair,
            pileupsNormalConpair = normalPileupsConpair,
            markerTxtFile = markerTxtFile,
            pairName = pairName,
            memoryGb = memoryGb,
            threads = threads,
            concordanceAllPath = "~{qcDir}/~{pairName}.concordance.all.conpair.txt"
    }

    call qc.VerifyConcordanceHomoz {
        input:
            pileupsTumorConpair = tumorPileupsConpair,
            pileupsNormalConpair = normalPileupsConpair,
            markerTxtFile = markerTxtFile,
            pairName = pairName,
            memoryGb = memoryGb,
            threads = threads,
            concordanceHomozPath = "~{qcDir}/~{pairName}.concordance.homoz.conpair.txt"
    }

    call qc.Contamination {
        input:
            pileupsTumorConpair = tumorPileupsConpair,
            pileupsNormalConpair = normalPileupsConpair,
            markerTxtFile = markerTxtFile,
            pairName = pairName,
            memoryGb = memoryGb,
            threads = threads,
            contaminationPath = "~{qcDir}/~{pairName}.contamination.conpair.txt"
    }

    output {
        File concordanceAll = VerifyConcordanceAll.concordanceAll
        File concordanceHomoz = VerifyConcordanceHomoz.concordanceHomoz
        File contamination = Contamination.contamination
    }
}
