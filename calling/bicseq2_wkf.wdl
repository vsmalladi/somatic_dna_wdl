version 1.0

import "calling.wdl" as calling
import "../wdl_structs.wdl"

workflow BicSeq2 {
    # command 
    #   run BicSeq2 caller
    input {
        String tumor
        String normal
        String pairName
        Array[String] listOfChromsFull
        
        Bam normalFinalBam
        Bam tumorFinalBam
        Int readLength
        Int coordReadLength
        Map[Int, Map[String, File]] uniqCoords
        
        File tumorConfigFile
        File normalConfigFile
        Int tumorMedianInsertSize = 400
        Int normalMedianInsertSize = 400
        Map[String, File] chromFastas
        
        IndexedReference referenceFa
        
        File segConfigFile
        
        Int lambda = 4
    }
    
    scatter (chrom in listOfChromsFull) {
            String tempNormalSeq = "~{normal}_~{chrom}.seq"
            String tempTumorSeq = "~{tumor}_~{chrom}.seq"
            File chromFastasFile = chromFastas[chrom]
        }
    Array[String] tempNormalSeqsPaths = tempNormalSeq
    Array[String] tempTumorSeqsPaths = tempTumorSeq
    Array[File] chromFastasFiles = chromFastasFile
    
    call calling.UniqReads as uniqReadsNormal {
        input:
            sampleId = normal,
            finalBam = normalFinalBam,
            tempSeqsPaths = tempNormalSeqsPaths,
            memoryGb = 8
    }
    
    call calling.UniqReads as uniqReadsTumor {
        input:
            sampleId = tumor,
            finalBam = tumorFinalBam,
            tempSeqsPaths = tempTumorSeqsPaths,
            memoryGb = 8
    }
    
    scatter(chrom in listOfChromsFull) {
        String tempTumorNormFile = "~{tumor}/~{tumor}_~{chrom}.norm.bin.txt"
    }
    Array[String]+ tempTumorNormPaths = tempTumorNormFile
    
    # prep chrom resources
    scatter(chrom in listOfChromsFull) {
        File uniqCoordFile = uniqCoords[coordReadLength][chrom]
    }
    Array[File] uniqCoordsFiles = uniqCoordFile
    
    call calling.Bicseq2Norm as tumorBicseq2Norm {
        input:
            sampleId = tumor,
            tempSeqs = uniqReadsTumor.tempSeqs,
            configFile = tumorConfigFile,
            sampleId = tumor,
            readLength = readLength,
            medianInsertSize = tumorMedianInsertSize,
            tempNormPaths = tempTumorNormPaths,
            uniqCoordsFiles=uniqCoordsFiles,
            chromFastasFiles=chromFastasFiles,
            memoryGb = 8
    }
    
    scatter(chrom in listOfChromsFull) {
        String tempNormalNormFile = "~{normal}/~{normal}_~{chrom}.norm.bin.txt"
    }
    Array[String]+ tempNormalNormPaths = tempNormalNormFile
    
    call calling.Bicseq2Norm as normalBicseq2Norm {
        input:
            sampleId = normal,
            tempSeqs = uniqReadsNormal.tempSeqs,
            configFile = normalConfigFile,
            sampleId = normal,
            readLength = readLength,
            medianInsertSize = normalMedianInsertSize,
            tempNormPaths = tempNormalNormPaths,
            uniqCoordsFiles=uniqCoordsFiles,
            chromFastasFiles=chromFastasFiles,
            memoryGb = 8
    }
    
    call calling.Bicseq2Wgs {
        input:
            pairName = pairName,
            tempTumorNorms = tumorBicseq2Norm.tempNorm,
            tempNormalNorms = normalBicseq2Norm.tempNorm,
            segConfigFile = segConfigFile,
            lambda = lambda,
            memoryGb = 8
    }
    
    output {
        File bicseq2Png = Bicseq2Wgs.bicseq2Png
        File bicseq2 = Bicseq2Wgs.bicseq2
    }
}
    