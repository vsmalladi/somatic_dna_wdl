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
        Array[String] listOfChroms
        
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
    
    scatter (chrom in listOfChroms) {
            String tempNormalSeq = "~{normal}_~{chrom}.seq"
            String tempTumorSeq = "~{tumor}_~{chrom}.seq"
            File chromFastaFile = chromFastas[chrom]
        }
    Array[String] tempNormalSeqsPaths = tempNormalSeq
    Array[String] tempTumorSeqsPaths = tempTumorSeq
    Array[File] chromFastaFiles = chromFastaFile
    
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
    
    scatter(chrom in listOfChroms) {
        String tempTumorNormFile = "~{tumor}/~{tumor}_~{chrom}.norm.bin.txt"
    }
    Array[String]+ tempTumorNormPaths = tempTumorNormFile
    
    call calling.Bicseq2Norm as tumorBicseq2Norm {
        input:
            sampleId = tumor,
            tempSeqs = uniqReadsTumor.tempSeqs,
            configFile = tumorConfigFile,
            sampleId = tumor,
            readLength = readLength,
            medianInsertSize = tumorMedianInsertSize,
            tempNormPaths = tempTumorNormPaths,
            chromFastas = chromFastaFiles,
            uniqCoords = uniqCoords[coordReadLength],
            memoryGb = 8
    }
    
    scatter(chrom in listOfChroms) {
        File tempNormalNormFile = "~{normal}/~{normal}_~{chrom}.norm.bin.txt"
    }
    Array[File]+ tempNormalNormPaths = tempNormalNormFile
    
    call calling.Bicseq2Norm as normalBicseq2Norm {
        input:
            sampleId = normal,
            tempSeqs = uniqReadsNormal.tempSeqs,
            configFile = normalConfigFile,
            sampleId = normal,
            readLength = readLength,
            medianInsertSize = normalMedianInsertSize,
            tempNormPaths = tempNormalNormPaths,
            chromFastas = chromFastaFiles,
            uniqCoords = uniqCoords[coordReadLength],
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
    