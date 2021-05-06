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
        Map[Int, Array[File]] uniqCoords
        
        File tumorConfigFile
        File normalConfigFile
        Int tumorMedianInsertSize
        Int normalMedianInsertSize
        String tumorParamsPath
        String normalParamsPath
        Map[String, File] chromFastas
        
        IndexedReference referenceFa
        
        Int threads
        Int memoryGb
        String samtoolsDockerImage
        String bicseq2DockerImage
        String gatkDockerImage
        String pysamDockerImage
        File jsonLog
    }
    
    scatter (chrom in listOfChroms) {
            String tempNormalSeq = "~{normal}_~{chrom}.seq"
            String tempTumorSeq = "~{tumor}_~{chrom}.seq"
            File chromFastaFile = chromFastas[chrom]
        }
    Array[String] tempNormalSeqsPaths = tempNormalSeq
    Array[String] tempTumorSeqsPaths = tempTumorSeq
    Array[File] chromFastaFiles = chromFastaFile
    
    call calling.UniqReads uniqReadsNormal {
        input:
            sampleId = normal,
            finalBam = normalFinalBam,
            tempSeqsPaths = tempNormalSeqsPaths,
            memoryGb = 8
    }
    
    call calling.UniqReads uniqReadsTumor {
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
            tempSeqs = UniqReads.tempTumorSeqs,
            configFile = tumorConfigFile,
            sampleId = tumor,
            readLength = readLength,
            medianInsertSize = tumorMedianInsertSize,
            paramsPath = tumorParamsPath,
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
            tempSeqs = UniqReads.tempNormalSeqs,
            configFile = normalConfigFile,
            sampleId = normal,
            readLength = readLength,
            medianInsertSize = normalMedianInsertSize,
            paramsPath = normalParamsPath,
            tempNormPaths = tempNormalNormPaths,
            chromFastas = chromFastaFiles,
            uniqCoords = uniqCoords[coordReadLength],
            memoryGb = 8
    }
    
    call calling.Bicseq2Wgs {
        input:
            tempTumorNorms = tumorBicseq2Norm.tempNorm,
            tempNormalNorms = normalBicseq2Norm.tempNorm,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = bicseq2DockerImage
    }
    
    output {
        File bicseq2Png = Bicseq2Wgs.bicseq2Png
        File bicseq2 = Bicseq2Wgs.bicseq2
    }
}
    