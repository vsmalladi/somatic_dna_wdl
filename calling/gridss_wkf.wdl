version 1.0

import "calling.wdl" as calling
import "../wdl_structs.wdl"

workflow Gridss {
    # command 
    #   run Gridss caller
    input {
        String tumor
        String normal
        String pairName
        
        BwaReference gridssReferenceFa
        Array[File] gridssAdditionalReference
                
        Bam normalFinalBam
        Bam tumorFinalBam
        
        Int assembleChunks = 5
        
        String bsGenome
        File ponTarGz
        
        Int threads = 8
        Int memory_gb = 32
    }
    
    call calling.GridssPreprocess as tumorGridssPreprocess {
        input:
            threads = threads,
            memory_gb = memory_gb,
            gridssReferenceFa = gridssReferenceFa,
            gridssAdditionalReference = gridssAdditionalReference,
            finalBam = tumorFinalBam  
    }
    
    call calling.GridssPreprocess as normalGridssPreprocess {
        input:
            threads = threads,
            memory_gb = memory_gb,
            gridssReferenceFa = gridssReferenceFa,
            gridssAdditionalReference = gridssAdditionalReference,
            finalBam = normalFinalBam  
    }
    
    scatter(i in range(assembleChunks)) {
        call calling.GridssAssembleChunk {
            input:
                threads = threads,
                memory_gb = 100,
                pairName = pairName,
                gridssReferenceFa = gridssReferenceFa,
                gridssAdditionalReference = gridssAdditionalReference,
                jobIndex = i,
                assembleChunks = assembleChunks,
                tumorFinalBam = tumorFinalBam,
                normalFinalBam = normalFinalBam,
                normalSvBam = normalGridssPreprocess.svBam,
                normalCigarMetrics = normalGridssPreprocess.cigarMetrics,
                normalIdsvMetrics = normalGridssPreprocess.idsvMetrics,
                normalTagMetrics = normalGridssPreprocess.tagMetrics,
                normalMapqMetrics = normalGridssPreprocess.mapqMetrics,
                normalInsertSizeMetrics = normalGridssPreprocess.insertSizeMetrics,
                tumorSvBam = tumorGridssPreprocess.svBam,
                tumorCigarMetrics = tumorGridssPreprocess.cigarMetrics,
                tumorIdsvMetrics = tumorGridssPreprocess.idsvMetrics,
                tumorTagMetrics = tumorGridssPreprocess.tagMetrics,
                tumorMapqMetrics = tumorGridssPreprocess.mapqMetrics,
                tumorInsertSizeMetrics = tumorGridssPreprocess.insertSizeMetrics   
        }
    }
    
    call calling.GridssAssemble {
        input:
            threads = threads,
            memory_gb = 100,
            pairName = pairName,
            gridssReferenceFa = gridssReferenceFa,
            gridssAdditionalReference = gridssAdditionalReference,
            tumorFinalBam = tumorFinalBam,
            normalFinalBam = normalFinalBam,
            downsampled = GridssAssembleChunk.downsampled,
            excluded = GridssAssembleChunk.excluded,
            subsetCalled = GridssAssembleChunk.subsetCalled,
            
            normalSvBam = normalGridssPreprocess.svBam,
            normalCigarMetrics = normalGridssPreprocess.cigarMetrics,
            normalIdsvMetrics = normalGridssPreprocess.idsvMetrics,
            normalTagMetrics = normalGridssPreprocess.tagMetrics,
            normalMapqMetrics = normalGridssPreprocess.mapqMetrics,
            normalInsertSizeMetrics = normalGridssPreprocess.insertSizeMetrics,
            tumorSvBam = tumorGridssPreprocess.svBam,
            tumorCigarMetrics = tumorGridssPreprocess.cigarMetrics,
            tumorIdsvMetrics = tumorGridssPreprocess.idsvMetrics,
            tumorTagMetrics = tumorGridssPreprocess.tagMetrics,
            tumorMapqMetrics = tumorGridssPreprocess.mapqMetrics,
            tumorInsertSizeMetrics = tumorGridssPreprocess.insertSizeMetrics
            
    }
    
    call calling.GridssCalling {
        input:
            threads = threads,
            memory_gb = 100,
            pairName = pairName,
            gridssReferenceFa = gridssReferenceFa,
            gridssAdditionalReference = gridssAdditionalReference,
            tumorFinalBam = tumorFinalBam,
            normalFinalBam = normalFinalBam,
            downsampled = GridssAssembleChunk.downsampled,
            excluded = GridssAssembleChunk.excluded,
            subsetCalled = GridssAssembleChunk.subsetCalled,
            gridssassemblyBam = GridssAssemble.gridssassemblyBam,
            gridssassemblySvBam = GridssAssemble.gridssassemblySvBam,
            
            normalSvBam = normalGridssPreprocess.svBam,
            normalCigarMetrics = normalGridssPreprocess.cigarMetrics,
            normalIdsvMetrics = normalGridssPreprocess.idsvMetrics,
            normalTagMetrics = normalGridssPreprocess.tagMetrics,
            normalMapqMetrics = normalGridssPreprocess.mapqMetrics,
            normalInsertSizeMetrics = normalGridssPreprocess.insertSizeMetrics,
            tumorSvBam = tumorGridssPreprocess.svBam,
            tumorCigarMetrics = tumorGridssPreprocess.cigarMetrics,
            tumorIdsvMetrics = tumorGridssPreprocess.idsvMetrics,
            tumorTagMetrics = tumorGridssPreprocess.tagMetrics,
            tumorMapqMetrics = tumorGridssPreprocess.mapqMetrics,
            tumorInsertSizeMetrics = tumorGridssPreprocess.insertSizeMetrics
            
    }
    
    call calling.GridssFilter {
        input:
            threads = threads,
            memory_gb = memory_gb,
            pairName = pairName,
            bsGenome = bsGenome,
            ponTarGz = ponTarGz,
            gridssUnfilteredVcf = GridssCalling.gridssUnfilteredVcf
            
    }
    
    output {
        File gridssVcf = GridssFilter.gridssVcf
    }
}