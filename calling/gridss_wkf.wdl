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

        BwaReference bwaReference
        Array[File] gridssAdditionalReference

        Bam normalFinalBam
        Bam tumorFinalBam

        Int assembleChunks = 5

        String bsGenome
        File ponTarGz

        Int threads = 8
        
        Boolean highMem = false
    }
    
    Int preMemoryGb = 16
    Int tumorDiskSize = ceil( size(tumorFinalBam.bam, "GB") * 1.5 )
    Int normalDiskSize = ceil( size(normalFinalBam.bam, "GB") * 1.5 )
    
    if (highMem) {
        Int preMemoryGb = 32
        Int tumorDiskSize = ceil( size(tumorFinalBam.bam, "GB") * 2 ) + 20
        Int normalDiskSize = ceil( size(normalFinalBam.bam, "GB") * 2 ) + 20
    }
    
    call calling.GridssPreprocess as tumorGridssPreprocess {
        input:
            threads = threads,
            memoryGb = preMemoryGb,
            bwaReference = bwaReference,
            gridssAdditionalReference = gridssAdditionalReference,
            finalBam = tumorFinalBam,
            diskSize = tumorDiskSize
    }

    call calling.GridssPreprocess as normalGridssPreprocess {
        input:
            threads = threads,
            memoryGb = preMemoryGb,
            bwaReference= bwaReference,
            gridssAdditionalReference = gridssAdditionalReference,
            finalBam = normalFinalBam,
            diskSize = normalDiskSize
    }

    Int assembleMemoryGb = 48
    Int assemblediskSize = ceil( size(tumorFinalBam.bam, "GB") * 1.4 ) + ceil( size(normalFinalBam.bam, "GB")  * 1.4)
    
    if (highMem) {
        Int assembleMemoryGb = 100
        Int assemblediskSize = ceil( size(tumorFinalBam.bam, "GB") * 2 ) + ceil( size(normalFinalBam.bam, "GB")  * 2) + 20
    }
    
    scatter(i in range(assembleChunks)) {
        call calling.GridssAssembleChunk {
            input:
                threads = threads,
                memoryGb = assembleMemoryGb,
                pairName = pairName,
                bwaReference = bwaReference,
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
                tumorInsertSizeMetrics = tumorGridssPreprocess.insertSizeMetrics,
                diskSize = assemblediskSize
        }
    }

    call calling.GridssAssemble {
        input:
            threads = threads,
            memoryGb = assembleMemoryGb,
            pairName = pairName,
            bwaReference = bwaReference,
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
            tumorInsertSizeMetrics = tumorGridssPreprocess.insertSizeMetrics,
            diskSize = assemblediskSize

    }

    call calling.GridssCalling {
        input:
            threads = threads,
            memoryGb = assembleMemoryGb,
            pairName = pairName,
            bwaReference = bwaReference,
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
            tumorInsertSizeMetrics = tumorGridssPreprocess.insertSizeMetrics,
            diskSize = assemblediskSize

    }
    
    Int filterDiskSize = 4
    
    if (highMem) {
        Int filterDiskSize = 30
    }

    call calling.GridssFilter {
        input:
            threads = threads,
            memoryGb = preMemoryGb,
            pairName = pairName,
            bsGenome = bsGenome,
            ponTarGz = ponTarGz,
            gridssUnfilteredVcf = GridssCalling.gridssUnfilteredVcf,
            diskSize = filterDiskSize

    }

    output {
        IndexedVcf gridssVcf = GridssFilter.gridssVcf
    }
}
