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
    
    Int lowPreMemoryGb = 16
    Int lowTumorDiskSize = ceil( size(tumorFinalBam.bam, "GB") * 1.5 )
    Int lowNormalDiskSize = ceil( size(normalFinalBam.bam, "GB") * 1.5 )
    
    if (highMem) {
        Int highPreMemoryGb = 32
        Int highTumorDiskSize = ceil( size(tumorFinalBam.bam, "GB") * 2 ) + 20
        Int highNormalDiskSize = ceil( size(normalFinalBam.bam, "GB") * 2 ) + 20
    }
    
    Int preMemoryGb = select_first([highPreMemoryGb, lowPreMemoryGb])
    Int tumorDiskSize = select_first([highTumorDiskSize, lowTumorDiskSize])
    Int normalDiskSize = select_first([highNormalDiskSize, lowNormalDiskSize])
    
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

    Int lowAssembleMemoryGb = 48
    Int lowAssembleDiskSize = ceil( size(tumorFinalBam.bam, "GB") * 1.4 ) + ceil( size(normalFinalBam.bam, "GB")  * 1.4)
    
    if (highMem) {
        Int highAssembleMemoryGb = 100
        Int highAssembleDiskSize = ceil( size(tumorFinalBam.bam, "GB") * 2 ) + ceil( size(normalFinalBam.bam, "GB")  * 2) + 20
    }
    
    Int assembleMemoryGb = select_first([highAssembleMemoryGb, lowAssembleMemoryGb])
    Int assembleDiskSize = select_first([highAssembleDiskSize, lowAssembleDiskSize])
   
   
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
                diskSize = assembleDiskSize
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
            diskSize = assembleDiskSize

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
            diskSize = assembleDiskSize

    }
    
    Int lowFilterDiskSize = 4
    
    if (highMem) {
        Int highFilterDiskSize = 30
    }
    
    Int filterDiskSize = select_first([highFilterDiskSize, lowFilterDiskSize])

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
