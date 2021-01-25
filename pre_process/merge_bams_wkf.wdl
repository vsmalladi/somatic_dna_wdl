version 1.0

import "merge_bams.wdl" as mergeBams
import "../wdl_structs.wdl"

workflow MergeBams {
    # command 
    #     merge lane level BAMs
    input {
        #    command merge flowcell
        Array[Bam] laneFixmateBams
        Array[String] laneFixmateBamPaths
        String sampleId
        IndexedVcf MillsAnd1000G
        IndexedVcf Indels
        IndexedVcf DbSnp
        File chromFile
        IndexedReference indexedReference
        # resources
        Int mem
        Int threads
        String gatkDockerImage
        String novosortDockerImage
        String samtoolsDockerImage
    }
    
    call mergeBams.NovosortMarkDup {
        input:
            laneBams = laneFixmateBams,
            laneBamsLists = laneFixmateBamPaths,
            sampleId = sampleId,
            dockerImage = novosortDockerImage,
            mem = mem,
            threads = threads
    }
    
    call mergeBams.IndexBam as novosortMarkDupIndexed {
        input:
            bam = NovosortMarkDup.mergedDedupBamOnly,
            dockerImage = samtoolsDockerImage
    }
    
    call mergeBams.Bqsr38 {
        input:
            mergedDedupBam = novosortMarkDupIndexed.indexedBam,
            MillsAnd1000G = MillsAnd1000G,
            indexedReference = indexedReference,
            Indels = Indels,
            DbSnp = DbSnp,
            chromFile = chromFile,
            sampleId = sampleId,
            mem = mem,
            threads = threads,
            dockerImage = gatkDockerImage
            
    }
    call mergeBams.PrintReads {
        input:
            indexedReference = indexedReference,
            mergedDedupBam = novosortMarkDupIndexed.indexedBam,
            recalGrp = Bqsr38.recalGrp,
            sampleId = sampleId,
            mem = mem,
            threads = threads,
            dockerImage = gatkDockerImage
            
    }
    
    output {
        Bam mergedDedupBam = novosortMarkDupIndexed.indexedBam
        Bam finalBam = PrintReads.finalBam
    }
}