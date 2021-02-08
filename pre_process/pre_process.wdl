version 1.0

import "align_fastq_wkf.wdl" as alignFastq
import "merge_bams_wkf.wdl" as mergeBams
import "../wdl_structs.wdl"

workflow Preprocess {
    # command 
    #   align FASTQ files
    #   and 
    #   merge lane level BAMs
    input {
        Array[Fastqs]+ listOfFastqPairs
        BwaReference bwaReference
        #    command merge flowcell
        String sampleId
        IndexedVcf MillsAnd1000G
        IndexedVcf Indels
        IndexedVcf DbSnp
        File chromFile
        IndexedReference indexedReference
        # resources
        #    prep flowcell
        Int mem
        Int threads
        String bwaDockerImage
        String shortAlignDockerImage
        String gatkDockerImage
        #    merge flowcell
        String novosortDockerImage
        String samtoolsDockerImage

    }
    
    call alignFastq.AlignFastq {
        input:
            listOfFastqPairs = listOfFastqPairs,
            bwaReference = bwaReference,
            bwaMem = mem,
            threads = threads,
            bwaDockerImage = bwaDockerImage,
            shortAlignDockerImage = shortAlignDockerImage,
            gatkDockerImage = gatkDockerImage,
    }
    
    call mergeBams.MergeBams {
        input:
            laneFixmateBams = AlignFastq.laneFixmateBam,
            sampleId = sampleId,
            MillsAnd1000G = MillsAnd1000G,
            Indels = Indels,
            DbSnp = DbSnp,
            chromFile = chromFile,
            indexedReference = indexedReference,
            mem = mem,
            threads = threads,
            gatkDockerImage = gatkDockerImage,
            novosortDockerImage = novosortDockerImage,
            samtoolsDockerImage = samtoolsDockerImage
    }
    
    output {
        Bam mergedDedupBam = MergeBams.mergedDedupBam
        Bam finalBam = MergeBams.finalBam
    }
}