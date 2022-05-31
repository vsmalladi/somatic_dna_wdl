version 1.0

import "./wdl_structs.wdl"
import "tasks/bam_cram_conversion.wdl" as cramConversion


workflow BamToCramWorkflow {
    input {
        IndexedReference referenceFa
        Array[pairInfo]+ pairInfos
        Array[SampleBamInfo]+ normalSampleBamInfos
    }

    # need to find UNIQUE bams (don't convert if part of more than one pair)
    call cramConversion.UniqueBams as uniqueBams {
        input:
            pairInfosJson = write_json(pairInfos)
    }

    scatter(bamInfo in uniqueBams.uniqueBams) {
        call cramConversion.SamtoolsBamToCram as bamToCram {
            input:
                inputBam = bamInfo.finalBam,
                referenceFa = referenceFa,
                sampleId = bamInfo.sampleId,
                diskSize = (ceil(size(bamInfo.finalBam.bam, "GB") * 1.5)) + 50
        }
    }

    # then create a new map of sample cram infos (to use as input where crams used instead of bams)
    # treating is as new object. Other option is to make a new struct with bam+cram together
    call cramConversion.UpdateCramInfos as updateCramInfo {
        input:
            pairInfosJson = write_json(pairInfos),
            normalInfosJson = write_json(normalSampleBamInfos),
            cramInfosJson = write_json(bamToCram.cramInfo)
    }
    #    # the output of this can be used like:
    #    Array normalSampleCramInfos = updateCramInfo.normalSampleCramInfos
    #    Array pairCramInfos = updateCramInfo.pairCramInfos

    output {
        Array[SampleCramInfo] crams = bamToCram.cramInfo

    }
}
