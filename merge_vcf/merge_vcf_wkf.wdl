version 1.0

import "prep_merge_vcf_wkf.wdl" as prepMergeVcf
import "merge_callers_wkf.wdl" as mergeCallers
import "merge_chroms_wkf.wdl" as mergeChroms
import "../wdl_structs.wdl"

workflow MergeVcf {
    # command
    #   Call variants in BAMs
    #   merge and filter raw VCFs
    #   annotate
    input {
        PairRawVcfInfo pairRawVcfInfo
        IndexedReference referenceFa
        Array[String]+ listOfChroms
        
        # merge callers
        File knownGeneBed
        File ponFile
        IndexedVcf germFile
    }
    
    call prepMergeVcf.PrepMergeVcf as filteredMantaSVPrepMergeVcf {
        input:
            callerVcf=pairRawVcfInfo.filteredMantaSV,
            tumor=pairRawVcfInfo.tumor,
            normal=pairRawVcfInfo.normal,
            tool='manta',
            pairName=pairRawVcfInfo.pairId,
            referenceFa=referenceFa
            
    }
    
    call prepMergeVcf.PrepMergeVcf as strelka2SnvPrepMergeVcf {
        input:
            callerVcf=pairRawVcfInfo.strelka2Snv,
            tumor=pairRawVcfInfo.tumor,
            normal=pairRawVcfInfo.normal,
            tool='strelka2',
            pairName=pairRawVcfInfo.pairId,
            referenceFa=referenceFa
            
    }
    
    call prepMergeVcf.PrepMergeVcf as strelka2IndelPrepMergeVcf {
        input:
            callerVcf=pairRawVcfInfo.strelka2Indel,
            tumor=pairRawVcfInfo.tumor,
            normal=pairRawVcfInfo.normal,
            tool='strelka2',
            pairName=pairRawVcfInfo.pairId,
            referenceFa=referenceFa
            
    }
    
    call prepMergeVcf.PrepMergeVcf as mutect2PrepMergeVcf {
        input:
            callerVcf=pairRawVcfInfo.mutect2,
            tumor=pairRawVcfInfo.tumor,
            normal=pairRawVcfInfo.normal,
            tool='mutect2',
            pairName=pairRawVcfInfo.pairId,
            referenceFa=referenceFa
            
    }
    
    call prepMergeVcf.PrepMergeVcf as lancetPrepMergeVcf {
        input:
            callerVcf=pairRawVcfInfo.lancet,
            tumor=pairRawVcfInfo.tumor,
            normal=pairRawVcfInfo.normal,
            tool='lancet',
            pairName=pairRawVcfInfo.pairId,
            referenceFa=referenceFa
            
    }
    
    call prepMergeVcf.PrepMergeVcf as svabaIndelPrepMergeVcf {
        input:
            callerVcf=pairRawVcfInfo.svabaIndel,
            tumor=pairRawVcfInfo.tumor,
            normal=pairRawVcfInfo.normal,
            tool='svaba',
            pairName=pairRawVcfInfo.pairId,
            referenceFa=referenceFa
            
    }
    
    call mergeCallers.MergeCallers {
        input:
            tumor=pairRawVcfInfo.tumor,
            normal=pairRawVcfInfo.normal,
            pairName=pairRawVcfInfo.pairId,
            listOfChroms=listOfChroms,
            knownGeneBed=knownGeneBed,
            referenceFa=referenceFa,
            normalFinalBam=pairRawVcfInfo.normalFinalBam,
            tumorFinalBam=pairRawVcfInfo.tumorFinalBam,
            ponFile=ponFile,
            germFile=germFile,
            allVcfCompressed=[filteredMantaSVPrepMergeVcf.preppedVcf, 
                strelka2SnvPrepMergeVcf.preppedVcf,
                strelka2IndelPrepMergeVcf.preppedVcf,
                mutect2PrepMergeVcf.preppedVcf,
                lancetPrepMergeVcf.preppedVcf,
                svabaIndelPrepMergeVcf.preppedVcf]
        
    }
    
    call mergeChroms.MergeChroms {
        input:
            tumor=pairRawVcfInfo.tumor,
            normal=pairRawVcfInfo.normal,
            pairName=pairRawVcfInfo.pairId,
            referenceFa=referenceFa,
            finalChromVcf=MergeCallers.finalChromVcf,
            
    }
    
    output {
        File mergedVcf = MergeChroms.unannotatedVcf
    }
}