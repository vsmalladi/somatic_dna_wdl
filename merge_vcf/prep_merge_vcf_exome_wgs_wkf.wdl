version 1.0

import "merge_vcf.wdl" as merge_vcf
import "../wdl_structs.wdl"

workflow PrepMergeVcf {
    input {
        String tumorId
        String normalId
        String tool
        File callerVcf
        String pairName
        IndexedReference referenceFa
    }
    
    call merge_vcf.RenameExomeWgsMetadata {
        input:
            callerVcf = callerVcf,
            tool=tool,
            pairName=pairName
    }
    
    call merge_vcf.MergePrepExomeWgs {
        input:
            pairName = pairName,
            tool = tool,
            renameMetaVcf = RenameExomeWgsMetadata.renameMetaVcf
    }

    call merge_vcf.RenameVcf {
        input:
            pairName = pairName,
            tumorId = tumorId,
            normalId = normalId,
            tool = tool,
            prepCallerVcf = MergePrepExomeWgs.prepCallerVcf
    }
    
    call merge_vcf.CompressVcf as renameCompressVcf {
        input:
            vcf = RenameVcf.renameVcf,
            memoryGb = 4
    }
    
    call merge_vcf.IndexVcf as renameIndexVcf {
        input:
            vcfCompressed = renameCompressVcf.vcfCompressed
    }
    
    call merge_vcf.SplitMultiAllelic as prepSplitMultiAllelic {
        input:
            pairName = pairName,
            referenceFa = referenceFa,
            vcfCompressedIndexed = renameIndexVcf.vcfCompressedIndexed,
            splitVcfPath = sub(basename(renameIndexVcf.vcfCompressedIndexed.vcf), ".rename.vcf.gz$", ".split.vcf")
    }
    
    call merge_vcf.SplitMnv {
        input:
            tool = tool,
            mnvVcfPath = sub(basename(renameIndexVcf.vcfCompressedIndexed.vcf), ".rename.vcf.gz$", ".split.vcf"),
            splitVcf = prepSplitMultiAllelic.splitVcf
            
    }
    
    call merge_vcf.Gatk4MergeSortVcf {
        input:
            tempVcfs = [SplitMnv.mnvVcf],
            sortedVcfPath = sub(basename(SplitMnv.mnvVcf), "$", ".gz"),
            referenceFa = referenceFa,
            gzipped = true,
            threads = 4,
            memoryGb = 8,
            diskSize = 20
            
    }
    
    output {
        IndexedVcf preppedVcf = Gatk4MergeSortVcf.sortedVcf
    }
}
