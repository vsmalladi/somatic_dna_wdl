version 1.0

import "merge_vcf.wdl" as merge_vcf
import "../wdl_structs.wdl"

# note that we will still need to match dictionaries or provide matching reference files

workflow MergeCallers {
    input {
        String tumor
        String normal
        String pairName
        Array[String]+ listOfChroms
        Array[IndexedVcf]+ allVcfCompressed
        
        File knownGeneBed
        
        IndexedReference referenceFa
        Bam normalFinalBam
        Bam tumorFinalBam
        
        File ponFile
        
        IndexedVcf germFile
        
        String bedtoolsDockerImage
        String bcftoolsDockerImage
        String pysamDockerImage
        String lancetDockerImage
        String gatkDockerImage
        String bgzipDockerImage
        Int threads
        Int memoryGb
    }
    
    scatter(vcfCompressed in allVcfCompressed) {
            File allVcfCompressedFile = vcfCompressed.vcf
        }
    Array[File]+ allVcfCompressedList = allVcfCompressedFile
    
    scatter(chrom in listOfChroms) {
        call merge_vcf.MergeCallers as allMergeCallers {
            input:
                chrom = chrom,
                pairName = pairName,
                allVcfCompressed = allVcfCompressed,
                allVcfCompressedList = allVcfCompressedList,
                memoryGb = memoryGb,
                threads = threads,
                dockerImage = bcftoolsDockerImage
        }
        
        call merge_vcf.StartCandidates { 
        input:
            chrom = chrom,
            pairName = pairName,
            knownGeneBed = knownGeneBed,
            mergedChromVcf = allMergeCallers.mergedChromVcf,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = bedtoolsDockerImage
            
        }
        
        call merge_vcf.GetCandidates {
        input:
            chrom = chrom,
            pairName = pairName,
            startChromVcf = StartCandidates.startChromVcf,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = pysamDockerImage
            
        }
        
        call merge_vcf.VcfToBed {
        input:
            chrom = chrom,
            pairName = pairName,
            candidateChromVcf = GetCandidates.candidateChromVcf,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = bedtoolsDockerImage
            
        }
        
        call merge_vcf.LancetConfirm {
            input:
                chrom = chrom,
                pairName = pairName,
                referenceFa = referenceFa,
                normalFinalBam = normalFinalBam,
                tumorFinalBam = tumorFinalBam,
                candidateChromBed = VcfToBed.candidateChromBed,
                memoryGb = memoryGb,
                threads = threads,
                dockerImage = lancetDockerImage  
        }
        
        # LancetConfirm.lancetChromVcf
        call merge_vcf.CompressVcf as lancetCompressVcf {
            input:
                vcf = LancetConfirm.lancetChromVcf,
                memoryGb = memoryGb,
                threads = threads,
                dockerImage = bgzipDockerImage
        }
        
        call merge_vcf.IndexVcf as lancetIndexVcf {
            input:
                vcfCompressed = lancetCompressVcf.vcfCompressed,
                memoryGb = memoryGb,
                threads = threads,
                dockerImage = gatkDockerImage
        }
        
        # GetCandidates.candidateChromVcf
        call merge_vcf.CompressVcf as candidateCompressVcf {
            input:
                vcf = GetCandidates.candidateChromVcf,
                memoryGb = memoryGb,
                threads = threads,
                dockerImage = bgzipDockerImage
        }
        
        call merge_vcf.IndexVcf as candidateIndexVcf {
            input:
                vcfCompressed = candidateCompressVcf.vcfCompressed,
                memoryGb = memoryGb,
                threads = threads,
                dockerImage = gatkDockerImage
        }
        
        call merge_vcf.IntersectVcfs { 
            input:
                chrom = chrom,
                pairName = pairName,
                vcfCompressedLancet = lancetIndexVcf.vcfCompressedIndexed,
                vcfCompressedCandidate = candidateIndexVcf.vcfCompressedIndexed,
                memoryGb = memoryGb,
                threads = threads,
                dockerImage = bcftoolsDockerImage
                
        }
        
        #  =================================================================
        #                   Prep supporting Lancet calls
        #  =================================================================
        
        call merge_vcf.RenameMetadata {
            input:
                callerVcf = IntersectVcfs.vcfConfirmedCandidate,
                memoryGb = memoryGb,
                threads = threads,
                dockerImage = pysamDockerImage
        }
        
        call merge_vcf.MergePrepSupport {
            input:
                tool = "lancet",
                renameMetaVcf = RenameMetadata.renameMetaVcf,
                memoryGb = memoryGb,
                threads = threads,
                dockerImage = pysamDockerImage
        }
        
        call merge_vcf.RenameVcf {
            input:
                pairName = pairName,
                tumor = tumor,
                normal = normal,
                tool = "lancet",
                prepCallerVcf = MergePrepSupport.prepCallerVcf,
                memoryGb = memoryGb,
                threads = threads,
                dockerImage = pysamDockerImage
        }
        
        call merge_vcf.CompressVcf as confirmedCompressVcf {
            input:
                vcf = RenameVcf.renameVcf,
                memoryGb = memoryGb,
                threads = threads,
                dockerImage = bgzipDockerImage
        }
        
        call merge_vcf.IndexVcf as confirmedIndexVcf {
            input:
                vcfCompressed = confirmedCompressVcf.vcfCompressed,
                memoryGb = memoryGb,
                threads = threads,
                dockerImage = gatkDockerImage
        }
        
        call merge_vcf.SplitMultiAllelic as confirmSplitMultiAllelic {
            input:
                vcfCompressedIndexed = confirmedIndexVcf.vcfCompressedIndexed,
                splitVcfPath = sub(basename(confirmedIndexVcf.vcfCompressedIndexed.vcf), ".rename.vcf.gz$", ".split.vcf"),
                memoryGb = memoryGb,
                threads = threads,
                dockerImage = bcftoolsDockerImage
        }
        
        call merge_vcf.SplitMnv {
            input:
                tool = "lancet",
                mnvVcfPath = sub(basename(confirmedIndexVcf.vcfCompressedIndexed.vcf), ".rename.vcf.gz$", ".split.vcf"),
                splitVcf = confirmSplitMultiAllelic.splitVcf,
                memoryGb = memoryGb,
                threads = threads,
                dockerImage = pysamDockerImage     
        }
        
        call merge_vcf.RemoveContig {
            input:
                memoryGb = memoryGb,
                threads = threads,
                dockerImage = pysamDockerImage, 
                mnvVcfPath = sub(basename(confirmedIndexVcf.vcfCompressedIndexed.vcf), ".rename.vcf.gz$", ".split.vcf"),
                removeChromVcf = SplitMnv.mnvVcf
        }
        
        call merge_vcf.Gatk4MergeSortVcf {
            input:
                tempVcfs = [RemoveContig.removeContigVcf],
                sortedVcfPath = sub(basename(select_first([RemoveContig.removeContigVcf, SplitMnv.mnvVcf])), "$", ".gz"),
                referenceFa = referenceFa,
                memoryGb = memoryGb,
                threads = threads,
                dockerImage = gatkDockerImage 
        }
        
        call merge_vcf.CompressVcf as callerCompressVcf {
            input:
                vcf = allMergeCallers.mergedChromVcf,
                memoryGb = memoryGb,
                threads = threads,
                dockerImage = bgzipDockerImage
        }
    
        call merge_vcf.IndexVcf as callerIndexVcf {
            input:
                vcfCompressed = callerCompressVcf.vcfCompressed,
                memoryGb = memoryGb,
                threads = threads,
                dockerImage = gatkDockerImage
        }
        
        call merge_vcf.MergeCallers as lancetMergeCallers {
            input:
                chrom = chrom,
                pairName = pairName,
                allVcfCompressed = [callerIndexVcf.vcfCompressedIndexed, Gatk4MergeSortVcf.sortedVcf],
                allVcfCompressedList = [callerIndexVcf.vcfCompressedIndexed.vcf, Gatk4MergeSortVcf.sortedVcf.vcf],
                memoryGb = memoryGb,
                threads = threads,
                dockerImage = bcftoolsDockerImage
        }
        
        #  =================================================================
        #                     Merge columns
        #  =================================================================
        
        call merge_vcf.MergeColumns {
            input:
                chrom = chrom,
                tumor = tumor,
                normal = normal,
                pairName = pairName,
                supportedChromVcf = lancetMergeCallers.mergedChromVcf,
                memoryGb = memoryGb,
                threads = threads,
                dockerImage = pysamDockerImage
        }
        
        call merge_vcf.AddNygcAlleleCountsToVcf {
            input:
                chrom = chrom,
                pairName = pairName,
                columnChromVcf = MergeColumns.columnChromVcf,
                normalFinalBam = normalFinalBam,
                tumorFinalBam = tumorFinalBam,
                memoryGb = memoryGb,
                threads = threads,
                dockerImage = pysamDockerImage
        }
        
        call merge_vcf.AddFinalAlleleCountsToVcf {
            input:
                chrom = chrom,
                pairName = pairName,
                preCountsChromVcf = AddNygcAlleleCountsToVcf.preCountsChromVcf,
                memoryGb = memoryGb,
                threads = threads,
                dockerImage = pysamDockerImage
        }
        
        call merge_vcf.FilterPon {
            input:
                chrom = chrom,
                pairName = pairName,
                countsChromVcf = AddFinalAlleleCountsToVcf.countsChromVcf,
                ponFile = ponFile,
                memoryGb = memoryGb,
                threads = threads,
                dockerImage = pysamDockerImage
        }
        
        call merge_vcf.FilterVcf {
            input:
                chrom = chrom,
                pairName = pairName,
                ponOutFile = FilterPon.ponOutFile,
                germFile = germFile,
                memoryGb = memoryGb,
                threads = threads,
                dockerImage = pysamDockerImage
        }
        
        call merge_vcf.SnvstomnvsCountsbasedfilterAnnotatehighconf {
            input:
                chrom = chrom,
                pairName = pairName,
                filteredOutFile = FilterVcf.filteredOutFile,
                memoryGb = memoryGb,
                threads = threads,
                dockerImage = pysamDockerImage
        }
        
    }
    
    output {
            Array[File] finalChromVcf = SnvstomnvsCountsbasedfilterAnnotatehighconf.finalChromVcf
        }
}