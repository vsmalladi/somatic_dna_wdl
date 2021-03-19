version 1.0

import "merge_vcf.wdl" as mergeVcf
import "../calling/calling.wdl" as calling
import "../wdl_structs.wdl"

# note that we will still need to match dictionaries or provide matching reference files

workflow MergeCallers {
    input {
        String tumor
        String normal
        String pairName
        Array[String]+ listOfChroms
        Array[IndexedVcf]+ allVcfCompressed
        
        File intervalListBed
        
        IndexedReference referenceFa
        Bam normalFinalBam
        Bam tumorFinalBam
        
        File ponFile
        
        IndexedVcf germFile
    }
    
    scatter(vcfCompressed in allVcfCompressed) {
            File allVcfCompressedFile = vcfCompressed.vcf
        }
    Array[File]+ allVcfCompressedList = allVcfCompressedFile
    
    scatter(chrom in listOfChroms) {
        call mergeVcf.MergeCallers as allMergeCallers {
            input:
                chrom = chrom,
                pairName = pairName,
                allVcfCompressed = allVcfCompressed,
                allVcfCompressedList = allVcfCompressedList
        }
        
        call mergeVcf.StartCandidates { 
        input:
            chrom = chrom,
            pairName = pairName,
            intervalListBed = intervalListBed,
            mergedChromVcf = allMergeCallers.mergedChromVcf
            
        }
        
        call mergeVcf.GetCandidates {
        input:
            chrom = chrom,
            pairName = pairName,
            startChromVcf = StartCandidates.startChromVcf
            
        }
        
        call mergeVcf.VcfToBed {
        input:
            chrom = chrom,
            pairName = pairName,
            candidateChromVcf = GetCandidates.candidateChromVcf
            
        }
        
        call calling.LancetWGSRegional as lancetConfirm {
            input:
                chrom = chrom,
                chromBed = VcfToBed.candidateChromBed,
                referenceFa = referenceFa,
                normalFinalBam = normalFinalBam,
                tumorFinalBam = tumorFinalBam,
                pairName = pairName,
                lancetChromVcfPath = "~{pairName}.lancet.merged.v6.~{chrom}.vcf",
                threads = 8,
                memoryGb = 40,
                diskSize = (ceil( size(tumorFinalBam.bam, "GB") + size(normalFinalBam.bam, "GB")) ) + 20
        }
        
        # lancetConfirm.lancetChromVcf
        call mergeVcf.CompressVcf as lancetCompressVcf {
            input:
                vcf = lancetConfirm.lancetChromVcf,
                memoryGb = 4
        }
        
        call mergeVcf.IndexVcf as lancetIndexVcf {
            input:
                vcfCompressed = lancetCompressVcf.vcfCompressed
        }
        
        # GetCandidates.candidateChromVcf
        call mergeVcf.CompressVcf as candidateCompressVcf {
            input:
                vcf = GetCandidates.candidateChromVcf,
                memoryGb = 4
        }
        
        call mergeVcf.IndexVcf as candidateIndexVcf {
            input:
                vcfCompressed = candidateCompressVcf.vcfCompressed
        }
        
        call mergeVcf.IntersectVcfs { 
            input:
                chrom = chrom,
                pairName = pairName,
                vcfCompressedLancet = lancetIndexVcf.vcfCompressedIndexed,
                vcfCompressedCandidate = candidateIndexVcf.vcfCompressedIndexed
                
        }
        
        #  =================================================================
        #                   Prep supporting Lancet calls
        #  =================================================================
        
        call mergeVcf.RenameMetadata {
            input:
                pairName = pairName,
                callerVcf = IntersectVcfs.vcfConfirmedCandidate,
                tool = "lancet"
        }
        
        call mergeVcf.MergePrepSupport {
            input:
                pairName = pairName,
                tool = "lancet",
                renameMetaVcf = RenameMetadata.renameMetaVcf
        }
        
        call mergeVcf.RenameVcf {
            input:
                pairName = pairName,
                tumor = tumor,
                normal = normal,
                tool = "lancet",
                prepCallerVcf = MergePrepSupport.prepCallerVcf
        }
        
        call mergeVcf.CompressVcf as confirmedCompressVcf {
            input:
                vcf = RenameVcf.renameVcf,
                memoryGb = 4  
        }
        
        call mergeVcf.IndexVcf as confirmedIndexVcf {
            input:
                vcfCompressed = confirmedCompressVcf.vcfCompressed
        }
        
        call mergeVcf.SplitMultiAllelic as confirmSplitMultiAllelic {
            input:
                pairName = pairName,
                vcfCompressedIndexed = confirmedIndexVcf.vcfCompressedIndexed,
                splitVcfPath = sub(basename(confirmedIndexVcf.vcfCompressedIndexed.vcf), ".rename.vcf.gz$", ".split.vcf"),
                referenceFa = referenceFa
        }
        
        call mergeVcf.SplitMnv {
            input:
                tool = "lancet",
                mnvVcfPath = sub(basename(confirmedIndexVcf.vcfCompressedIndexed.vcf), ".rename.vcf.gz$", ".split.vcf"),
                splitVcf = confirmSplitMultiAllelic.splitVcf  
        }
        
        call mergeVcf.RemoveContig {
            input:
                mnvVcfPath = sub(basename(confirmedIndexVcf.vcfCompressedIndexed.vcf), ".rename.vcf.gz$", ".split.vcf"),
                removeChromVcf = SplitMnv.mnvVcf
        }
        
        call mergeVcf.Gatk4MergeSortVcf {
            input:
                tempVcfs = [RemoveContig.removeContigVcf],
                sortedVcfPath = sub(basename(select_first([RemoveContig.removeContigVcf, SplitMnv.mnvVcf])), "$", ".gz"),
                referenceFa = referenceFa,
                gzipped = true,
                memoryGb = 16,
                diskSize = 60
        }
        
        call mergeVcf.CompressVcf as callerCompressVcf {
            input:
                vcf = allMergeCallers.mergedChromVcf,
                memoryGb = 4
        }
    
        call mergeVcf.IndexVcf as callerIndexVcf {
            input:
                vcfCompressed = callerCompressVcf.vcfCompressed
        }
        
        call mergeVcf.MergeCallers as lancetMergeCallers {
            input:
                chrom = chrom,
                pairName = pairName,
                allVcfCompressed = [callerIndexVcf.vcfCompressedIndexed, Gatk4MergeSortVcf.sortedVcf],
                allVcfCompressedList = [callerIndexVcf.vcfCompressedIndexed.vcf, Gatk4MergeSortVcf.sortedVcf.vcf]
        }
        
        #  =================================================================
        #                     Merge columns
        #  =================================================================
        
        call mergeVcf.MergeColumns {
            input:
                chrom = chrom,
                tumor = tumor,
                normal = normal,
                pairName = pairName,
                supportedChromVcf = lancetMergeCallers.mergedChromVcf
        }
        
        call mergeVcf.AddNygcAlleleCountsToVcf {
            input:
                chrom = chrom,
                pairName = pairName,
                columnChromVcf = MergeColumns.columnChromVcf,
                normalFinalBam = normalFinalBam,
                tumorFinalBam = tumorFinalBam
        }
        
        call mergeVcf.AddFinalAlleleCountsToVcf {
            input:
                chrom = chrom,
                pairName = pairName,
                preCountsChromVcf = AddNygcAlleleCountsToVcf.preCountsChromVcf
        }
        
        call mergeVcf.FilterPon {
            input:
                chrom = chrom,
                pairName = pairName,
                countsChromVcf = AddFinalAlleleCountsToVcf.countsChromVcf,
                ponFile = ponFile
        }
        
        call mergeVcf.FilterVcf {
            input:
                chrom = chrom,
                pairName = pairName,
                ponOutFile = FilterPon.ponOutFile,
                germFile = germFile
        }
        
        call mergeVcf.SnvstomnvsCountsbasedfilterAnnotatehighconf {
            input:
                chrom = chrom,
                pairName = pairName,
                filteredOutFile = FilterVcf.filteredOutFile
        }
        
    }
    
    output {
            Array[File] finalChromVcf = SnvstomnvsCountsbasedfilterAnnotatehighconf.finalChromVcf
        }
}