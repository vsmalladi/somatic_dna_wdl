version 1.0


struct Bam {
    File bam
    File bamIndex
    String? md5sum
}

struct BwaReference {
    File fasta
    File sa
    File pac
    File bwt
    File ann
    File? alt
    File amb
    File? dict
    File? index
}

struct IndexedReference {
    File fasta
    File? dict
    File index
}

struct Fastqs {
    File fastqR1
    String? md5sumR1
    File fastqR2
    String? md5sumR2
    String sampleId
    String readgroupId
    String rgpu
}

struct sampleInfo {
    String sampleId
    Array[Fastqs] listOfFastqPairs
    Float expectedCoverage
}

struct SampleBamInfo {
    String sampleId
    Bam finalBam
}

struct PairRawVcfInfo {
    String pairId
    File? mergedVcf
    File? mainVcf
    File? supplementalVcf
    File filteredMantaSV
    File strelka2Snv
    File strelka2Indel
    File mutect2
    File lancet
    File svabaSv
    File svabaIndel
    File gridssVcf
    File bicseq2Png
    File bicseq2
    String tumor
    String normal
    Bam tumorFinalBam
    Bam normalFinalBam
}

struct MergedPairVcfInfo {
    String pairId
    String tumor
    String normal
    File unannotatedVcf
}

struct PairVcfInfo {
    String pairId
    String tumor
    String normal
    File mainVcf
    File supplementalVcf
    File vcfAnnotatedTxt
    File maf
}

struct FinalPairInfo {
    String pairId
    String tumor
    String normal
    File mainVcf
    File supplementalVcf
    File filteredMantaSV
    File strelka2Snv
    File strelka2Indel
    File mutect2
    File lancet
    File svabaSv
    File svabaIndel
    File gridssVcf
    File bicseq2Png
    File bicseq2
    Bam tumorFinalBam
    Bam normalFinalBam
}

struct PairRelationship {
    String pairId
    String tumor
    String normal
}

struct pairInfo {
    String pairId
    Bam tumorFinalBam
    Bam normalFinalBam
    String tumor
    String normal
}

struct IndexedVcf {
    File vcf
    File index
}

struct IndexedTable {
    File table
    File index
}
