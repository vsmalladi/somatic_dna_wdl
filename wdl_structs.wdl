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
    File alt
    File amb    
}

struct IndexedReference {
    File fasta
    File dict?
    File index
}

struct Fastqs {
    File fastqR1
    String bamBase
    String rgInfo
    String? md5sumR1
    File fastqR2
    String? md5sumR1
}

struct IndexedVcf {
    File vcf
    File index
}

struct IndexedTable {
    File table
    File index
}

