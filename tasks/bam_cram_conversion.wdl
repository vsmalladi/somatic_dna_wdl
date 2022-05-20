version 1.0

import "../wdl_structs.wdl"

task SamtoolsBamToCram {
    input {
        # command
        Bam inputBam
        IndexedReference referenceFa
        String sampleId
        String cramPath = "~{sampleId}.final.cram"
        String indexPath = "~{sampleId}.final.cram.crai"
        # resources
        Int diskSize
        Int threads = 8
        Int memoryGb = 8
    }

    command {
        set -e -o pipefail

        samtools \
        view \
        -C \
        -T ~{referenceFa.fasta} \
        -t ~{referenceFa.index} \
        -o ~{cramPath} ~{inputBam.bam} \
        --verbosity=8 \
        --threads ~{threads}

        samtools index ~{cramPath} ~{indexPath}
    }

    output {
        Cram finalCram = object {
                             cram : "~{cramPath}",
                             cramIndex : "~{indexPath}"
                         }
    }

    runtime {
        docker : "gcr.io/nygc-public/samtools@sha256:32f29fcd7af01b3941e6f93095e8d899741e81b50bcc838329bd8df43e120cc3"
        disks: "local-disk " + diskSize + " LOCAL"
        memory: memoryGb + " GB"
        cpu: threads
    }
}

task SamtoolsCramToBam {
    input {
        # command
        Cram inputCram
        IndexedReference referenceFa
        String sampleId
        String bamPath = "~{sampleId}.final.bam"
        String indexPath = "~{sampleId}.final.bai"
        # resources
        Int diskSize
        Int threads = 8
        Int memoryGb = 8
    }

    command {
        set -e -o pipefail

        samtools \
        view \
        -b \
        -T ~{referenceFa.fasta} \
        -t ~{referenceFa.index} \
        -o ~{bamPath} ~{inputCram.cram} \
        --verbosity=8 \
        --threads ~{threads}

        samtools index ~{bamPath} ~{indexPath}
    }

    output {
        Bam finalBam = object {
                           bam : "~{bamPath}",
                           bamIndex : "~{indexPath}"
                       }
    }

    runtime {
        docker : "gcr.io/nygc-public/samtools@sha256:32f29fcd7af01b3941e6f93095e8d899741e81b50bcc838329bd8df43e120cc3"
        disks: "local-disk " + diskSize + " LOCAL"
        memory: memoryGb + " GB"
        cpu: threads
    }
}


task UniqueBams {
    input {
        # command
        String pairInfosPath
    }

    command {
        python /unique_bams.py \
        --pair-infos ~{pairInfosPath}
    }

    output {
        Array[Bam] uniqueBams = read_json('unique_bams.json')
    }

    runtime {
        docker : "TBD"
    }
}

