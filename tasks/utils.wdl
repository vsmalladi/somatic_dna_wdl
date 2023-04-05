version 1.0

import "../wdl_structs.wdl"

task CreateBlankFile {
    input {
        String fileId
        String blankFilePath = "blank_file_~{fileId}.txt"
        Int diskSize = 1
        Int memoryGb = 1
    }

    command {
        echo "" > ~{blankFilePath}
    }

    output {
        File blankFile = "~{blankFilePath}"
    }

    runtime {
        mem: memoryGb + "G"
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-public/genome-utils@sha256:59603ab0aeda571c38811d6d1820d1b546a69fc342120bef75597bfd7905ea1f"
    }

}

# for wdl version 1.0

task GetIndex {
    input {
        String sampleId
        Array[String] sampleIds
    }

    command {
        python /get_index.py \
        --sample-id ~{sampleId} \
        --sample-ids ~{sep=' ' sampleIds}
    }

    output {
        Int index = read_int(stdout())
    }

    runtime {
        docker: "gcr.io/nygc-public/workflow_utils@sha256:40fa18ac3f9d9f3b9f037ec091cb0c2c26ad6c7cb5c32fb16c1c0cf2a5c9caea"
    }
}
