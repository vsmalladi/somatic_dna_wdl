#!/bin/bash

# USAGE: run.sh [-h] --options OPTIONS --wdl-file WDL_FILE
#               --url URL --log-dir LOG_DIR
#               --project PROJECT 
#               [--library {WGS,Exome}]
#               [--genome {Human_GRCh38_full_analysis_set_plus_decoy_hla, Human_GRCh38_tcga}]
#               [--read-length READ_LENGTH]
#               [--pairs-file PAIRS_FILE]
#               [--samples-file SAMPLES_FILE]
#               [--interval-list {SureSelect_V6plusCOSMIC.target.GRCh38_full_analysis_set_plus_decoy_hla}]
#               [--custom-inputs [CUSTOM_INPUTS [CUSTOM_INPUTS ...]]]
#               [--skip-validate]
#               [--dry-run]
# DESCRIPTION: validate workflow, create input json and submit workflow to cromwell.
# Script requires jq, cromwell-tools, gcloud to be in the path
# Script returns the workflow uuid.
# Script shows submission command and command to check status .
# in the STDERR stream

# ================== COPYRIGHT ================================================
# New York Genome Center
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2021) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
#
#    Jennifer M Shelton (jshelton@nygenome.org)
#    Nico Robine (nrobine@nygenome.org)
#    Minita Shah (mshah@nygenome.org)
#    Timothy Chu (tchu@nygenome.org)
#    Will Hooper (whooper@nygenome.org)
#
# ================== /COPYRIGHT ===============================================

help_top="run.sh [-h] --options OPTIONS --wdl-file WDL_FILE
               --url URL --log-dir LOG_DIR
               --project PROJECT 
               [--library {WGS,Exome}]
               [--genome {Human_GRCh38_full_analysis_set_plus_decoy_hla, Human_GRCh38_tcga}]
               [--read-length READ_LENGTH]
               [--pairs-file PAIRS_FILE]
               [--samples-file SAMPLES_FILE]
               [--interval-list {SureSelect_V6plusCOSMIC.target.GRCh38_full_analysis_set_plus_decoy_hla}]
               [--custom-inputs [CUSTOM_INPUTS [CUSTOM_INPUTS ...]]]
               [--skip-validate]
               [--dry-run]
"

help_long="-h, --help            show this help message and exit
  --url URL             Cromwell server URL (required)
  --log-dir LOG_DIR     Output directory for all logs and reports
                        related to this workflow UUID (required)
  --options OPTIONS     Options json file for cromwell (required)
  --wdl-file WDL_FILE   WDL workflow. An input JSON that matches this
                        WDL workflow will be created (required)
  --library {WGS,Exome}
                        Sequence library type.
  --genome {Human_GRCh38_full_analysis_set_plus_decoy_hla, Human_GRCh38_tcga}
                        Genome key to use for pipeline.
  --project PROJECT     Project name associated with account.
  --read-length READ_LENGTH     Required only for steps like BiqSeq2 that 
                        use read_length as input.
  --pairs-file PAIRS_FILE
                        CSV file with items that are required to have
                        \"tumor\", \"normal\" and \"pair_id\" in the columns.
  --samples-file [SAMPLES_FILE]
                        Not generally required. If tasks only require
                        sample_id and do not use pairing information sample
                        info can be populated with a CSV file. The CSV file
                        requires a columns named [\"sampleId\"].
  --interval-list {SureSelect_V6plusCOSMIC.target.GRCh38_full_analysis_set_plus_decoy_hla}
                        File basename for interval list.If not supplied the
                        default (the SureSelect interval list for your genome)
                        will be used (only needed for future Exome workflows)
  --custom-inputs [CUSTOM_INPUTS]
                        Optional JSON file with custom input variables. The
                        name of the variable in the input file must match the
                        name of the variable in the WDL workflow. It is not
                        required that the input specify the workflow. By
                        default the input will be added to the top-level
                        workflow. Any variable defined in in this JSON will
                        overwrite any reference variable in the the config 
                        directory or workflow default.
  --skip-validate       Skip the step where input files are validated.
                        Otherwise all gs//: URIs will be checked to see that a
                        file exists. Disable with caution.Cromwell will launch
                        instances and run without checking. Test a small pairs
                        file to ensure all references exist and at least some
                        sample input files can be read by the current user.
  --dry-run             Skip the step where the job is submitted to cromwell-tools.
"

print_help() {
    echo "${help_top}"
    echo "DESCRIPTION: validate workflow, create input json and submit workflow to cromwell."
    echo "Script requires jq, cromwell-tools, gcloud to be in the path."
    echo "Script shows submission command and command to check status"
    echo "in the STDERR stream."
    echo "${help_long}"
    exit 1
}


# Defaults
library="WGS"
genome="Human_GRCh38_full_analysis_set_plus_decoy_hla"

for arg in "$@"; do
    case $arg in
        -u|--url)
        url="$2"
        shift # Remove argument name from processing
        shift # Remove argument value from processing
        ;;
        -d|--log-dir)
        log_dir="$2"
        shift # Remove argument name from processing
        shift # Remove argument value from processing
        ;;
        -l|--library)
        library="$2"
        shift # Remove argument name from processing
        shift # Remove argument value from processing
        ;;
        -g|--genome)
        genome="$2"
        shift # Remove argument name from processing
        shift # Remove argument value from processing
        ;;
        -n|--project)
        project_id="$2"
        shift # Remove argument name from processing
        shift # Remove argument value from processing
        ;;
        -r|--read-length)
        read_length="$2"
        shift # Remove argument name from processing
        shift # Remove argument value from processing
        ;;
        -o|--options)
        options="$2"
        shift # Remove argument name from processing
        shift # Remove argument value from processing
        ;;
        -p|--pairs-file)
        pairs_file="$2"
        shift # Remove argument name from processing
        shift # Remove argument value from processing
        ;;
        -s|--samples-file)
        samples_file="$2"
        shift # Remove argument name from processing
        shift # Remove argument value from processing
        ;;
        -i|--interval-list)
        interval_list="$2"
        shift # Remove argument name from processing
        shift # Remove argument value from processing
        ;;
        -w|--wdl-file)
        workflow="$2"
        shift # Remove argument name from processing
        shift # Remove argument value from processing
        ;;
        -c|--custom-inputs)
        custom_inputs="$2"
        shift # Remove argument name from processing
        shift # Remove argument value from processing
        ;;
        -v|--skip-validate)
        skip_validate=1
        shift # Remove --initialize from processing
        ;;
        -x|--dry-run)
        dry_run=1
        shift # Remove --initialize from processing
        ;;
        -h|--help)
        print_help
        shift # Remove --initialize from processing
        ;;
    esac
done
        
        
# Required:
if [ -z "$log_dir" ]; then
    echo "Error: Missing required value for -l log dir (for run logs)" >&2
    print_help
    exit 1
fi

if [ -z "$url" ]; then
    echo "Error: Missing required value for -u cromwell server URL" >&2
    print_help
    exit 1
fi

if [ -z "$project_id" ]; then
    echo "Error: Missing required value for -p project_id" >&2
    print_help
    exit 1
fi

set -e 
set -o pipefail

script_dir=$(dirname "$0")

echo "Validate workflow..." >&2
cd ${log_dir}
womtool \
validate \
${workflow} \
--list-dependencies


# create input json
echo "Create input json and confirm files exist..." >&2
meta_command="python ${script_dir}/tools/meta.py \
    --project ${project_id} \
    --library ${library} \
    --genome ${genome} \
    --wdl-file ${workflow} \
    --options ${options}"
if [ ! -z "$read_length" ]; then
    meta_command="${meta_command} \
    --read-length ${read_length}"
fi
if [ ! -z "$custom_inputs" ]; then
    meta_command="${meta_command} \
    --custom-inputs ${custom_inputs}"
fi
if [ ! -z "$pairs_file" ]; then
    meta_command="${meta_command} \
    --pairs-file ${pairs_file}"
fi
if [ ! -z "$samples_file" ]; then
    meta_command="${meta_command} \
    --samples-file ${samples_file}"
fi
if [ ! -z "$interval_list" ]; then
    meta_command="${meta_command} \
    --interval-list ${interval_list}"
fi
if [ ! -z "$skip_validate" ]; then
    meta_command="${meta_command} \
    --skip-validate"
fi

eval ${meta_command}

# zip dependencies
echo "Zip dependencies..." >&2
cd ${script_dir}
zip dependencies.zip wdl_structs.wdl */*.wdl
cd -

echo "Precheck input json..." >&2
workflow_name=$( basename ${workflow} | sed 's/\.wdl//' )
womtool \
validate \
--inputs ${workflow_name}Input.json \
${workflow}

if [ -z "$dry_run" ]; then
    # start run:
    echo "Submit run and write log..." >&2
    cd ${log_dir}
    uuid=$( bash ${script_dir}/tools/submit.sh \
        -u ${url} \
        -w ${workflow} \
        -o ${options} \
        -d ${script_dir}/dependencies.zip \
        -p ${log_dir}/${project_id}_projectInfo.json \
        -i ${workflow_name}Input.json )
fi

echo "Done" >&2
