#!/bin/bash

# USAGE: run.sh [-h] --options OPTIONS --wdl-file WDL_FILE
#               --url URL --log-dir LOG_DIR
#               --project-name PROJECT_NAME 
#               [--library {WGS,Exome}]
#               [--genome {Human_GRCh38_full_analysis_set_plus_decoy_hla, Human_GRCh38_tcga}]
#               [--pairs-file PAIRS_FILE]
#               [--samples-file SAMPLES_FILE]
#               [--labels-file LABELS_FILE]
#               [--interval-list {SureSelect_V6plusCOSMIC.target.GRCh38_full_analysis_set_plus_decoy_hla}]
#               [--custom-inputs CUSTOM_INPUTS]
#               [--skip-validate]
#               [--dry-run]
#               [--local]
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
               --project-name PROJECT_NAME 
               [--library {WGS,Exome}]
               [--genome {Human_GRCh38_full_analysis_set_plus_decoy_hla, Human_GRCh38_tcga}]
               [--pairs-file PAIRS_FILE]
               [--samples-file SAMPLES_FILE]
               [--labels-file LABELS_FILE]
               [--interval-list {SureSelect_V6plusCOSMIC.target.GRCh38_full_analysis_set_plus_decoy_hla}]
               [--custom-inputs CUSTOM_INPUTS]
               [--skip-validate]
               [--dry-run]
               [--local]
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
  --project-name PROJECT_NAME Project name associated with account.
  --pairs-file PAIRS_FILE
                        CSV file with items that are required to have
                        \"tumor\", \"normal\" and \"pairId\" in the columns.
                        Optionally, include \"tumorBam\", \"normalBam\" columns to create 
                        \"pairInfos\" and \"normalSampleBamInfos\" automatically.
  --samples-file [SAMPLES_FILE]
                        Not generally required. If tasks only require
                        sampleId and do not use pairing information sample
                        info can be populated with a CSV file. The CSV file
                        requires a columns named [\"sampleId\"].
  --interval-list {SureSelect_V6plusCOSMIC.target.GRCh38_full_analysis_set_plus_decoy_hla}
                        File basename for interval list.If not supplied the
                        default (the SureSelect interval list for your genome)
                        will be used (only needed for future Exome workflows)
  --custom-inputs [CUSTOM_INPUTS]
                        Optional comma separated list of JSON files with custom 
                        input variables. The name of the variable in the input 
                        file must match the name of the variable in the 
                        WDL workflow. It is not required that the input 
                        specify the workflow. By default the input will 
                        be added to the top-level workflow. Any variable 
                        defined in in this JSON will overwrite any reference 
                        variable in the the config directory or workflow default.
  --labels-file         Labels json file for cromwell (not required)
  --skip-validate       Skip the step where input files are validated.
                        Otherwise all gs//: URIs will be checked to see that a
                        file exists. Disable with caution. Cromwell will launch
                        instances and run without checking. Test a small pairs
                        file to ensure all references exist and at least some
                        sample input files can be read by the current user.
  --dry-run             Skip the step where the job is submitted to cromwell-tools.
  --local               Submit to local cromwell server.
"

print_help() {
    echo "${help_top}"
    echo "DESCRIPTION: validate workflow, create input json and submit workflow to cromwell."
    echo "Script requires jq, cromwell-tools, gcloud to be in the path."
    echo "Script shows submission command and command to check status"
    echo "in the STDERR stream."
  
    echo 'Creation of input JSON:'  
    echo 'The WDL is used to determine which variables are required. '
    echo  'Required or optional variables are defined from custom inputs JSON. '
    echo  'Any required variable not defined in the custom inputs JSON will be defined from the '
    echo  'reference JSONs in the config directory (as long as the variable names are identical).'
    echo  'The pairing/sample info CSVs (--pairs-file/--samples-file) are used to create pairRelationships and '
    echo  '(if columns named tumorBam and normalBam exist) map BAMs to pairs. '
    echo  'If "production" or "external" inputs are true then validation of NYGC internal-only files is skipped '
    echo  'The pipelines are written to skip tasks that localize these files if "production" or "external" are true '
    echo  'so any inability to read these files will not negatively affect the run.'
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
        -n|--project-name)
        project_name="$2"
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
        -z|--labels-file)
        labels_file="$2"
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
        -y|--local)
        local=1
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

if [ -z "$project_name" ]; then
    echo "Error: Missing required value for -p project_name" >&2
    print_help
    exit 1
fi

log_dir=$( realpath ${log_dir})
mkdir -p ${log_dir}

set -e 
set -o pipefail


init_script_dir=$(dirname "$0")
script_dir=$( realpath ${init_script_dir})
workflow_dir=$(dirname ${workflow})
workflow_dir=$( realpath ${workflow_dir})

workflow_name=$( basename ${workflow} | sed 's/\.wdl//' )
submission_inputs=${log_dir}/${workflow_name}Input.json
project_info=${log_dir}/${project_name}_projectInfo.json

echo "Validate workflow..." >&2
womtool \
validate \
${workflow} \
--list-dependencies


# create input json
echo "Create input json and confirm files exist..." >&2
meta_command="python ${script_dir}/tools/meta.py \
    --project-name ${project_name} \
    --custom-inputs-out ${submission_inputs} \
    --file-out ${project_info} \
    --library ${library} \
    --genome ${genome} \
    --wdl-file ${workflow} \
    --options ${options}"
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

if [ ! -z "$local" ]; then
    meta_command="${meta_command} \
    --local"
fi
eval ${meta_command}

# zip dependencies
echo "Zip dependencies..." >&2
cd ${workflow_dir}
zip dependencies.zip *.wdl */*.wdl
cd -

echo "Precheck input json..." >&2

womtool \
validate \
--inputs ${submission_inputs} \
${workflow}

if [ -z "$dry_run" ]; then
    # start run:
    echo "Submit run and write log..." >&2
    
    submit_command="bash ${script_dir}/tools/submit.sh \
        -u ${url} \
        -w ${workflow} \
        -o ${options} \
        -d ${workflow_dir}/dependencies.zip \
        -p ${project_info} \
        -i ${submission_inputs}"
        
    if [ ! -z "$labels_file" ]; then
        submit_command="${submit_command} \
        -l ${labels_file}"
    fi
    uuid=$( ${submit_command} )
fi

echo "Done" >&2
