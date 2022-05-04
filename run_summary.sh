#!/bin/bash
# USAGE: run_summary.sh -u URL -n PROJECT_NAME -g GCP_PROJECT -d log_dir [-r RUNINFO_JSON] [-i UUID] [-b BILLING EXPORT] [-p PAIRS_FILE] [-s SAMPLES_FILE]
# DESCRIPTION: submit workflow to cromwell.
# Use -l for large jobs that may timeout.
# Script requires jq, cromwell-tools, gcloud to be in the path
# Script generates summary of any workflow.
# If no *.RunInfo.json file is provided then
# the most recent file in the log dir will be used

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
#    James Roche (jroche@nygenome.org)
#
# ================== /COPYRIGHT ===============================================

init_wdl_dir=$(dirname "$0")
cd ${init_wdl_dir}
wdl_dir=$( pwd )

    
print_help() {
  echo "USAGE: run_summary.sh -u URL -n PROJECT_NAME -g GCP_PROJECT -d log_dir [-r RUNINFO_JSON] [-i UUID] [-b BILLING EXPORT] [-p PAIRS_FILE] [-s SAMPLES_FILE]" >&2
  echo "DESCRIPTION: submit workflow to cromwell." >&2
  echo "Use -l for large jobs that may timeout." >&2
  echo "Script requires jq, cromwell-tools, gcloud to be in the path." >&2
  echo "Script also requires a default_credentials_JSON to be created by the user" >&2
  echo "Run the following to generate a default credentials file " >&2
  echo "$ gcloud auth application-default login" >&2
  echo "# Script generates summary of any workflow." >&2
  exit 1
}

print_usage() {
  echo "USAGE: run_summary.sh -u URL -n PROJECT_NAME -g GCP_PROJECT -d log_dir [-r RUNINFO_JSON] [-i UUID] [-b BILLING EXPORT] [-p PAIRS_FILE] [-s SAMPLES_FILE]" >&2
  exit 1
}



for arg in "$@"; do
  case $arg in
    -b|--billing-export)
        billing_export="$2"
        shift # Remove argument name from processing
        shift # Remove argument value from processing
        ;;
    -u|--url)
        url="$2"
        shift # Remove argument name from processing
        shift # Remove argument value from processing
        ;;
    -n|--project-name)
        project_name="$2"
        shift # Remove argument name from processing
        shift # Remove argument value from processing
        ;;
    -g|--gcp-project)
        gcp_project="$2"
        shift # Remove argument name from processing
        shift # Remove argument value from processing
        ;;
    -r|--run-info-json)
        run_info_json="$2"
        shift # Remove argument name from processing
        shift # Remove argument value from processing
        ;;
    -i|--uuid)
        uuid="$2"
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
    -d|--log-dir)
        log_dir="$2"
        shift # Remove argument name from processing
        shift # Remove argument value from processing
        ;;
    -h|--help)
        print_help
        shift # Remove --initialize from processing
        ;;
  esac
done

if [ -z "$log_dir" ]; then
    echo "Error: Missing required value for -l log dir (for run logs)" >&2
    print_usage
    exit 1
fi

if [ -z "$url" ]; then
    echo "Error: Missing required value for -u cromwell server URL" >&2
    print_usage
    exit 1
fi

if [ -z "$project_name" ]; then
    echo "Error: Missing required value for -p project_name" >&2
    print_usage
    exit 1
fi

# runs from either a run_info_json or a uuid
if [ -z "$uuid" ]; then
    if [ -z "$run_info_json" ]; then
        # get most recent file if none provided
        echo "Looking up the most recent RunInfo.json file..."
        cd ${log_dir}
        run_info_json=$( find . -name "*.RunInfo.json" -print0 | xargs -r -0 ls -1 -t | head -1)
    fi
fi

start_from_uuid="True"
if [ -z "$uuid" ]; then
    if [ -z "$run_info_json" ]; then
            echo "Error: Missing value for -i uuid and no RunInfo.json file is specified or in the current directory" >&2
            print_usage
    else
        start_from_uuid="False"
    fi
    
fi

set -e 
set -o pipefail

mkdir -p ${log_dir}/tmp/
cd ${log_dir}
# Define variables
if [[ $start_from_uuid == 'True' ]]; then
    workflow_uuid=${uuid}
    command="python ${wdl_dir}/tools/create_run_info.py --uuid ${uuid} --project-name ${project_name}"
    if [ ! -z "$pairs_file" ]; then
        command="${command} \
        --pairs-file ${pairs_file}"
    fi
    if [ ! -z "$samples_file" ]; then
        command="${command} \
        --samples-file ${samples_file}"
    fi
    echo "Make relevant input json..."
    eval ${command}
    run_info_json="${project_name}.${workflow_uuid}.RunInfo.json"
else
    workflow_uuid=$( cat ${run_info_json} | jq .workflow_uuid | sed 's/"//g')
fi

output_info_file="${log_dir}/${project_name}.${workflow_uuid}_outputInfo.json"
metrics_file="${log_dir}/${project_name}.${workflow_uuid}_outputMetrics.csv"
non_retried_metrics_file="${log_dir}/${project_name}.${workflow_uuid}_outputMetrics.non_retried.csv"
plot_file="${log_dir}/${project_name}.${workflow_uuid}_outputMetrics.html"
nav="${wdl_dir}/pandoc/nav_wgs_v7"
pandoc_dir="${wdl_dir}/pandoc/"
md="${log_dir}/${project_name}.${workflow_uuid}_outputMetrics.md"
html="${log_dir}/${project_name}.${workflow_uuid}_outputMetrics.html"
header="${log_dir}/${project_name}.${workflow_uuid}_outputMetrics.header.txt"
cost_request_file="${log_dir}/tmp/cost_request.csv"

echo "Collect output files..."
collect_command="time python ${wdl_dir}/tools/collect.py \
--run-data ${run_info_json} \
--url ${url} \
--output-info-file ${output_info_file} \
--gcp-project ${gcp_project}"

eval ${collect_command}

if [[ ! -f $output_info_file ]]; then
    echo 'No output files were created by workflow so no outputInfo.json exists for further steps' >&2
    exit 0
fi

echo "Gather usage metrics..."
runtime_command="time python ${wdl_dir}/tools/runtime.py \
    --output-info-file ${output_info_file} \
    --metrics-file-prefix ${log_dir}/${project_name}. \
    --url ${url} \
    --gcp-project ${gcp_project}"

eval ${runtime_command}

#skip if monitoring_image not in options
monitored=$( cat ${output_info_file} | jq ".options.monitoring_image")
if [ -z "$monitored" ]; then
    echo "Plot usage metrics..."
    python ${wdl_dir}/tools/plot_runtime.py \
        --name ${project_name} \
        --output-info ${output_info_file} \
        --metrics ${metrics_file} \
        --non-retry-metrics ${non_retried_metrics_file} \
        --plot ${plot_file}
        
    time bash ${wdl_dir}/tools/html_printer.sh \
        ${md} \
        ${html} \
        ${header} \
        ${nav} \
        ${pandoc_dir}
else
    echo "Skipping plot usage metrics because no monitoring image was declared in options json"
fi

    
if [ ! -z "$billing_export" ]; then
    echo "Adding cost per instance id to metrics..."
    python ${wdl_dir}/tools/cost.py \
    --output-metrics "${metrics_file}" \
    --gcp-project "${gcp_project}" \
    --url "${url}"  \
    --billing-export "${billing_export}" \
    --out-file-prefix "${log_dir}/${project_name}.${workflow_uuid}"
    
    python ${wdl_dir}/tools/join.py \
    ${uuid} \
    ${costs_file} \
    ${metrics_file} \
    "${log_dir}/${project_name}.${workflow_uuid}"
else
    echo "Skipping adding cost per instance id to metrics because no --billing-export flag was used"


fi

