#!/bin/bash
# USAGE: run_post.sh -u URL -p PROJECT_ID -d LOG_DIR [-r RUNINFO_JSON]
# DESCRIPTION: submit workflow to cromwell.
# Script requires jq, cromwell-tools, gcloud to be in the path
# Script generates summary of any workflow.
# If no *.RunInfo.json file is provided then
# the most recent file in the log dir will be used

init_wdl_dir=$(dirname "$0")
cd ${init_wdl_dir}
wdl_dir=$( pwd )

print_help() {
  echo "USAGE: run_post.sh -u URL -p PROJECT_ID -d LOG_DIR [-r RUNINFO_JSON]"
  echo "DESCRIPTION: submit workflow to cromwell."
  echo "Script requires jq, cromwell-tools, gcloud to be in the path."
  echo "Script also requires GOOGLE_APPLICATION_CREDENTIALS to be defined"
  echo "# Script generates summary of any workflow."
  exit 1
}

print_usage() {
  echo "USAGE: run_post.sh -u URL -p PROJECT_ID -d log_dir [-r RUNINFO_JSON] " >&2
  exit 1
}

while getopts 'u:p:d:r:h' flag; do
  case "${flag}" in
    u) url="${OPTARG}" ;;
    p) project_id="${OPTARG}" ;;
    r) run_info_json="${OPTARG}" ;;
    d) log_dir="${OPTARG}" ;;
    h) print_help ;;
    \?) print_usage; echo "Unknown option: $OPTARG" >&2 ;;
    *) print_usage; echo "Unimplemented option: $OPTARG" >&2 ;;
  esac
done


if [ -z "$GOOGLE_APPLICATION_CREDENTIALS" ]; then
    echo "Error: Define GOOGLE_APPLICATION_CREDENTIALS in order to run script" >&2
    print_usage
    exit 1
fi

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

if [ -z "$project_id" ]; then
    echo "Error: Missing required value for -p project_id" >&2
    print_usage
    exit 1
fi


if [ -z "$run_info_json" ]; then
    # get most recent file if none provided
    cd ${log_dir}
    run_info_json=$( find . -name "*.RunInfo.json" -print0 | xargs -r -0 ls -1 -t | head -1)
fi

set -e 
set -o pipefail

cd ${log_dir}
workflow_uuid=$( cat ${run_info_json} | jq .workflow_uuid | sed 's/"//g' )

output_info_file="${log_dir}/${project_id}.${workflow_uuid}_outputInfo.json"
metrics_file="${log_dir}/${project_id}.${workflow_uuid}_outputMetrics.csv"
non_retried_metrics_file="${log_dir}/${project_id}.${workflow_uuid}_outputMetrics.non_retried.csv"
plot_file="${log_dir}/${project_id}.${workflow_uuid}_outputMetrics.html"
nav="${wdl_dir}/pandoc/nav_wgs_v7"
pandoc_dir="${wdl_dir}/pandoc/"
md="${log_dir}/${project_id}.${workflow_uuid}_outputMetrics.md"
html="${log_dir}/${project_id}.${workflow_uuid}_outputMetrics.html"
header="${log_dir}/${project_id}.${workflow_uuid}_outputMetrics.header.txt"


echo "Collect subworkflow uuids..."
time python ${wdl_dir}/tools/collect.py \
--run-data ${run_info_json} \
--url ${url}

if [[ ! -s $output_info_file ]]; then
    echo 'No output files were created by workflow so no outputInfo.json exists for further steps' >&2
    exit 0
fi

echo "Gather usage metrics..."
time python ${wdl_dir}/tools/runtime.py \
    ${output_info_file}

echo "Plot usage metrics..." 
python ${wdl_dir}/tools/plot_runtime.py \
    --project-id ${project_id} \
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

