#!/bin/bash
# USAGE: run_summary.sh -u URL -n PROJECT_NAME -g GCP_PROJECT -d log_dir [-r RUNINFO_JSON] [-i UUID] [-b BILLING EXPORT] [--from-billing] [-p PAIRS_FILE] [-s SAMPLES_FILE]
# DESCRIPTION: monitor or summarize cromwell workflow.
# Script requires jq to be in the path
# Script generates summary of any workflow.
# If no *.RunInfo.json OR uuid file is provided then
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
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

init_wdl_dir=$(dirname "$0")
cd ${init_wdl_dir}
wdl_dir=$( pwd )


help_long="-h, --help          Show this help message and exit
  --billing-export BILLING     The name of the table for the SQL query
                               of the billing table.
  --from-billing               Skip regenerating runtime metrics. Use existing
                               metrics files; add the cost values from the billing
                               table; and plot the results (including billing).
  --gcp-project GCP_PROJECT    GCP project id for api queries.
  --run-info-json RUNINFO_JSON Optional file that includes any sampleIds and the
                               main workflow UUID.
  --url URL                    Cromwell server URL (required)
  --log-dir LOG_DIR            Output directory for all logs and reports
                               related to this workflow UUID (required)
  --project-name PROJECT_NAME  Project name associated with account.
  --pairs-file PAIRS_FILE
                               Optional, CSV file with items that are required to have
                               \"tumor\", \"normal\" and \"pairId\" in the columns.
                               Optionally, include \"tumorBam\", \"normalBam\" columns to create
                               \"pairInfos\" and \"normalSampleBamInfos\" automatically.
  --samples-file [SAMPLES_FILE]
                               Not generally required. If tasks only require
                               sampleId and do not use pairing information sample
                               info can be populated with a CSV file. The CSV file
                               requires a columns named [\"sampleId\"].

"

print_help() {
  echo "USAGE: run_summary.sh -u URL -n PROJECT_NAME -g GCP_PROJECT -d log_dir [-r RUNINFO_JSON] [-i UUID] [-b BILLING EXPORT] [--from-billing] [-p PAIRS_FILE] [-s SAMPLES_FILE]" >&2
  echo "DESCRIPTION: monitor or summarize cromwell workflow." >&2
  echo "Script requires jq to be in the path." >&2
  echo "If no *.RunInfo.json file OR uuid is provided then" >&2
  echo "the most recent file in the log dir will be used" >&2
  echo "Script generates summary of any workflow." >&2
  echo "${help_long}" >&2
  exit 1
}

print_usage() {
  echo "USAGE: run_summary.sh -u URL -n PROJECT_NAME -g GCP_PROJECT -d log_dir [-r RUNINFO_JSON] [-i UUID] [-b BILLING EXPORT] [--from-billing] [-p PAIRS_FILE] [-s SAMPLES_FILE]" >&2
  exit 1
}



for arg in "$@"; do
  case $arg in
    -b|--billing-export)
        billing_export="$2"
        shift # Remove argument name from processing
        shift # Remove argument value from processing
        ;;
    -z|--from-billing)
        start_from_billing=1
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
    echo "Error: Missing required value for -d log dir (for run logs)" >&2
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
costs_file="${log_dir}/${project_name}.${workflow_uuid}_outputCosts.csv"

if [ -z "$start_from_billing" ]; then
    echo "Collect output files..."
    collect_command="time python ${wdl_dir}/tools/collect.py \
    --run-data ${run_info_json} \
    --url ${url} \
    --output-info-file ${output_info_file} \
    --gcp-project ${gcp_project}"

    eval ${collect_command}
fi

if [[ ! -f $output_info_file ]]; then
    echo 'No output files were created by workflow so no outputInfo.json exists for further steps' >&2
    exit 0
fi

if [ -z "$start_from_billing" ]; then
    echo "Gather usage metrics..."
    runtime_command="time python ${wdl_dir}/tools/runtime.py \
        --output-info-file ${output_info_file} \
        --metrics-file-prefix ${log_dir}/${project_name}. \
        --url ${url} \
        --gcp-project ${gcp_project}"

    eval ${runtime_command}
fi

# billing (optional)
plot_metrics_file=${metrics_file}
if [ ! -z "$billing_export" ]; then
    echo "Adding cost per instance id to metrics..."
    python ${wdl_dir}/tools/cost.py \
    --output-metrics "${metrics_file}" \
    --workflow-uuid "${workflow_uuid}" \
    --gcp-project "${gcp_project}" \
    --url "${url}"  \
    --billing-export "${billing_export}" \
    --out-file-prefix "${log_dir}/${project_name}.${workflow_uuid}"
    
    count=$(cat ${costs_file} |  wc -l)
    
    if [ "$count" -gt 1 ]; then
        python ${wdl_dir}/tools/join.py \
        ${workflow_uuid} \
        ${costs_file} \
        ${metrics_file} \
        "${log_dir}/${project_name}.${workflow_uuid}"
        # make plots with table that includes cost
        plot_metrics_file="${log_dir}/${project_name}.${workflow_uuid}.outputMetrics.cost.csv"
    else
        echo -e "WARNING:${YELLOW}The billing-export query failed to find any matching records. Check to see if you are pointing to the correct table and using the correct gcp project.${NC}"
    fi
else
    echo -e "WARNING:${YELLOW}Skipping adding cost per instance id to metrics because no --billing-export flag was used.${NC}"
fi

# plotting (optional)
#skip if monitoring_image not in options
monitored=$( cat ${output_info_file} | jq ".options.monitoring_image")
if [ ! "$monitored" = "null" ]; then
    echo "Plot usage metrics..."
    python ${wdl_dir}/tools/plot_runtime.py \
        --name ${project_name} \
        --output-info ${output_info_file} \
        --metrics ${plot_metrics_file} \
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


