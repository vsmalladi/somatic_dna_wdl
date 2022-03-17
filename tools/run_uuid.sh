#!/bin/bash

url=$1
log_dir=$2
gcp_project=$3
run_info=$4
workflow_uuid=$5
cost_file=$6
wdl_dir=$7

echo "Run metadata and runtime gathering..."
time bash ${wdl_dir}/run_post.sh \
-u ${url} \
-d ${log_dir} \
-p ${gcp_project} \
-r ${run_info}

metrics_file="${log_dir}/${gcp_project}.${workflow_uuid}_outputMetrics.csv"

if [[ -s ${cost_file} ]]; then
    echo "Run cost joining..."
    python ${wdl_dir}/tools/join.py \
    ${workflow_uuid} \
    ${cost_file} \
    ${metrics_file} \
    ${log_dir}/${gcp_project}.${workflow_uuid}
fi
