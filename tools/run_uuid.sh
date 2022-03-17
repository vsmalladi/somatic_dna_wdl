#!/bin/bash

url=$1
log_dir=$2
project_name=$3
gcp_project=$4
run_info=$5
workflow_uuid=$6
cost_file=$7
wdl_dir=$8

echo "Run metadata and runtime gathering..."
time bash ${wdl_dir}/run_post.sh \
-u ${url} \
-d ${log_dir} \
-p ${project_name} \
-g ${gcp_project} \
-r ${run_info}

metrics_file="${log_dir}/${gcp_project}.${workflow_uuid}_outputMetrics.csv"

if [[ -s ${cost_file} ]]; then
    echo "Run cost joining..."
    python ${wdl_dir}/tools/join.py \
    ${workflow_uuid} \
    ${cost_file} \
    ${metrics_file} \
    ${log_dir}/${project_name}.${workflow_uuid}
fi
