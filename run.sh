#!/bin/bash
# USAGE: uuid=$( run.sh -u CROMWELL_URL -w WORKFLOW -o OPTIONS_JSON -i INPUTS... )
# DESCRIPTION: submit workflow to cromwell. Script requires jq to be in the path
# script returns the workflow uuid
# script shows submission command and command to check status 
# in the STDERR stream


print_usage() {
  printf "Usage: uuid=$( run.sh -u URL -w WORKFLOW -o OPTIONS_JSON -i INPUTS... )"
}

while getopts 'u:w:o:d:p:i:' flag; do
  case "${flag}" in
    u) url="${OPTARG}" ;;
    w) workflow="${OPTARG}" ;;
    o) options="${OPTARG}" ;;
    d) dependencies="${OPTARG}" ;;
    p) project_data="${OPTARG}" ;;
    i) inputs="${OPTARG}" ;;
    *) print_usage
       exit 1 ;;
  esac
done

script_dir=$(dirname "$0")


module load gcloud cromwell-tools jq
# compose submission...
cmd="cromwell-tools submit --url ${url} \
--username $(gcloud secrets versions access latest --secret="cromwell_username") \
--password $(gcloud secrets versions access latest --secret="cromwell_password") \
-w ${workflow} -i ${inputs} -o ${options} -d ${dependencies}"

# run workflow...
echo $cmd >&2
response=$( eval $cmd )

# report status command and return UUID
uuid=$( echo ${response}  | jq -r ".id")

echo "cromwell-tools status --url ${url} \
--username $(gcloud secrets versions access latest --secret="cromwell_username") \
--password $(gcloud secrets versions access latest --secret="cromwell_password") \
--uuid ${uuid}" >&2

echo ${uuid}

python ${script_dir}/tools/log.py \
--project-data ${project_data} \
--uuid ${uuid} \
--inputs ${inputs}
