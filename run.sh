#!/bin/bash
# USAGE: uuid=$( run.sh URL WORKFLOW OPTIONS_JSON INPUTS... )
# DESCRIPTION: submit workflow to cromwell. Script requires jq to be in the path
# script returns the workflow uuid
# script shows submission command and command to check status 
# in the STDERR stream
url=$1
workflow=$2
options=$3
shift
shift
shift
inputs=$@


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
