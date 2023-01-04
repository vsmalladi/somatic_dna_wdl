
from cromwell_tools import api
import logging as log
import pandas as pd
import numpy as np
import datetime
import json

log.basicConfig(format='%(levelname)s:  %(message)s', level=log.INFO)

class Runtime:
    def __init__(self,
                 cromwell_connection,
                 bigquery_connection,
                 uuid,
                 billing_table_list
                 ):
        self.url = cromwell_connection.cromwell_server
        self.cromwell_auth = cromwell_connection.cromwell_auth
        self.gcp_project = cromwell_connection.gcp_project
        self.bqclient = bigquery_connection.bqclient
        self.bqstorageclient = bigquery_connection.bqstorageclient
        self.billing_tables = billing_table_list
        self.workflow_uuid = uuid

        self.loaded = self.load_metadata_api()
        self.metadata['main_workflow_id'] = self.workflow_uuid
        if not self.loaded:
            if not self.end_date:
                log.warning('No end time found for workflow: ' + uuid + '. It may still be active?')
                # drop all rows just to be safe and prevent duplication
                self.metadata.drop(self.metadata.index, inplace=True)
            if not self.run_date:
                log.warning('No start time found for workflow: ' + uuid + '. How could this happen?')
                self.metadata.drop(self.metadata.index, inplace=True)
        elif self.monitoring_image:
            self.load_runtime()

            # check for UUIDs that were missed in the metadata api call
            metadata_workflow_id = set(self.metadata.workflow_id.to_list())
            runtime_workflow_id = set(self.runtime.workflow_id.to_list())
            missing = runtime_workflow_id.difference(metadata_workflow_id)
            assert len(missing) == 0 , 'Some workflow uuids found in runtime but not in metadata ' + ' '.join(missing)

            # merge runtime info with metadata
            # self.runtime.columns = ['project_id', 'zone', 'instance_id', 'instance_name', 'preemptible',
            #                         'runtime_workflow_id', 'task_call_name', 'shard', 'attempt', 'cpu_count',
            #                         'cpu_platform', 'mem_total_gb', 'disk_mounts', 'disk_total_gb',
            #                         'vm_start_time']

            # start time from metadata doesn't really capture tasks that wait for a vm, runtime start_time more accurate
            self.runtime.rename(columns={'runtime_workflow_id': 'workflow_id',
                                         'start_time': 'vm_start_time'},
                                inplace=True)

            self.runtime = self.reset_instance_id(self.runtime)
            # rename the columns before merge, then replace where runtime NaN
            self.metadata.rename(columns={'mem_total_gb':'mem_total_gb_runtime',
                                          'disk_mounts': 'disk_mounts_runtime',
                                          'disk_total_gb': 'disk_total_gb_runtime',
                                          }, inplace=True)

            self.metadata = pd.merge(self.metadata, self.runtime[['instance_name', 'instance_id',
                                                                  'cpu_platform', 'mem_total_gb', 'disk_mounts',
                                                                  'disk_total_gb', 'vm_start_time']]
                                     .drop_duplicates(subset=['instance_id']),
                                     on='instance_name', how='left')
            # self.metadata['vm_summary'] = self.metadata.cpu_platform + ' ' + self.metadata.disk_type

            if not self.metadata['vm_start_time'].isnull().all():
                self.metadata['wait_time_m'] = (self.metadata['vm_start_time'] - self.metadata['start_time']) / pd.Timedelta(minutes=1)
                self.metadata['vm_runtime_m'] = (self.metadata['end_time'] - self.metadata['vm_start_time']) / pd.Timedelta(minutes=1)
                self.metadata['cpu_time_m'] = self.metadata['vm_runtime_m'] * self.metadata['cpu_count']
            else:
                self.metadata['wait_time_m'] = np.nan
                self.metadata['vm_runtime_m'] = np.nan
                self.metadata['cpu_time_m'] = np.nan
            self.metadata['mem_total_gb'].fillna(self.metadata['mem_total_gb_runtime'], inplace=True)
            self.metadata['disk_mounts'].fillna(self.metadata['disk_mounts_runtime'], inplace=True)
            self.metadata['disk_total_gb'].fillna(self.metadata['disk_total_gb_runtime'], inplace=True)
            self.metadata.drop(columns=['mem_total_gb_runtime', 'disk_mounts_runtime',
                                        'disk_total_gb_runtime'], inplace=True)
            self.metadata = self.reset_instance_id(self.metadata)
            self.load_metrics()
            if len(self.metrics) > 0:
                self.metrics = self.reset_instance_id(self.metrics)
                log.info('add custom cols')
                self.add_custom_cols()
            # remove non-run rows
            # log.info('remove non-run rows')
            # self.metadata_cache = self.metadata.loc[pd.isnull(self.metadata[['instance_name']]).any(axis=1)].copy()
            # self.metadata = self.metadata.dropna(subset=['instance_name']).copy()


            if self.billing_tables != ["Unknown"]:
                self.add_cost()
                self.join_cost()
        else:
            pass


    def load_metadata_api(self):
        '''
        project_id    zone    instance_name    preemptible    workflow_name
            workflow_id    task_call_name    shard    attempt    start_time
            end_time    execution_status    cpu_count    mem_total_gb
            disk_mounts    disk_total_gb    disk_types    docker_image    inputs

        convert metadata from api into table for pandas dataframe
            Includes start_time, end_time, execution_status, cpu_count,
            mem_total_gb, task_call_name, workflow_name

            inputs
        Returns False if no instance name is associated with an end time
        '''
        log.info('Loading Metadata...')
        metadata = self.load_from_api_no_expand(self.workflow_uuid)
        # whole-run level metadata
        main_workflow_name = metadata.get('workflowName', '')
        labels = metadata.get('submittedFiles', {}).get('labels', '{}') # dictionary in a string
        # get rid of this? just get it from the root workflow (or can labels be updated throughout workflow?)

        # also get some other stuff
        options = json.loads(metadata.get('submittedFiles', {}).get('options', '{}'))
        self.monitoring_image = options.get('monitoring_image', None)
        self.outputs_dir = options.get('final_workflow_outputs_dir', None)
        self.workflowLog = metadata.get('workflowLog', '')
        self.workflowRoot = metadata.get('workflowRoot', '')

        # instance-level results
        task_call_names = []
        wdl_task_names = []
        sub_workflow_names = []
        workflow_names = []
        execution_statuss = []
        mem_total_gbs = []
        cpu_counts = []
        disk_types = []
        disk_total_gbs = []
        docker_images = []
        disk_mounts = []
        attempts = []
        backend_statuss = []
        preemptibles = []
        instance_names = []
        zones = []
        project_ids = []
        return_codes = []
        workflow_ids = []
        shards = []
        start_times = []
        end_times = []
        localization_ms = []
        delocalization_ms = []
        user_action_ms = []
        vm_start_times = []
        inputs = []
        self.tasks = []
        self.call_lister(metadata['calls'],
                         main_workflow_name,
                         workflow_id=self.workflow_uuid
                         )
        # order   call_type, attempt, workflowName
        for workflow_id, task_call_name, call, workflow_name in self.tasks:
            task_call_names.append(task_call_name)
            wdl_task_names.append(call.get('backendLabels', {}).get('wdl-task-name', task_call_name))
            sub_workflow_names.append(call.get('backendLabels', {}).get('cromwell-sub-workflow-name', task_call_name.split('.')[0]))
            workflow_names.append(workflow_name.split('.')[-1])
            execution_statuss.append(call['executionStatus'])
            # replace with dict.get with default values
            if 'runtimeAttributes' in call:
                mem_total_gbs.append(float(call['runtimeAttributes'].get('memory', '0 GB').split()[0]))
                cpu_counts.append(int(call['runtimeAttributes'].get('cpu', 0)))
                disk_types.append(call['runtimeAttributes'].get('disks', '0 NA').split()[-1])
                disk_total_gbs.append(call['runtimeAttributes'].get('disks', '0 NA').split()[-2])
                docker_images.append(call['runtimeAttributes'].get('docker', ''))
            else:
                mem_total_gbs.append(np.nan)
                cpu_counts.append(np.nan)
                disk_types.append('')
                disk_total_gbs.append(np.nan)
                docker_images.append('')
            attempts.append(call['attempt'])
            disk_mounts.append('Unknown')
            if 'backendStatus' in call:
                backend_status = call['backendStatus'] # Preempted or Failed or Success
            # elif 'callCaching' in call and call['callCaching']['hit']:
            elif call.get('callCaching', {}).get('hit', False):
                backend_status = 'CacheHit'
            else:
                backend_status = 'Error'
            backend_statuss.append(backend_status)
            preemptibles.append(call.get('preemptible', np.nan))
            instance_names.append(call.get('jes', {}).get('instanceName', np.nan))
            zones.append(call.get('jes', {}).get('zone', np.nan))
            project_ids.append(call.get('jes', {}).get('googleProject', np.nan))
            # At least one case of workflow with {"backendStatus": "Success"} and no return code
            # command below is example subworkflow oldest parent f1e8d232-da38-4d88-9b9d-9cc7b6e1e4be
            # sh wdl_port/tools/get_metadata.sh -u "https://cromwell-compbio01.nygenome.org" -w "11fe69cc-b98c-4f1a-80d0-0c13f640cbaf" | jq '.calls | ."BicSeq2.uniqReadsTumor" | .[1]'
            return_codes.append(call.get('returnCode', np.nan))
            # get workflow uuid
            if 'id' in call:
                workflow_ids.append(call['id'])
            elif 'subWorkflowId' in call:
                workflow_ids.append(call['subWorkflowId'])
            elif workflow_id:
                workflow_ids.append(workflow_id)
            else:
                workflow_ids.append(np.nan)
            shards.append(call['shardIndex'])
            start_times.append(call.get('start', np.nan))
            end_times.append(call.get('end', np.nan))
            # get execution events
            execution_events = call.get("executionEvents", [])
            localization_ms.append(self.execution_time(execution_events, "Localization"))
            delocalization_ms.append(self.execution_time(execution_events, "Delocalization"))
            user_action_ms.append(self.execution_time(execution_events, "UserAction"))
            # vm_start_times.append(self.vm_start_time(execution_events))
            inputs.append(call.get('inputs', {}))
        # output
        self.metadata = pd.DataFrame({'task_call_name' : task_call_names,
                                      'wdl_task_name': wdl_task_names,
                                      'sub_workflow_name': sub_workflow_names,
                                      'workflow_name' : workflow_names,
                                      'execution_status' : execution_statuss,
                                      'mem_total_gb' : mem_total_gbs,
                                      'cpu_count' : cpu_counts,
                                      'disk_types' : disk_types,
                                      'disk_total_gb' : disk_total_gbs,
                                      'docker_image' : docker_images,
                                      'disk_mounts' : disk_mounts,
                                      'attempt'  : attempts,
                                      'backend_status' : backend_statuss,
                                      'preemptible' : preemptibles,
                                      'instance_name' : instance_names,
                                      'zone' : zones,
                                      'project_id' : project_ids,
                                      'return_code' : return_codes,
                                      'workflow_id' : workflow_ids,
                                      'shard' : shards,
                                      'start_time' : start_times,
                                      'end_time' : end_times,
                                      'localization_m' : localization_ms,
                                      'delocalization_m' : delocalization_ms,
                                      'user_action_time_m' : user_action_ms,
                                      # 'vm_start_time' : vm_start_times,
                                      'inputs' : inputs
                                      })
        self.metadata['main_workflow_name'] = main_workflow_name
        if len(self.metadata) > 0: # below failed with empty df
            self.metadata['disk_type'] = self.metadata.apply(lambda row: row.disk_types.replace('[\'', '').replace('\']', '') if str(row.disk_types) != 'nan' else '', axis=1)
        self.metadata['labels'] = labels
        self.metadata['inputs'] = self.metadata['inputs'].apply(json.dumps)
        self.deduplicate_metadata()
        self.metadata['start_time'] = pd.to_datetime(self.metadata['start_time'],
                                                     format="%Y-%m-%dT%H:%M:%S.%fZ",
                                                     utc=True)
        # self.metadata['vm_start_time'] = pd.to_datetime(self.metadata['vm_start_time'],
        #                                                 format="%Y-%m-%dT%H:%M:%S.%fZ",
        #                                                 utc=True)
        self.metadata['end_time'] = pd.to_datetime(self.metadata['end_time'],
                                                   format="%Y-%m-%dT%H:%M:%S.%fZ",
                                                   utc=True)

        # self.run_date = self.load_date(self.metadata.start_time.min())
        # self.end_date = self.load_date(self.metadata.end_time.max())
        # if not self.end_date:
        #     return False
        # why use call min/max instead of the more inclusive cromwell times. I guess to minimize the conversions...
        run_date = pd.to_datetime(metadata.get('start', ''),
                                  format="%Y-%m-%dT%H:%M:%S.%fZ",
                                  utc=True)
        end_date = pd.to_datetime(metadata.get('end', ''),
                                  format="%Y-%m-%dT%H:%M:%S.%fZ",
                                  utc=True)
        self.run_date = self.load_date(run_date)
        self.end_date = self.load_date(end_date)
        if not (self.run_date and self.end_date):
            return False

        # this is important for big query to work properly.
        # Example: uuid 16f5c5b3-cfcb-4a32-8c0d-6f3c72bd2842 returns no costs
        # Just padding the end should be fine but doing both for safety
        self.padded_start = (run_date - datetime.timedelta(days=1)).strftime(format='%Y-%m-%d')
        self.padded_end = (end_date + datetime.timedelta(days=7)).strftime(format='%Y-%m-%d')

        # if there are no calls we need some kind of placeholder data to fill tasks table to prevent re-querying
        if len(self.metadata) == 0:
            filldata = {'main_workflow_id': self.workflow_uuid,
                        'workflow_id': self.workflow_uuid,
                        'main_workflow_name': main_workflow_name,
                        'workflow_name': main_workflow_name,
                        'start_time': run_date,
                        'end_time': end_date,
                        'labels': labels,
                        'inputs': dict(), # needs to be json-like string to be able to load into table. Null or NaN fail
                        }
            tdf = pd.DataFrame(filldata, index=[0])
            self.metadata = pd.concat([self.metadata, tdf], ignore_index=True)
        return True

    def load_from_api_no_expand(self, uuid):
        # replacing with cromwell api
        response = api.metadata(uuid=uuid,
                                auth=self.cromwell_auth,
                                expandSubWorkflows=False)
        if response.status_code == 200:
            log.info('Retrieved metadata from api for workflow: ' + uuid + ' (Root workflow: ' + self.workflow_uuid + ')')
            data = response.json()
            for k, v in data['calls'].items():
                for i, t in enumerate(data['calls'][k]):
                    if 'subWorkflowId' in t:
                        if k == 'MergeVcf.MergeCallers': # this is problematic and fails need to make response smaller
                            subdata = self.load_from_api_with_prune(t['subWorkflowId'])
                        else:
                            subdata = self.load_from_api_no_expand(t['subWorkflowId'])
                        if subdata:
                            data['calls'][k][i]['subWorkflowMetadata'] = subdata
                            del data['calls'][k][i]['subWorkflowId']
                        else:
                            data['calls'][k][i] = t
                    else:
                        data['calls'][k][i] = t
            return data
        elif response.status_code == 503:  # Typically production, No timely response (full 5 minute wait). Returns string-like bytes
            log.error('Failed to get metadata from api for workflow: ' + uuid + ' (Root workflow: ' + self.workflow_uuid + ')')
            log.error('Response code: ' + str(response.status_code))
            log.error('Response: ' + response.content.decode('utf-8'))
            return {}
        elif response.status_code == 500:  # Typically compbio/dlp, rows served exceeds limit (quick failure). Returns json-like bytes
            log.error('Failed to get metadata from api for workflow: ' + uuid + ' (Root workflow: ' + self.workflow_uuid + ')')
            log.error('Response code: ' + str(response.status_code))
            log.error('Response: ' + response.content.decode('utf-8'))
            return {}
        else:  # placeholder for other errors, can assess if they arise (ie. 40X for uuid doesn't exist)
            log.error('Failed to get metadata from api for workflow: ' + uuid + ' (Root workflow: ' + self.workflow_uuid + ')')
            log.warning('Response code: ' + str(response.status_code))
            log.warning('Response: ' + response.content.decode('utf-8'))
            return {}

    def load_from_api_with_prune(self, uuid):
        '''
        This isn't called recursively like load_from_api_no_expand
        MergeVcf.MergeCallers doesn't have any subworkflows so should be okay but not ideal
        this is structured poorly, is there a better way. Don't really want to use while loop
        Also, each failure takes 5 minutes before server response of 503, any way to speed up?
        '''
        response = api.metadata(uuid=uuid,
                                auth=self.cromwell_auth,
                                expandSubWorkflows=False)
        if response.status_code != 200: # too large, try again without unnecessary keys
            exclude = ['stdout', 'stderr', 'commandLine', 'backendLogs', 'outputs']
            log.warning('MergeCallers metadata query failed. Retrying excluding keys ' + ', '.join(exclude))
            response = api.metadata(uuid=uuid,
                                    auth=self.cromwell_auth,
                                    expandSubWorkflows=False,
                                    excludeKey=exclude)
            if response.status_code != 200:  # too large, try again without executionEvents
                exclude.append('executionEvents')
                log.warning('MergeCallers metadata query failed. Retrying excluding keys ' + ', '.join(exclude))
                response = api.metadata(uuid=uuid,
                                        auth=self.cromwell_auth,
                                        expandSubWorkflows=False,
                                        excludeKey=exclude)
                if response.status_code != 200:  # too large, try again without 'callCaching'
                    exclude.append('callCaching')
                    log.warning('MergeCallers metadata query failed. Retrying excluding keys ' + ', '.join(exclude))
                    response = api.metadata(uuid=uuid,
                                            auth=self.cromwell_auth,
                                            expandSubWorkflows=False,
                                            excludeKey=exclude)
                    if response.status_code != 200:
                        log.error('Failed to get metadata from api for workflow: ' + uuid)
                        return {}
                    else:
                        log.info('Retrieved partial metadata from api for workflow: ' + uuid)
                        return response.json()
                else:
                    log.info('Retrieved partial metadata from api for workflow: ' + uuid)
                    return response.json()
            else:
                log.info('Retrieved partial metadata from api for workflow: ' + uuid)
                return response.json()
        else:
            log.info('Retrieved metadata from api for workflow: ' + uuid)
            return response.json()

    def execution_time(self, events, description):
        '''takes list of execution events and returns the time in minutes for any description
        note that there can be duplicated descriptions (ie. waiting) this returns first in list but not necessarily chronologically orderd
        '''
        try:
            for event in events:
                if event['description'] == description:
                    start = event['startTime']
                    end = event['endTime']
                    return self.sec_between(start, end) / 60.0
        except KeyError:
            pass
        return np.nan

    def call_lister(self, call_dictionary, workflowName, workflow_id=False):
        '''takes as input the metadata "calls" section
        and yields all task calls nested in metadata
        '''
        # workflow_id was being updated and then used in the else. Leaf will be assigned previous call's workflow id
        # just use a different variable name for the subworkflows
        for call_type in call_dictionary:
            for i, attempt in enumerate(call_dictionary[call_type]): # subnode is a list within the dictionary
                if 'subWorkflowMetadata' in attempt:
                    if 'id' in attempt['subWorkflowMetadata']:
                        sub_workflow_id = attempt['subWorkflowMetadata']['id']
                    # subWorkflowId should never appear since no longer using expandSubWorkflows, leaving for now though
                    elif 'subWorkflowId' in attempt['subWorkflowMetadata']:
                        sub_workflow_id = attempt['subWorkflowMetadata']['subWorkflowId']
                    self.call_lister(attempt['subWorkflowMetadata']['calls'],
                                     attempt['subWorkflowMetadata']['workflowName'],
                                     workflow_id=sub_workflow_id)
                else:
                    self.tasks.append([workflow_id, call_type, attempt, workflowName])

    def sec_between(self, d1,d2):
        d1 = datetime.datetime.strptime(d1, "%Y-%m-%dT%H:%M:%S.%fZ")
        d2 = datetime.datetime.strptime(d2, "%Y-%m-%dT%H:%M:%S.%fZ")
        return abs((d2 - d1).seconds)

    def load_date(self, date_record):
        try:
            return date_record.strftime(format='%Y-%m-%d')
        except ValueError:
            return False

    def vm_start_time(self, events):
        '''Grab the vm start time from cromwell API'''
        # this doesn't work at all with repetitive, unordered lists. Worker released vs. Assigned
        try:
            for event in events:
                if event['description'].startswith('Worker '):
                    # get the time a vm worker was assigned to the task
                    vm_start_time = event['startTime']
                    return vm_start_time
        except KeyError:
            pass
        return np.nan

    def deduplicate_metadata(self):
        '''deduplicate and add in sample ids'''
        self.metadata['sample_id'] = self.metadata['inputs'].apply(self.assign_samples)
        hashable = ['instance_name', 'workflow_id', 'attempt', 'shard']
        self.metadata = self.metadata.drop_duplicates(subset=hashable)

    def assign_samples(self, input):
        '''  whole new approach since there is no RunInfo.json
        do a half-hearted attempt where known input keys exist, can add to this later after a few tries
        '''
        x = json.loads(input)
        keys = x.keys()
        if 'sampleId' in keys:
            return x['sampleId']
        elif 'pairName' in keys:
            return x['pairName']
        elif 'name' in keys:
            return x['name']
        # above are obvious and self-explanatory, below are task-specific and likely unwieldy long-term
        elif 'fastqs' in keys: # Skewer
            return x['fastqs'].get('sampleId', 'Unknown')
        elif 'fastqsAlign' in keys: # AlignBwaMem
            return x['fastqsAlign'].get('sampleId', 'Unknown')
        elif 'contaminationFile' in keys: # SomaticQcCheck
            return x['contaminationFile'].split('/')[-1].split('.')[0]
        elif 'wgsMetricsFile' in keys: # BamQcCheck
            return x['wgsMetricsFile'].split('/')[-1].split('.')[0]
        elif 'vcf' in keys: # genotypedFilteredCompressVcf, haplotypecallerCompressVcf
            return x['vcf'].split('/')[-1].split('.')[0]
        elif 'sortedVcfPath' in keys: # genotypedFilteredMergeSortVcf
            return x['sortedVcfPath'].split('.')[0]
        elif 'vcfCompressed' in keys: # IndexVcf
            return x['vcfCompressed'].split('/')[-1].split('.')[0]
        elif 'outprefix' in keys: # Angsd, FastNgsAdmix
            return x['outprefix']
        elif 'insertSizeMetrics' in keys: # IndexVcf
            return x['insertSizeMetrics'].split('/')[-1].split('.')[0]
        elif 'outputTablePath' in keys: # IndexVcf
            return x['outputTablePath'].split('.')[0]
        elif 'orderedVcfPath' in keys: # IndexVcf
            return x['orderedVcfPath'].split('.')[0]
        elif 'laneBam' in keys: # fixmate, shortalignmark
            return x['laneBam'].split('/')[-1].split('.')[0].split('_')[0]
        elif 'bamBase' in keys and x['bamBase']: # Gridss
            return x['bamBase'].split('.')[0].split('_')[0]
        # and an extra messy way to handle callers
        elif 'inVcf' in keys:
            inVcf = x['inVcf']
            if '/call-Strelka2/' in inVcf:
                return inVcf.split('/call-Strelka2/')[-1].split('.')[0]
            elif '/call-SvabaWgs/' in inVcf:
                return inVcf.split('/call-SvabaWgs/')[-1].split('.')[0]
            elif '/call-MantaWgs/' in inVcf:
                return inVcf.split('/call-MantaWgs/')[-1].split('.')[0]
            elif '/call-Gatk4MergeSortVcf/' in inVcf:
                return inVcf.split('/call-Gatk4MergeSortVcf/')[-1].split('.')[0]
            elif '/call-unfilteredGatk4MergeSortVcf/' in inVcf:
                return inVcf.split('/call-unfilteredGatk4MergeSortVcf/')[-1].split('.')[0]
            elif '/call-filteredGatk4MergeSortVcf/' in inVcf:
                return inVcf.split('/call-filteredGatk4MergeSortVcf/')[-1].split('.')[0]
            else:
                return 'Unknown'
        else:
            return 'Unknown'

    def load_metrics(self, chunk_size=15):
        '''convert big query metrics table into pandas dataframe
        Includes instance_id, cpu_used_percent, mem_used_gb, disk_used_gb,
        disk_read_iops, disk_write_iops.
        Use to find max and plot resource usage over time.
        '''
        log.info('Loading Metrics...')
        self.instance_ids = self.runtime.dropna(subset=['instance_id']).instance_id.unique().tolist()

        metrics_dfs = []
        chunk_list = lambda all_instance_ids, chunk: [all_instance_ids[i:i+chunk]
                                                      for i in range(0, len(all_instance_ids), chunk)]
        chunked_list = chunk_list(self.instance_ids, chunk_size)
        for instance_ids in chunked_list:
            instance_id_list = self.format_query_list(instance_ids, integer=True)
            query_string = '''
            SELECT * FROM `''' + self.gcp_project + '''.cromwell_monitoring.metrics`
                WHERE DATE(timestamp) >= "''' + self.padded_start + '''"
                AND DATE(timestamp) <= "''' + self.padded_end + '''"
                AND instance_id IN ''' + instance_id_list
            metrics = (
                self.bqclient.query(query_string)
                    .result()
                    .to_dataframe()
            )
            metrics_dfs.append(metrics)
        if metrics_dfs:
            self.metrics = pd.concat(metrics_dfs)
        else:
            self.metrics = pd.DataFrame()
        log.info('Done loading Metrics')

    def load_runtime(self):
        '''convert big query runtime table into pandas dataframe
        Includes both instance_id and workflow_id.
        Use to map workflow_id to instance_ids

        Add complete list

        '''
        log.info('Loading Runtime...')
        workflow_uuids = self.metadata.dropna(subset=['workflow_id']).workflow_id.unique().tolist()
        self.sub_workflow_uuid_list = self.format_query_list(workflow_uuids)
        instance_names = self.metadata.dropna(subset=['instance_name']).instance_name.unique().tolist()
        self.instance_name_list = self.format_query_list(instance_names)
        query_string = '''
        SELECT * FROM `''' + self.gcp_project + '''.cromwell_monitoring.runtime`
            WHERE DATE(start_time) >= "''' + self.padded_start + '''"
            AND DATE(start_time) <= "''' + self.padded_end + '''"
            AND (workflow_id IN ''' + self.sub_workflow_uuid_list + '''
                OR instance_name IN ''' + self.instance_name_list  + ')'
        self.runtime = (
            self.bqclient.query(query_string)
                .result()
                .to_dataframe()
        )
        # de-listify some columns
        self.runtime['disk_mounts'] = self.runtime['disk_mounts'].apply(lambda x: x[0])
        self.runtime['disk_total_gb'] = self.runtime['disk_total_gb'].apply(lambda x: x[0])
        log.info('Done loading Runtime')

    def add_custom_cols(self):
        '''Add per instance id values from metrics to metadata'''
        self.metrics['row_max_cpu_percent'] = self.metrics['cpu_used_percent'].apply(max)
        self.metrics['max_cpu_used_percent'] = self.metrics.groupby(['instance_id']).row_max_cpu_percent.transform(max)
        self.metrics['row_max_disk_used_gb'] = self.metrics['disk_used_gb'].apply(max)
        self.metrics['max_disk_used_gb'] = self.metrics.groupby(['instance_id']).row_max_disk_used_gb.transform(max)
        self.metadata = pd.merge(self.metadata, self.metrics[['instance_id',
                                                              'max_cpu_used_percent',
                                                              'max_disk_used_gb']].drop_duplicates(subset=['instance_id']),
                                 on='instance_id', how='left')
        self.metrics['max_mem_used_gb'] = self.metrics.groupby(['instance_id']).mem_used_gb.transform(max)
        self.metadata = pd.merge(self.metadata, self.metrics[['instance_id',
                                                              'max_mem_used_gb']].drop_duplicates(subset=['instance_id']),
                                 on='instance_id', how='left')

    def format_query_list(self, sub_workflow_uuids, integer=False):
        if integer:
            sub_workflow_uuid_list = '(' + ', '.join([str(id) for id in sub_workflow_uuids]) + ')'
        else:
            sub_workflow_uuid_list = '("' + '", "'.join([str(id) for id in sub_workflow_uuids]) + '")'
        return sub_workflow_uuid_list

    def reset_instance_id(self, data):
        data['instance_id'] = data['instance_id'].fillna(-1)
        data['instance_id'] = data['instance_id'].astype(int)
        data['instance_id'] = data['instance_id'].astype(str)
        data['instance_id'] = data['instance_id'].replace('-1', np.nan)
        return data

    def get_flow_runtime(self, metadata, ids, col='sample_subworkflow_run_time_h'):
        ''' For sub workflow and workflows this will also include wait time and
        return the end-to-end wallclock time for a particular id.
        For tasks this will only return the longest runtime.
        '''
        metadata['grouped_end_time'] = metadata.groupby(ids).end_time.transform(max)
        metadata['grouped_start_time'] = metadata.groupby(ids).start_time.transform(min)
        metadata[col] = (metadata['grouped_end_time'] - metadata['grouped_start_time']).dt.total_seconds() / 3600.0
        metadata.drop(columns=['grouped_end_time', 'grouped_start_time'], inplace=True)
        return metadata

    def load_plot_metrics(self, metadata):
        '''add cpu_time, max_mem, wall_clock, cpu_used_percent, mem_efficiency, mem_used_gb
        also add cpu_used_percent", "disk_total_gb",  "disk_used_gb",  "cpu_count",   "cpu_platform", "preemptible"

        '''
        log.info('Calculating runtime metrics...')
        metadata['run_time'] = metadata['end_time'] - metadata['start_time']
        metadata['run_time_m'] = metadata['run_time'].dt.total_seconds() / 60.0
        # get cpu_time
        metadata['cpu_time_m'] = metadata['run_time_m'] * metadata['cpu_count']
        # ============================
        #     Tasks w/in subworkflow
        # ============================
        log.info('Tasks w/in subworkflow...')
        # sub task wallclock time from start to finish
        metadata['sample_task_run_time_h'] = metadata['run_time_m'] / 60.0
        metadata['sample_task_core_h'] = metadata['cpu_time_m'] /60.0

        # ============================
        #       sub workflows
        # ============================
        log.info('Subworkflows...')
        grouped = metadata.groupby(['sample_id', 'workflow_name'])
        sub_metadata = grouped.agg(sample_subworkflow_core_h=pd.NamedAgg(column='sample_task_core_h',
                                                                         aggfunc=sum)).reset_index()
        metadata = pd.merge(metadata, sub_metadata, on=['sample_id', 'workflow_name'], how='left')

        # sub workflows wallclock time from start to finish
        metadata = self.get_flow_runtime(metadata, ids=['sample_id', 'workflow_name'],
                                         col='sample_subworkflow_run_time_h').copy()
        if 'max_mem_used_gb' in metadata.columns:
            metadata['subworkflow_max_mem_g'] = metadata.groupby(['sample_id',
                                                                  'workflow_name']).max_mem_used_gb.transform(max)

        # ============================
        #       workflows
        # ============================
        # full workflow wallclock time from start to finish for an id
        log.info('Workflow...')
        grouped = metadata.groupby('sample_id')
        sub_metadata = grouped.agg(sample_workflow_core_h=pd.NamedAgg(column='sample_task_core_h',
                                                                      aggfunc=sum)).reset_index()
        metadata = pd.merge(metadata, sub_metadata, on=['sample_id'], how='left')
        metadata = self.get_flow_runtime(metadata, ids=['sample_id'], col='sample_workflow_run_time_h').copy()
        if 'max_mem_used_gb' in metadata.columns:
            metadata['workflow_max_mem_g'] = metadata.groupby(['sample_id']).max_mem_used_gb.transform(max)
        log.info('Done calculating runtime metrics...')

        return metadata

    def add_cost(self):
        # need to loop when >1 billing table
        cost_dfs = []
        log.info('Loading Cost: ' + self.workflow_uuid + '...')
        for b in self.billing_tables:
            query_string = '''SELECT
        (SELECT value from UNNEST(labels) where key = 'cromwell-workflow-id' limit 1) as cromwell_workflow_id,
        (SELECT value from UNNEST(labels) where key = 'cromwell-sub-workflow-name' limit 1) as sub_workflow_name,
        (SELECT value from UNNEST(labels) where key = 'wdl-task-name' limit 1) as wdl_task_name,
        sum(cost) as cost,
        sku.description as service_type
        FROM `''' + b + '''` 
        WHERE DATE(_PARTITIONTIME) >= "''' + self.padded_start + '''"
            AND DATE(_PARTITIONTIME) <= "''' + self.padded_end + '''"
            AND project.id = "''' + self.gcp_project + '''"
            AND EXISTS (SELECT 1 FROM UNNEST(labels) WHERE key = 'cromwell-workflow-id' and value in ("cromwell-''' + self.workflow_uuid + '''"))
        group by sub_workflow_name, wdl_task_name, cromwell_workflow_id, sku.description
            '''
            df = (
                self.bqclient.query(query_string)
                    .result()
                    .to_dataframe()
            )
        cost_dfs.append(df)
        self.costs = pd.concat(cost_dfs)
        log.info('Done loading Cost')

    def join_cost(self):

        self.metadata['cromwell_workflow_id'] = 'cromwell-' + self.metadata['main_workflow_id']
        unique_roots = self.metadata.drop_duplicates(['main_workflow_name', 'cromwell_workflow_id']).set_index(
            'cromwell_workflow_id')
        self.costs['sub_workflow_name'] = np.where(self.costs['sub_workflow_name'].isna(),
                                                   self.costs['cromwell_workflow_id'].map(unique_roots['main_workflow_name']),
                                                   self.costs['sub_workflow_name']
                                                   )
        self.costs[['cost_type', 'machine_type', 'preemptible']] = self.costs.apply(mapper, axis=1, result_type="expand")
        grouper = ['cromwell_workflow_id', 'sub_workflow_name', 'wdl_task_name']
        wide_costs = pd.pivot_table(self.costs,
                                    values='cost',
                                    index=grouper,
                                    columns='cost_type',
                                    aggfunc='sum',
                                    fill_value=0
                                    ).reset_index()
        self.metadata = pd.merge(self.metadata, wide_costs, on=grouper, how='outer')
        cost_opts = ['capacity_cost', 'core_cost', 'ram_cost', 'egress_cost', 'unknown_cost']
        cost_cols = [c for c in cost_opts if c in self.metadata.columns]
        # the costs are still for the full group, divide them evenly
        self.metadata['group_count'] = self.metadata.groupby(grouper)['task_call_name'].transform('count')
        for c in cost_cols:
            self.metadata[c] = self.metadata[c] / self.metadata['group_count']
        self.metadata['total_cost'] = self.metadata[cost_cols].sum(axis=1)
        self.metadata['group_cost'] = self.metadata.groupby(grouper)['total_cost'].transform('sum')
        self.metadata['group_runtime'] = self.metadata.groupby(grouper)['vm_runtime_m'].transform('sum')
        self.metadata['runtime_scaled_total_cost'] = self.metadata['group_cost'] * \
                                                     (self.metadata['vm_runtime_m'] / self.metadata['group_runtime'])
        # add the service type, treat as set. Almost always going to be a single type
        core_costs = self.costs[self.costs['cost_type'] == 'core_cost'].copy()
        service_types = core_costs.groupby(grouper).agg({'machine_type': set}).reset_index()
        service_types['machine_type'] = service_types.machine_type.str.join(';')
        self.metadata = pd.merge(self.metadata, service_types, on=grouper, how='outer')
        missing = self.costs[self.costs['cost_type'] == 'unknown_cost']['service_type'].unique().tolist()
        if len(missing) > 0:
            log.warning('Following service_types are not yet described in mapper(). Add them to join.py: ')
            log.warning(', '.join(missing))
        self.metadata.drop(columns=['group_count', 'group_cost', 'group_runtime', 'cromwell_workflow_id'], inplace=True)
        cost_cols_map = {col : 'avg_' + col for col in cost_cols + ['total_cost']}
        self.metadata = self.metadata.rename(columns=cost_cols_map)
        self.metadata['vm_summary'] = self.metadata.cpu_platform + ' ' + self.metadata.machine_type + ' ' + self.metadata.disk_types


def is_preemptible(service):
    if 'preemptible' in service:
        return True
    else:
        return False


def formatter(service, drop=[], truncate=False, suffix=False):
    '''takes a string service type, split to list of words, drop some useless words,
    drop up to truncate word and suffix. can't just keep it due to 'storage'-> 'capacity' '''
    word_list = service.split()
    clean_words = [a for a in word_list if a not in drop]
    if truncate:
        try:
            truncate_index = clean_words.index(truncate)
        except ValueError:  # storage services are called capacity
            truncate_index = -1
        words = clean_words[:truncate_index]
    else:
        words = clean_words
    if suffix:
        words.append(suffix)
    return '_'.join(words)


def mapper(row):
    '''cost_type & preemptible are pretty clearly defined. machine_type not so much'''
    service_type = row['service_type'].lower()
    if 'storage' in service_type or 'capacity' in service_type:
        cost_type = 'capacity_cost'
        preemptible = np.nan
        machine_type = formatter(service_type,
                                 drop=['backed', 'storage', 'to', 'preemptible', 'vms'],
                                 truncate='capacity',
                                 suffix='capacity'
                                 )
    elif 'egress' in service_type:
        cost_type = 'egress_cost'
        preemptible = np.nan
        machine_type = formatter(service_type,
                                 truncate='egress',
                                 suffix='egress'
                                 )
    elif 'core' in service_type:
        cost_type = 'core_cost'
        preemptible = is_preemptible(service_type)
        machine_type = formatter(service_type,
                                 drop=['preemptible'],
                                 truncate='instance'
                                 )
    elif 'ram' in service_type:
        cost_type = 'ram_cost'
        preemptible = is_preemptible(service_type)
        machine_type = formatter(service_type,
                                 drop=['preemptible'],
                                 truncate='instance'
                                 )
    else:
        cost_type = 'unknown_cost'
        machine_type = 'unknown'
        preemptible = np.nan
    return cost_type, machine_type, preemptible
