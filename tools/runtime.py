import sys
from google.cloud import bigquery
from google.cloud import bigquery_storage
from google_auth_oauthlib import flow
import google.auth
import json
import datetime
import logging as log
import pandas as pd
import argparse
import os 
import numpy as np
import subprocess
import pprint
import ast


try:
    import Colorer
except ImportError:
    pass

pd.set_option('display.max_columns', 500)
log.basicConfig(format='%(levelname)s:  %(message)s', level=log.INFO)


class Runtime():
    def __init__(self,
                 output_info,
                 url,
                 metrics_file_prefix,
                 large_run=False,
                 limit=1000000,
                 gcp_project=False,
                 test=False):
        self.large_run = large_run
        self.test= test
        self.parent_dir = os.path.abspath(os.path.dirname(__file__))
        self.output_info = output_info
        self.workflow_uuid = output_info['workflow_uuid']
        self.url = url
        # find existing project id info
        self.prep_ids()
        self.metrics_file = metrics_file_prefix + self.workflow_uuid + '_outputMetrics.csv'
        self.non_retried_metrics_file = metrics_file_prefix + self.workflow_uuid + '_outputMetrics.non_retried.csv'
        if os.path.isfile(self.metrics_file):
            self.loaded = True
        else:
            # find relevant project info
            self.limit = str(limit)
            self.metrics_limit = str(limit * 4)
            # login
            self.credentials, self.gcp_project = self.login()
            if gcp_project:
                self.default_project = False
                self.gcp_query_project = gcp_project
            else:
                self.default_project = True
                self.gcp_query_project = self.gcp_project
            self.get_sample_info()
            # gather big cloud tables
            self.bqclient = bigquery.Client(project=self.gcp_project , 
                                            credentials=self.credentials)
            self.bqstorageclient = bigquery_storage.BigQueryReadClient(credentials=self.credentials)
            self.loaded = self.load_metadata_api()
            if not self.loaded:
                log.warning('No endtime found for workflow: ' + output_info['workflow_uuid'])
            else:
                self.load_runtime()
                # merge runtime info with metadata
                self.runtime.columns = ['project_id', 'zone', 'instance_id', 'instance_name', 'preemptible',
                               'runtime_workflow_id', 'task_call_name', 'shard', 'attempt', 'cpu_count',
                               'cpu_platform', 'mem_total_gb', 'disk_mounts', 'disk_total_gb',
                               'start_time']
                self.runtime['instance_id']
                self.runtime['instance_id'] = self.runtime['instance_id'].fillna(-1)
                self.runtime['instance_id'] = self.runtime['instance_id'].astype(int)
                self.runtime['instance_id'] = self.runtime['instance_id'].astype(str)
                self.runtime['instance_id'] = self.runtime['instance_id'].replace('-1', np.nan)
                cols_to_replace = ['mem_total_gb', 'disk_mounts',  'disk_total_gb']
                cols_to_keep = [col for col in self.metadata.columns if not col in cols_to_replace]
                self.metadata = pd.merge(self.metadata[cols_to_keep], self.runtime[['instance_name', 'instance_id',
                                                              'cpu_platform', 'mem_total_gb', 'disk_mounts', 
                                                              'disk_total_gb']].drop_duplicates(subset=['instance_id']), 
                                         on='instance_name', how='left')
                print(self.metadata.shape)
                self.metadata['instance_id'] = self.metadata['instance_id'].fillna(-1)
                self.metadata['instance_id'] = self.metadata['instance_id'].astype(int)
                self.metadata['instance_id'] = self.metadata['instance_id'].astype(str)
                self.metadata['instance_id'] = self.metadata['instance_id'].replace('-1', np.nan)
                self.load_metrics()
                print(self.metadata.shape)
                self.metrics['instance_id'] = self.metrics['instance_id'].fillna(-1)
                self.metrics['instance_id'] = self.metrics['instance_id'].astype(int)
                self.metrics['instance_id'] = self.metrics['instance_id'].astype(str)
                self.metrics['instance_id'] = self.metrics['instance_id'].replace('-1', np.nan)
                print('add custom cols')
                self.add_custom_cols()
                print(self.metadata.shape)
                # remove non-run rows
                print('remove non-run rows')
                self.metadata_cache = self.metadata.loc[pd.isnull(self.metadata[['instance_name']]).any(axis=1)].copy()
                self.metadata = self.metadata.dropna(subset=['instance_name']).copy()
                self.non_retry_metadata = self.metadata[self.metadata.backend_status == 'Success'].copy()
                print(self.metadata.shape)
                print('write metrics')
                self.metadata.to_csv(self.metrics_file + '.intermediate.csv', index=False, float_format='{:f}'.format, encoding='utf-8')
                # if no tasks are finished
                if self.metadata[self.metadata.execution_status != 'RUNNING'].empty:
                    self.metadata = pd.concat([self.metadata, self.metadata_cache], ignore_index=True)
                    self.metadata.to_csv(self.metrics_file, index=False)
                else:
                    print('make plot metrics')
                    self.metadata_full = self.metadata.copy()
                    self.metadata = self.load_plot_metrics(self.metadata)
                    print('make non retry plot metrics')
                    self.non_retry_metadata = self.load_plot_metrics(self.non_retry_metadata)
                    self.metadata = pd.concat([self.metadata, self.metadata_cache], ignore_index=True)
                    print('write metrics')
                    self.metadata.to_csv(self.metrics_file, index=False, float_format='{:f}'.format, encoding='utf-8')
                    self.non_retry_metadata.to_csv(self.non_retried_metrics_file, index=False, float_format='{:f}'.format, encoding='utf-8')
                    
    def login(self):
        '''Run the following to generate a default credentials file
    
        $ gcloud auth application-default login
    
        https://google-auth.readthedocs.io/en/latest/reference/google.auth.html#google.auth.default.'''
        credentials, gcp_project = google.auth.default()
        return credentials, gcp_project
                    
    def add_custom_cols(self):
        '''Add per instance id values from metrics to metadata'''
        self.metadata['max_cpu_used_percent'] = self.metadata.apply(lambda row: self.get_cpu_used(row.instance_id), axis=1)
        self.metadata['max_disk_used_gb'] = self.metadata.apply(lambda row: self.get_disk_used_gb(row.instance_id), axis=1)
        self.metrics['max_mem_used_gb'] = self.metrics.groupby(['instance_id']).mem_used_gb.transform(max)
        self.metadata = pd.merge(self.metadata, self.metrics[['instance_id', 'max_mem_used_gb']].drop_duplicates(subset=['instance_id']), on='instance_id', how='left')
        
    @staticmethod
    def read(file):
        with open(file) as output_info_file:
            output_info = json.load(output_info_file)
            return output_info
        
    def get_task_core_h(self, grouped, row, ids=['id', 
                                                 'task_call_name', 
                                                 'workflow_name']):
        ''' get runtime in core hours for the task (within a workflow) '''
        return (row.cpu_count * float(row.run_time_m)) / 60
            
    def get_flow_runtime(self, metadata, ids, col='sample_subworkflow_run_time_h'):
        ''' For sub workflow and workflows this will also include wait time and
        return the end-to-end wallclock time for a particular id.
        For tasks this will only return the longest runtime.
        '''
        metadata['grouped_end_time'] = metadata.groupby(ids).end_time.transform(max)
        metadata['grouped_start_time'] = metadata.groupby(ids).start_time.transform(min)
        metadata[col] = metadata.apply(lambda row: (row.grouped_end_time - row.grouped_start_time).total_seconds() / 60.0 / 60.0, axis=1)
        cols = [col for col in metadata.columns if not col in ['grouped_end_time', 'grouped_start_time']]
        return metadata[cols]
        
    def load_plot_metrics(self, metadata):
        '''add cpu_time, max_mem, wall_clock, cpu_used_percent, mem_efficiency, mem_used_gb
        also add cpu_used_percent", "disk_total_gb",  "disk_used_gb",  "cpu_count",   "cpu_platform", "preemptible"
        
        '''
        log.info('Calculating runtime metrics...')
        metadata['run_time'] = metadata['end_time'] - metadata['start_time']
        metadata['run_time_m'] = metadata.apply(lambda row:
                                                     (row.run_time.total_seconds() / 60.0), axis=1)
        # get cpu_time
        metadata['cpu_time_m'] = metadata.apply(lambda row:
                                                (row.run_time_m
                                                 * row.cpu_count),
                                                axis=1)
        # ==============
        #     Instances
        # ==============
#         log.info('Instances...')
#         grouped = metadata.groupby('task_call_name')
#         metadata['mean_task_core_h'] = metadata.apply(lambda row:
#                                                            (row.cpu_count
#                                                             * float(grouped.mean().run_time_m[row.task_call_name])) / 60,
#                                                            axis=1)
#         # task wallclock time from start to finish
#         metadata['mean_task_run_time_h'] = metadata.apply(lambda row:
#                                                           float(grouped.mean().run_time_m[row.task_call_name]) / 60,
#                                                           axis=1)
        
        # ============================
        #     Tasks w/in subworkflow
        # ============================
        log.info('Tasks w/in subworkflow...')
        grouped = metadata.groupby(['id', 'task_call_name', 'workflow_name'])
        # sub task wallclock time from start to finish
        metadata['sample_task_run_time_h'] = metadata.apply(lambda row: row.run_time_m / 60.0, axis=1)
        metadata['sample_task_core_h'] = metadata.apply(lambda row: row.cpu_time_m / 60.0, axis=1)
        # ============================
        #       sub workflows
        # ============================
        log.info('Subworkflows...')
        grouped = metadata.groupby(['id', 'workflow_name'])
        sub_metadata = grouped.agg(sample_subworkflow_core_h=pd.NamedAgg(column='sample_task_core_h', aggfunc=sum)).reset_index()
        metadata['sample_subworkflow_core_h'] = metadata.apply(lambda row:
                                                               sub_metadata[(sub_metadata.id == row.id)
                                                                            & (sub_metadata.workflow_name == row.workflow_name)].sample_subworkflow_core_h.tolist()[0], 
                                                               axis=1)
        # sub workflows wallclock time from start to finish
        metadata = self.get_flow_runtime(metadata, ids=['id', 'workflow_name'], col='sample_subworkflow_run_time_h')
        metadata['subworkflow_max_mem_g'] = metadata.groupby(['id', 'workflow_name']).max_mem_used_gb.transform(max)

        # ============================
        #       workflows
        # ============================
        # full workflow wallclock time from start to finish for an id
        log.info('Workflow...')
        grouped = metadata.groupby('id')
        sub_metadata = grouped.agg(sample_workflow_core_h=pd.NamedAgg(column='sample_task_core_h', aggfunc=sum)).reset_index()
        metadata['sample_workflow_core_h'] = metadata.apply(lambda row: 
                                                            sub_metadata[(sub_metadata.id == row.id)].sample_workflow_core_h.tolist()[0],
                                                            axis=1)
        metadata = self.get_flow_runtime(metadata, ids=['id'], col='sample_workflow_run_time_h')
        metadata['workflow_max_mem_g'] = metadata.groupby(['id']).max_mem_used_gb.transform(max)
        log.info('Done calculating runtime metrics...')
        return metadata
    
    def load_end_date(self, date_record):
        try:
            return date_record.strftime(format='%Y-%m-%d')
        except ValueError:
            return False     

    def modify_date(self, run_date):
        return '-'.join(run_date.split('-')[0:3])
    
    def get_sample_info(self):
        if 'sampleIds' in self.output_info['project_data']:
            self.sample_count = len(self.output_info['project_data']['sampleIds'])
        else:
            self.sample_count = 0
        if 'listOfPairRelationships' in self.output_info['project_data']:
            self.pair_count = len([pair['pairId'] for pair in self.output_info['project_data']['listOfPairRelationships']])
        else:
            self.pair_count = 0

    def format_query_list(self, sub_workflow_uuids, integer=False):
        if integer:
            sub_workflow_uuid_list = '(' + ', '.join([str(id) for id in sub_workflow_uuids]) + ')'
        else:
            sub_workflow_uuid_list = '("' + '", "'.join([str(id) for id in sub_workflow_uuids]) + '")'
        return sub_workflow_uuid_list
    
    def add_ids(self):
        '''deduplicate and add in sample ids '''
        self.metadata['id'] = self.metadata.apply(lambda row: self.match_id(row.inputs).replace('\n', ' '), axis=1)
#         unhashable = ['disk_mounts', 'disk_total_gb', 'disk_used_gb', 'disk_types', 'inputs']
#         hashable = ['instance_name', 'workflow_id', 'attempt', 'shard']
#         self.metadata = self.metadata.drop_duplicates(subset=hashable)
        
    def get_localization(self, event):
        '''Handle changes in output organization'''
        start = False
        end = False
        if 'description' in event:
            if event['description'] == 'Localization' and \
                    'startTime' in event and \
                    'endTime' in event:
                start = event['startTime']
                end = event['endTime']
        return start, end
    
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
        '''
        log.info('Loading Metadata...')
        if self.large_run:
            all_metadata = [self.load_from_api()]
            if len(all_metadata) == 0:
                all_metadata = self.load_from_large_api()
        else:
            all_metadata = self.load_from_large_api()
        # whole-run level metadata
        main_workflow_name = all_metadata[0]['workflowName']
        labels = all_metadata[0]['submittedFiles']['labels'] # dictionary in a string        
        # instance-level results
        task_call_names = []
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
        inputs = []
        for metadata in all_metadata:
            self.tasks = []
            self.call_lister(metadata['calls'], 
                          metadata['workflowName'])
            # order   call_type, attempt, workflowName
            for workflow_id, task_call_name, call, workflow_name in self.tasks:
                with open('example_call.json', 'w') as project_info_file:
                    json.dump(call, project_info_file, indent=4)
                task_call_names.append(task_call_name)
                workflow_names.append(workflow_name.split('.')[-1])
                execution_statuss.append(call['executionStatus'])
                if 'runtimeAttributes' in call:
                    mem_total_gbs.append(int(call['runtimeAttributes']['memory'].split()[0]))
                    cpu_counts.append(int(call['runtimeAttributes']['cpu']))
                    disk_types.append(call['runtimeAttributes']['disks'].split()[-1])
                    disk_total_gbs.append(call['runtimeAttributes']['disks'].split()[-2])
                    docker_images.append(call['runtimeAttributes']['docker'])
                else:
                    mem_total_gbs.append(0)
                    cpu_counts.append(0)
                    disk_types.append('')
                    disk_total_gbs.append(0)
                    docker_images.append('')
                attempts.append(call['attempt'])
                disk_mounts.append('Unknown')
                
                if 'backendStatus' in call:
                    with open('example_backend.json', 'w') as project_info_file:
                        json.dump(call, project_info_file, indent=4)
                    backend_status = call['backendStatus'] # Preempted or Failed or Success
                elif 'callCaching' in call and call['callCaching']['hit']:
                    backend_status = 'CacheHit'
                else:
                    backend_status = 'CacheHit'
#                     backend_status = 'Error'
                backend_statuss.append(backend_status)
                if backend_status in ['CacheHit']:
                    preemptibles.append(np.nan)
                    instance_names.append(np.nan)
                    zones.append(np.nan)
                    project_ids.append(np.nan)
                elif backend_status in ['Running']:
                    preemptibles.append(call['preemptible'])
                    instance_names.append(np.nan)
                    zones.append(np.nan)
                    project_ids.append(call['jes']['googleProject'])
                else:
                    preemptibles.append(call['preemptible'])
                    instance_names.append(call['jes']['instanceName'])
                    zones.append(call['jes']['zone'])
                    project_ids.append(call['jes']['googleProject'])
                if backend_status in ['Preempted', 'Running', 'CacheHit']:
                    return_codes.append(np.nan)
                elif backend_status in ['CacheHit'] and \
                        not 'returnCode' in call:
                    # handle PAPI errors
                    return_codes.append(np.nan)
                elif not 'returnCode' in call:
                    print(backend_status)
                    # handle PAPI errors
                    return_codes.append(np.nan)
                else:
                    return_codes.append(call['returnCode'])
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
                start_times.append(call['start'])
                if 'end' in call:
                    # get localization time
                    localization_m = np.nan
                    for event in call["executionEvents"]:
                        start_possible, end_possible = self.get_localization(event)
                        if start_possible:
                            start = start_possible
                            end = end_possible
                            localization_m = self.sec_between(start, end) / 60.0
                    localization_ms.append(localization_m)
                    end_times.append(call['end'])
                else:
                    localization_ms.append(np.nan)
                    end_times.append(np.nan)
                inputs.append(call['inputs'])
        # output
        self.metadata = pd.DataFrame({'task_call_name' : task_call_names,
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
                                'inputs' : inputs
        })
        self.metadata = self.prune()
        self.metadata['main_workflow_name'] = main_workflow_name
        self.metadata['labels'] = labels
        self.add_ids()
        self.metadata['start_time'] = pd.to_datetime(self.metadata['start_time'],
                                                     format="%Y-%m-%dT%H:%M:%S.%fZ")
        self.metadata['end_time'] = pd.to_datetime(self.metadata['end_time'],
                                                     format="%Y-%m-%dT%H:%M:%S.%fZ")
        self.run_date = self.load_end_date(self.metadata.start_time.min())
        self.end_date = self.load_end_date(self.metadata.end_time.max())
        if not self.end_date:
            return False
        return True
    
    def prune(self):
        '''skip results from subsubworkflow uuids that may return less information if a better result is available'''
        dfs = []
        for instance_name in self.metadata.instance_name.unique().tolist():
            match = self.metadata[(self.metadata.instance_name == instance_name) &
                                  ([str(m) != 'nan' for m in self.metadata.workflow_id])].copy()
            null_count = np.Inf
            if not match.empty:
                for index, row in match.iterrows():
                    current_null_count = len([True for col in row if str(col) == 'nan'])
                    if current_null_count < null_count:
                        min_row = row
                dfs.append(pd.DataFrame([min_row]))
        return pd.concat(dfs, ignore_index=True)

    def get_cromwell_time(self, timestamp):
        return datetime.datetime.strptime(timestamp, "%Y-%m-%dT%H:%M:%S.%fZ")
    
    def sec_between(self, d1, d2):
        d1 = self.get_cromwell_time(d1)
        d2 = self.get_cromwell_time(d2)
        return abs((d2 - d1).seconds)

    def call_lister(self, call_dictionary, workflowName, statement='main', workflow_id=False):
        '''takes as input the metadata "calls" section 
        and yields all task calls nested in metadata
        '''
        for call_type in call_dictionary:
            for i, attempt in enumerate(call_dictionary[call_type]): # subnode is a list within the dictionary
                if 'subWorkflowMetadata' in attempt:
                    if 'id' in attempt['subWorkflowMetadata']:
                        workflow_id = attempt['subWorkflowMetadata']['id']
                    elif 'subWorkflowId' in attempt['subWorkflowMetadata']:
                        workflow_id = attempt['subWorkflowMetadata']['subWorkflowId']
                    self.call_lister(attempt['subWorkflowMetadata']['calls'], 
                                     attempt['subWorkflowMetadata']['workflowName'],
                                     statement='sub',
                                     workflow_id= workflow_id)
                else:
                    self.tasks.append([workflow_id, call_type, attempt, workflowName])

    def load_from_large_api(self):
        '''
            return json from api
        '''
        metadata_list = []
        # get main workflow results
        metadata = self.load_uuid_api(uuid=self.workflow_uuid, script='/get_main_metadata.sh')
        if metadata['status'] != 'Succeeded':
            log.error('Call to api failed')
            print(metadata)
            sys.exit(1)
        metadata_list.append(metadata)
        # get subworkflow results
        for uuid in self.output_info['sub_workflow_uuids']:
            if not uuid == self.workflow_uuid:
                metadata = self.load_uuid_api(uuid=uuid, script='/get_metadata.sh')
                if metadata['status'] != 'Succeeded':
                    log.error('Call to api failed')
                    print(metadata)
                    sys.exit(1)
                metadata_list.append(metadata)
        return metadata_list
    
    def load_uuid_api(self, uuid, script='/get_metadata.sh'):
        '''
            return json from api
        '''
        script = self.parent_dir + script
        if self.default_project:
            args = ['bash', script,
                    '-u', self.url,
                    '-w', uuid]
        else:
            args = ['bash', script,
                    '-u', self.url,
                    '-w', uuid,
                    '-g', self.gcp_query_project]
        try:
            result = subprocess.run(args,
                                    check=True,
                                    stdout=subprocess.PIPE).stdout.decode('utf-8')
        except subprocess.CalledProcessError as err:
            log.warning(err.output.decode('utf-8'))
            log.error('Failed to get metadata from api for workflow: ' + uuid)
            return {}
#         with open('api_example.json', 'w') as input_info_file:
#             json.dump(json.loads(result), input_info_file, indent=4)
        metadata = json.loads(result)
        return metadata

    def load_from_api(self, script='/get_metadata.sh'):
        '''
            return json from api
        '''
        script = self.parent_dir + script
        metadata = self.load_uuid_api(uuid=self.workflow_uuid, script=script)
        if metadata['status'] != 'Succeeded':
            self.large_run = True
            metadata = []
        return metadata
    
    def format_instance_id(self):
        instance_id_list = '(' + ', '.join([str(int(i)) for i in self.instance_ids]) + ')'
        return instance_id_list
        
    def load_runtime(self):
        '''convert big query runtime table into pandas dataframe
        Includes both instance_id and workflow_id.
        Use to map workflow_id to instance_ids
        '''
        log.info('Loading Runtime...')
        workflow_uuids = self.metadata.dropna(subset=['workflow_id']).workflow_id.unique().tolist()
        self.sub_workflow_uuid_list = self.format_query_list(workflow_uuids)
        instance_names = self.metadata.dropna(subset=['instance_name']).instance_name.unique().tolist()
        self.instance_name_list = self.format_query_list(instance_names)
        query_string = '''
        SELECT * FROM `''' + self.gcp_query_project + '''.cromwell_monitoring.runtime`
            WHERE DATE(start_time) >= "''' + self.run_date + '''"
            AND DATE(start_time) <= "''' + self.end_date + '''"
            AND workflow_id IN ''' + self.sub_workflow_uuid_list + '''
            AND instance_name IN ''' + self.instance_name_list + '''
            LIMIT ''' + self.limit + '''
        '''
        self.runtime = (
            self.bqclient.query(query_string)
            .result()
            .to_dataframe(bqstorage_client=self.bqstorageclient)
        )
        log.info('Done loading Runtime')
    
    def load_metrics(self, chunk_size=15):
        '''convert big query metrics table into pandas dataframe
        Includes instance_id, cpu_used_percent, mem_used_gb, disk_used_gb, 
        disk_read_iops, disk_write_iops.
        Use to find max and plot resource usage over time.
        '''
        log.info('Loading Metrics...')
        if self.test:
            self.metrics = pd.read_csv('metrics.csv')
        else:
            self.instance_ids = self.runtime.dropna(subset=['instance_id']).instance_id.unique().tolist()
    
            metrics_dfs = []
            chunk_list = lambda all_instance_ids, chunk: [all_instance_ids[i:i+chunk] 
                                                          for i in range(0, len(all_instance_ids), chunk)]
            chunked_list = chunk_list(self.instance_ids, chunk_size)
            for instance_ids in chunked_list:
                instance_id_list = self.format_query_list(instance_ids, integer=True)
                query_string = '''
                SELECT * FROM `''' + self.gcp_query_project + '''.cromwell_monitoring.metrics`
                    WHERE DATE(timestamp) >= "''' + self.run_date + '''"
                    AND DATE(timestamp) <= "''' + self.end_date + '''"
                    AND instance_id IN ''' + instance_id_list
                metrics = (
                    self.bqclient.query(query_string)
                    .result()
                    .to_dataframe(bqstorage_client=self.bqstorageclient)
                )
                metrics_dfs.append(metrics)
            self.metrics = pd.concat(metrics_dfs)
        log.info('Done loading Metrics')
    
    def get_cpu_used(self, instance_id):
        metrics_instance = self.metrics[self.metrics.instance_id == instance_id].copy()
        cpu_used_percents = []
        for cpu_used_percent_list in metrics_instance.cpu_used_percent.tolist():
            if isinstance(cpu_used_percent_list, str):
                cpu_used_percent_list = ast.literal_eval(cpu_used_percent_list)
            for cpu_used_percent in cpu_used_percent_list:
                cpu_used_percents.append(cpu_used_percent)
        if len(cpu_used_percents) > 0:
            return max(cpu_used_percents)
        return np.nan
    
    def get_disk_used_gb(self, instance_id):
        metrics_instance = self.metrics[self.metrics.instance_id == instance_id].copy()
        disk_used_gbs = []
        for disk_used_gb_list in metrics_instance.disk_used_gb.tolist():
            if isinstance(disk_used_gb_list, str):
                disk_used_gb_list = ast.literal_eval(disk_used_gb_list)
            for disk_used_gb in disk_used_gb_list:
                disk_used_gbs.append(disk_used_gb)
        if len(disk_used_gbs) > 0:
            return max(disk_used_gbs)
        return np.nan
    
    def get_disk_total_gb(self, instance_id):
        runtime_instance = self.runtime[self.runtime.instance_id == instance_id].copy()
        disk_total_gbs = []
        for disk_total_gb_list in runtime_instance.disk_total_gb.tolist():
            for disk_total_gb in disk_total_gb_list:
                disk_total_gbs.append(disk_total_gb)
        if len(disk_total_gbs) > 0:
            return max(disk_total_gbs)
        return np.nan
                
    def id_matches(self, identifier, id):
        '''test if it is a likely substring of the WDL input string'''
        if (id + '.' in identifier \
                or id == identifier \
                or id + '_' in identifier):
            return True
        return False
    
    def search_sample_ids(self, possible_indentifiers):
        '''sample_id '''
        matched_sample_ids = []
        if self.sample_ids:
            for sample_id in self.sample_ids:
                pair_in_name = False
                group_in_name = False
                match = False
                #look for longer samples id that the current sample id may be a substring of
                other_file_sample_ids = [other_sample_id for other_sample_id in self.sample_ids if
                                         (not other_sample_id == sample_id)]
                other_file_sample_ids = [other_sample_id for other_sample_id in other_file_sample_ids
                                         if sample_id in other_sample_id]
                self.other_file_sample_ids = list(other_file_sample_ids)
                self.other_file_sample_ids += [other_sample_id + '.' for other_sample_id in other_file_sample_ids]
                self.other_file_sample_ids += [other_sample_id + '_' for other_sample_id in other_file_sample_ids]
                self.other_file_sample_ids += [other_sample_id + '/' for other_sample_id in other_file_sample_ids]
                for identifier in possible_indentifiers:
                    # if id is in input string
                    if self.id_matches(identifier, sample_id):
                        # but no other id is also in the input string
                        other_possible_match = any([self.id_matches(identifier, id) for id in self.other_file_sample_ids])
                        if not other_possible_match:
                            # filters if pair id is longer and sample_id is a substring
                            if self.pair_relationships:
                                pair_in_name = any([pair_id in identifier for pair_id in [pair['pairId'] \
                                                                                          for pair in self.pair_relationships \
                                                                                          if len(pair['pairId']) > len(sample_id)]])
                            # filters if group id is longer and sample_id is a substring
                            if self.analysis_groups:
                                group_in_name = any([group_id in identifier for group_id in [group['analysisGroupId'] \
                                                                                             for group in self.analysis_groups \
                                                                                             if len(group['analysisGroupId']) > len(sample_id)]])
                            if not pair_in_name or group_in_name:
                                matched_sample_ids.append(sample_id)
        return list(set(matched_sample_ids))
    
    def search_group_ids(self, possible_indentifiers,
                         all_matched_sample_ids):
        '''match pair_id to WDL input string'''
        matched_group_ids = []
        if self.analysis_groups:
            analysis_group_ids = [group['analysisGroupId'] for group in self.analysis_groups]
            for analysis_group in self.analysis_groups:
                group_id = analysis_group['analysisGroupId']
                #look for longer group id that the current group id may be a substring of
                other_file_group_ids = [other_group_id for other_group_id in analysis_group_ids if
                                         (not other_group_id == group_id)]
                other_file_group_ids = [other_group_id for other_group_id in other_file_group_ids
                                         if group_id in other_group_id]
                self.other_file_group_ids = list(other_file_group_ids)
                self.other_file_group_ids += [other_group_id + '.' for other_group_id in other_file_group_ids]
                self.other_file_group_ids += [other_group_id + '_' for other_group_id in other_file_group_ids]
                self.other_file_group_ids += [other_group_id + '/' for other_group_id in other_file_group_ids]
                # look to see if all samples in the analysis group have been identified in the inputs
                for identifier in possible_indentifiers:
                    if set(analysis_group['sampleIds']).issubset(set(all_matched_sample_ids)) \
                            or self.id_matches(identifier, group_id):
                        other_possible_match = any([self.id_matches(identifier, id) for id in self.other_file_group_ids])
                        if not other_possible_match:
                            matched_group_ids.append(group_id)
        return list(set(matched_group_ids))
    
    def search_pair_ids(self, possible_indentifiers,
                        matched_sample_ids):
        '''match pair_id to WDL input string'''
        all_matched_sample_ids = matched_sample_ids
        matched_pair_ids = []
        if self.pair_relationships:
            pair_ids = [pair['pairId'] for pair in self.pair_relationships]
            for pair_relationship in self.pair_relationships:
                pair_id = pair_relationship['pairId']
                #look for longer pair id that the current pair id may be a substring of
                other_file_pair_ids = [other_pair_id for other_pair_id in pair_ids if
                                         (not other_pair_id == pair_id)]
                other_file_pair_ids = [other_pair_id for other_pair_id in other_file_pair_ids
                                         if pair_id in other_pair_id]
                self.other_file_pair_ids = list(other_file_pair_ids)
                self.other_file_pair_ids += [other_pair_id + '.' for other_pair_id in other_file_pair_ids]
                self.other_file_pair_ids += [other_pair_id + '_' for other_pair_id in other_file_pair_ids]
                self.other_file_pair_ids += [other_pair_id + '/' for other_pair_id in other_file_pair_ids]
                
                for identifier in possible_indentifiers:
                    # if both ids in a pair are in the command or if the pair
                    # or if the id is in the string
                    if (pair_relationship['tumor'] in matched_sample_ids \
                            and pair_relationship['normal'] in matched_sample_ids) \
                            or self.id_matches(identifier, pair_id):
                        # but no other id is also in the input string
                        other_possible_match = any([self.id_matches(identifier, id) for id in self.other_file_pair_ids])
                        if not other_possible_match:
                                matched_pair_ids.append(pair_id)
                                all_matched_sample_ids.append(pair_relationship['tumor'])
                                all_matched_sample_ids.append(pair_relationship['normal'])
        return list(set(matched_pair_ids)), all_matched_sample_ids
 
    def search_ids(self, possible_indentifiers):
        ''' search for matches to:
        multi-patient analysis group, 
        patient-centered cohort, 
        pair, 
        or sample_id '''
        matched_sample_ids = self.search_sample_ids(possible_indentifiers)
        matched_pair_ids, all_matched_sample_ids = self.search_pair_ids(possible_indentifiers, matched_sample_ids)
        matched_group_ids = self.search_group_ids(possible_indentifiers, all_matched_sample_ids)
        if len(matched_group_ids) > 0 :
            if len(matched_group_ids) > 1:
                log.warning('multiple group ids matched the input: ' + ' '.join(matched_group_ids))
                return 'Multi sample'
            return matched_group_ids[0]
        elif len(matched_pair_ids) > 0 :
            if len(matched_pair_ids) > 1:
                log.warning('multiple pair ids matched the input: ' + ' '.join(matched_pair_ids))
                return 'Multi sample'
            return matched_pair_ids[0]
        elif len(matched_sample_ids) > 0 :
            if len(matched_sample_ids) > 1:
                log.warning('multiple sample ids matched the input: ' + ' '.join(matched_sample_ids))
                return 'Multi sample'
            return matched_sample_ids[0]
        return 'Unknown'

    def prep_ids(self):
        ''' gather sample relationship information:
        multi-patient analysis groups, 
        patient-centered cohorts, 
        pairs, 
        or sample_ids
        '''
        self.pair_relationships = False
        self.analysis_groups = False
        self.sample_ids = False
        if 'listOfPairRelationships' in self.output_info['project_data'] \
                and len(self.output_info['project_data']['listOfPairRelationships']) > 0:
            self.pair_relationships = self.output_info['project_data']['listOfPairRelationships']
        if 'analysisGroups' in self.output_info['project_data'] \
                and len(self.output_info['project_data']['analysisGroups']) > 0:
            self.analysis_groups = self.output_info['project_data']['analysisGroups']
        if 'sampleIds' in self.output_info['project_data'] \
                and len(self.output_info['project_data']['sampleIds']) > 0:
            self.sample_ids = self.output_info['project_data']['sampleIds']

    def match_id_str(self, inputs):
        ''' convert inputs to dictionary (problematic but current best options?)
        '''
        possible_indentifiers = []
        inputs = inputs.replace(" None}", " 'None'}").replace("'}  {'", "'},  {'").replace("'", '"')
        inputs = inputs.replace('["', '[').replace('"]', ']')
        cleaned_inputs = []
        for input in inputs.split():
            if '\\"' not in input:
                cleaned_inputs.append(input)
            else:
                if input.endswith('"},'):
                    cleaned_inputs.append('""},')
                elif input.endswith('"],'):
                    cleaned_inputs.append('""],')
        inputs = ' '.join(cleaned_inputs)
        if inputs != '':
            inputs_dicts = json.loads(' '.join(cleaned_inputs))
            for inputs_dict in inputs_dicts:
                if inputs_dict['type'] == 'string':
                    possible_indentifiers.append(inputs_dict['value'])
            if len(possible_indentifiers) > 0:
                top_match = self.search_ids(possible_indentifiers)
                return top_match
        return 'No possible identifiers'
    
    def get_possible(self, dictionary):
        '''search for input strings that make reference 
            to sample id name'''
        if isinstance(dictionary, dict):
            for key in dictionary:
                value = dictionary[key]
                if isinstance(value, str):
                    self.possible_indentifiers.append(value)
                elif isinstance(value, dict):
                    self.get_possible(value)
                elif isinstance(value, list):
                    for pos in value:
                        if isinstance(pos, str):
                            self.possible_indentifiers.append(pos)
                        elif isinstance(pos, dict):
                            self.get_possible(pos)
    
    def match_id(self, inputs):
        ''' convert inputs to dictionary (problematic but current best options?)
        '''
        self.possible_indentifiers = []
        self.get_possible(inputs)
        if len(self.possible_indentifiers) > 0:
            top_match = self.search_ids(self.possible_indentifiers)
            return top_match
        return 'No possible identifiers'

    
def get_args():
    '''Parse input flags
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--metrics-file-prefix',
                        help='Prefix for output metrics file',
                        required=True
                        )
    parser.add_argument('--output-info-file',
                        help='JSON file with outputInfo. Includes workflow uuid, '
                        'options JSON dictionary submitted to cromwell, '
                        'project pairing, sample, '
                        'genome build, library, output files, '
                        'subworkflow uuids, id sorted output files'
                        'and interval list information',
                        required=False
                        )
    parser.add_argument('--large-run',
                        help='Query one subworkflow at a time to avoid timeouts.',
                        required=False,
                        action='store_true'
                        )
    parser.add_argument('--multi',
                        help='Records for multiple runs are stored in the JSON file with outputInfo.',
                        required=False,
                        action='store_true'
                        )
    parser.add_argument('--name',
                        default='nygc_pipeline',
                        help='Use with the --multi flag to name the output.',
                        required=False
                        )
    parser.add_argument('--manifest',
                        help='Use with the --multi flag to list outputInfo JSON files.',
                        required=False
                        )
    parser.add_argument('--gcp-project',
                        help='GCP project id if this is not the set project '
                        'assumes you have read only access',
                        default=False,
                        required=False
                        )
    parser.add_argument('--url',
                        help='URL for cromwell server.',
                        required=True
                        )
    args_namespace = parser.parse_args()
    return args_namespace.__dict__


def main():
    args = get_args()
    output_infos = []
    manifest_files = []
    if args['multi']:
        with open(args['manifest']) as manifest:
            for line in manifest:
                output_infos.append(Runtime.read(line.strip()))
                manifest_files.append(line.strip())
    else:
        output_info = Runtime.read(args['output_info_file'])
        output_infos = [output_info]
        manifest_files.append(args['output_info_file'])
    non_retried_metrics = []
    metrics = []
    final_manifest_files = []
    for i, output_info in enumerate(output_infos):
        runtime = Runtime(limit=100000,
                          url=args['url'],
                          output_info=output_info,
                          large_run=args['large_run'],
                          metrics_file_prefix=args['metrics_file_prefix'],
                          gcp_project=args['gcp_project'])
        if runtime.loaded:
            non_retried_metrics.append(runtime.non_retried_metrics_file)
            metrics.append(runtime.metrics_file)
            final_manifest_files.append(manifest_files[i])
    if args['multi']:
        manifest = pd.DataFrame({'non_retried_metrics' : non_retried_metrics,
                                 'metrics' : metrics,
                                 'output_info' : final_manifest_files})
        manifest.to_csv(args['name'] + '_MetricsInfoManifest.csv', index=False)
    

if __name__ == "__main__":
    main()       
    