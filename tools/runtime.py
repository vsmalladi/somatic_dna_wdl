import sys
from google.cloud import bigquery
from google.cloud import bigquery_storage
import json
import datetime
import logging as log
import pandas as pd
import Colorer
import make_auth
import argparse
import os 
import numpy as np

pd.set_option('display.max_columns', 500)
log.basicConfig(format='%(levelname)s:  %(message)s', level=log.INFO)


class Runtime():
    def __init__(self,
                 output_info,
                 limit=10000,
                 gcp_project=False,
                 metrics_file=False):
        self.output_info = output_info
        # find existing project id info
        self.prep_ids()
        if not metrics_file:
            self.metrics_file = output_info["project_data"]['project'].replace(' ', '_') + '.' +  output_info['workflow_uuid'] + '_outputMetrics.csv'
        else:
            self.metrics_file = metrics_file
        self.non_retried_metrics_file = self.metrics_file.replace('.csv', '.non_retried.csv')
        if os.path.isfile(self.metrics_file):
            self.loaded = True
        else:
            # find relevant project info
            self.limit = str(limit)
            self.metrics_limit = str(limit * 4)
            # login
            self.credentials, self.gcp_project = make_auth.login()
            if gcp_project:
                self.gcp_query_project = gcp_project
            else:
                self.gcp_query_project = self.gcp_project
            self.sub_workflow_uuids = self.output_info['sub_workflow_uuids']
            self.run_date = self.get_rundate()
            self.get_sample_info()
            # gather big cloud tables
            self.bqclient = bigquery.Client(project=self.gcp_project , 
                                            credentials=self.credentials)
            self.bqstorageclient = bigquery_storage.BigQueryReadClient(credentials=self.credentials)
            self.loaded = self.load_metadata()
            if not self.loaded:
                log.warning('No endtime found for workflow: ' + output_info['workflow_uuid'])
            else:
                self.load_runtime()
                self.start_instance_maps()
                self.load_metrics()
                self.add_to_instance_maps()
                hashable = ['instance_name', 'workflow_id', 'attempt', 'shard']
                self.metadata = self.metadata.drop_duplicates(subset=hashable)
                self.instance_id_map = self.get_max_mems(self.metadata, 
                                                         instance_id_map=self.instance_id_map, 
                                                         execution_status=False,
                                                         metrics_key='mem_used_gb')
                self.non_retry_instance_id_map = self.get_max_mems(self.metadata, 
                                                                   instance_id_map=self.non_retry_instance_id_map, 
                                                                   execution_status='Done',
                                                                   metrics_key='mem_used_gb')
                self.non_retry_metadata = self.metadata[self.metadata.execution_status == 'Done'].copy()
                if self.metadata[self.metadata.execution_status != 'RUNNING'].empty:
                    self.metadata.to_csv(self.metrics_file, index=False)
                else:
                    self.metadata = self.load_plot_metrics(self.metadata,
                                                           instance_id_map=self.instance_id_map)
                    self.non_retry_metadata = self.load_plot_metrics(self.non_retry_metadata,
                                                                     instance_id_map=self.non_retry_instance_id_map)
                    self.metadata.to_csv(self.metrics_file, index=False)
                    self.non_retry_metadata.to_csv(self.non_retried_metrics_file, index=False)
        
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
    
    def match_instance_id(self, instance_id_map, row, execution_status='Done'):
        for inst in instance_id_map:
            if (row.workflow_id in instance_id_map[inst]['workflow_ids'] \
                    and row.task_call_name == instance_id_map[inst]['task_call_name']):
                if execution_status:
                    if execution_status == row.execution_status:
                        yield inst
                else:
                    yield inst
            
    
    def get_max_mems(self, metadata, instance_id_map, 
                     execution_status='Done',
                     job_id_col='task_call_name', metrics_key='mem_used_gb'):
        '''Get max mem for a runs with the same sub-workflow_id and task_call_name'''
        runtime_results = metadata.drop_duplicates(['task_call_name', 'workflow_id', 'instance_name', 'shard', 'attempt'])
        max_values = {}
        for index, row in runtime_results.iterrows():
            instance_ids = []
            for inst in self.match_instance_id(instance_id_map, 
                                               row, 
                                               execution_status=execution_status):
                instance_ids.append(inst)
            instance_ids = list(set(instance_ids))
            try:
                max_value = self.metrics[self.metrics.instance_id.isin(instance_ids)][metrics_key].max()
            except ValueError:
                max_value = 0
            for inst in instance_ids:
                instance_id_map[inst][metrics_key] = max_value
        return instance_id_map
    
    def get_max_mem(self, grouped, row, ids, 
                    instance_id_map,
                    metadata, metrics_key='mem_used_gb'):
        '''Get max mem for a given subset of the runs'''
        runtime_results = metadata.copy()
        for id in ids:
            runtime_results = runtime_results[runtime_results[id] == row[id]].copy()
        max_values = []
        for inst in instance_id_map:
            if (row.workflow_id in instance_id_map[inst]['workflow_ids']) \
                    and (row.task_call_name == instance_id_map[inst]['task_call_name']):
                max_values.append(instance_id_map[inst][metrics_key])
        try:
            return max(max_values)
        except ValueError:
            return 0
        
    def get_flow_runtime(self, grouped, row, ids):
        ''' get max runtime in hours for group.
        For sub workflow and workflows this will also include wait time and
        return the end-to-end wallclock time for a particular id.
        For tasks this will only return the longest runtime.
        '''
        if 'task_call_name' in ids:
            runtime_results = grouped.run_time_m.agg(['max']).reset_index()
            for id in ids:
                runtime_results = runtime_results[runtime_results[id] == row[id]].copy()
            runtime = runtime_results['max'].tolist()[0]
            return runtime / 60.0
        else:
            end_time_results = grouped.end_time.agg(['max']).reset_index()
            for id in ids:
                end_time_results = end_time_results[end_time_results[id] == row[id]].copy()
            end_time = end_time_results['max'].tolist()[0]
            start_time_results = grouped.start_time.agg(['min']).reset_index()
            for id in ids:
                start_time_results = start_time_results[start_time_results[id] == row[id]].copy()
            start_time = start_time_results['min'].tolist()[0]
            runtime = end_time - start_time
        return runtime.total_seconds() / 60.0 / 60.0
    
    def get_from_map(self, instance_id_map, row, query):
        try:
            instance_id = [instance_id for instance_id in instance_id_map 
                           if instance_id_map[instance_id]['instance_name'] == row.instance_name][0]
        except IndexError:
            return None
        return instance_id_map[instance_id][query]
        
    def load_plot_metrics(self, metadata, instance_id_map):
        '''add cpu_time, max_mem, wall_clock, cpu_used_percent, mem_efficiency, mem_used_gb
        also add cpu_used_percent", "disk_total_gb",  "disk_used_gb",  "cpu_count",   "cpu_platform", "preemptible"
        
        '''
        metadata['max_cpu_used_percent'] = metadata.apply(lambda row: self.get_from_map(instance_id_map, 
                                                                                    row, 
                                                                                    query='max_cpu_used_percent'), axis=1)
        metadata['disk_total_gb'] = metadata.apply(lambda row: self.get_from_map(instance_id_map, 
                                                                                    row, 
                                                                                    query='disk_total_gb'), axis=1)
        metadata['disk_used_gb'] = metadata.apply(lambda row: self.get_from_map(instance_id_map, 
                                                                                    row, 
                                                                                    query='disk_used_gb'), axis=1)
        metadata['cpu_count'] = metadata.apply(lambda row: self.get_from_map(instance_id_map, 
                                                                                    row, 
                                                                                    query='cpu_count'), axis=1)
        metadata['cpu_platform'] = metadata.apply(lambda row: self.get_from_map(instance_id_map, 
                                                                                    row, 
                                                                                    query='cpu_platform') , axis=1)
        metadata['preemptible'] = metadata.apply(lambda row: self.get_from_map(instance_id_map, 
                                                                                    row, 
                                                                                    query='preemptible'), axis=1)
   
        # get wallclock
        metadata['start_time'] = pd.to_datetime(metadata.start_time, infer_datetime_format=True)
        metadata['end_time'] = pd.to_datetime(metadata.end_time, infer_datetime_format=True)
        metadata['run_time'] = metadata['end_time'] - metadata['start_time']
        metadata['run_time_m'] = metadata.apply(lambda row:
                                                     (row.run_time.total_seconds() / 60.0), axis=1)
        # get cpu_time
        metadata['cpu_time_m'] = metadata.apply(lambda row:
                                                (row.run_time_m
                                                 * row.cpu_count),
                                                axis=1)
        # ==============
        #     Tasks
        # ==============
        grouped = metadata.groupby('task_call_name')
        metadata['mean_task_core_h'] = metadata.apply(lambda row:
                                                           (row.cpu_count
                                                            * float(grouped.mean().run_time_m[row.task_call_name])) / 60,
                                                           axis=1)
        # task wallclock time from start to finish
        metadata['mean_task_run_time_h'] = metadata.apply(lambda row:
                                                          float(grouped.mean().run_time_m[row.task_call_name]) / 60,
                                                          axis=1)
        
        # ============================
        #     Tasks w/in subworkflow
        # ============================
        grouped = metadata.groupby(['id', 'task_call_name', 'workflow_name', 'instance_name'])
        # sub task wallclock time from start to finish
        metadata['sample_task_run_time_h'] = metadata.apply(lambda row: 
                                                            self.get_flow_runtime(grouped, row,
                                                                                  ids=['id', 
                                                                                       'instance_name',
                                                                                       'task_call_name',
                                                                                       'workflow_name']),
                                                            axis=1)
        metadata['max_mem_g'] = metadata.apply(lambda row: self.get_max_mem(grouped, row,
                                                                            metadata=metadata,
                                                                            metrics_key='mem_used_gb',
                                                                            instance_id_map=instance_id_map,
                                                                            ids=['id',
                                                                                 'instance_name',
                                                                                 'task_call_name',
                                                                                 'workflow_name']),
                                                                                 axis=1)
        metadata['sample_task_core_h'] = metadata.apply(lambda row:
                                                        self.get_task_core_h(grouped, row,
                                                                             ids=['id',
                                                                                  'instance_name',
                                                                                  'task_call_name',
                                                                                  'workflow_name']),
                                                        axis=1)
        # ============================
        #       sub workflows
        # ============================
        grouped = metadata.groupby(['id', 'workflow_name'])
        sub_metadata = grouped.agg(sample_subworkflow_core_h=pd.NamedAgg(column='sample_task_core_h', aggfunc=sum)).reset_index()
        metadata['sample_subworkflow_core_h'] = metadata.apply(lambda row:
                                                               sub_metadata[(sub_metadata.id == row.id)
                                                                            & (sub_metadata.workflow_name == row.workflow_name)].sample_subworkflow_core_h.tolist()[0], 
                                                               axis=1)
        # sub workflows wallclock time from start to finish
        metadata['sample_subworkflow_run_time_h'] = metadata.apply(lambda row:
                                                                   self.get_flow_runtime(grouped, row,
                                                                                         ids=['id', 'workflow_name']), 
                                                                   axis=1)
        metadata['subworkflow_max_mem_g'] = metadata.apply(lambda row: metadata[(metadata.id == row.id)
                                                                                & (metadata.workflow_name == row.workflow_name)].max_mem_g.max(),
                                                                                axis=1)
        # ============================
        #       workflows
        # ============================
        # full workflow wallclock time from start to finish for an id
        grouped = metadata.groupby('id')
        sub_metadata = grouped.agg(sample_workflow_core_h=pd.NamedAgg(column='sample_task_core_h', aggfunc=sum)).reset_index()
        metadata['sample_workflow_core_h'] = metadata.apply(lambda row: 
                                                            sub_metadata[(sub_metadata.id == row.id)].sample_workflow_core_h.tolist()[0],
                                                            axis=1)
        metadata['sample_workflow_run_time_h'] = metadata.apply(lambda row: 
                                                                self.get_flow_runtime(grouped, row, 
                                                                                      ids=['id']), 
                                                                axis=1)
        metadata['workflow_max_mem_g'] = metadata.apply(lambda row: metadata[(metadata.id == row.id)].max_mem_g.max(),
                                                        axis=1) 
        return metadata
    
    def load_end_date(self, date_record):
        try:
            return date_record.strftime(format='%Y-%m-%d')
        except ValueError:
            return False     

    def modify_date(self, run_date):
        return '-'.join(run_date.split('-')[0:3])
    
    def get_sample_info(self):
        self.sample_count = len(self.output_info['project_data']['sampleIds'])
        self.pair_count = len(self.output_info['project_data']['pairId'])
        
    def get_rundate(self):
        return self.modify_date(self.output_info['project_data']['run_date']) 
        
    def get_gcp_project(self):
        return self.output_info['project_data']['options']['monitoring_image'].split('/')[1]

    def format_sub_workflow(self, sub_workflow_uuids):
        sub_workflow_uuid_list = '("' + '", "'.join(sub_workflow_uuids) + '")'
        return sub_workflow_uuid_list
    
    def load_metadata(self):
        '''convert big query metadata table into pandas dataframe
            Includes start_time, end_time, execution_status, cpu_count,
            mem_total_gb, task_call_name, workflow_name
        '''
        log.info('Loading Metadata...')
        
        self.sub_workflow_uuid_list = self.format_sub_workflow(self.sub_workflow_uuids)
        
        query_string = '''
        SELECT * FROM `''' + self.gcp_query_project + '''.cromwell_monitoring.metadata`
            WHERE DATE(start_time) >= "''' + self.run_date + '''"
            AND workflow_id IN ''' + self.sub_workflow_uuid_list
            
#             + '''"
#             LIMIT ''' + self.limit + '''
#         '''
        self.metadata =  (
            self.bqclient.query(query_string)
            .result()
            .to_dataframe(bqstorage_client=self.bqstorageclient)
        )
        self.end_date = self.load_end_date(self.metadata.end_time.max())
        if not self.end_date:
            return False
        self.metadata['inputs'] = self.metadata.apply(lambda row: str(row.inputs).replace('\n', ' '), axis=1)
        self.metadata['id'] = self.metadata.apply(lambda row: self.match_id(row.inputs).replace('\n', ' '), axis=1)
        unhashable = ['disk_mounts', 'disk_total_gb', 'disk_used_gb', 'disk_types', 'inputs']
        hashable = ['instance_name', 'workflow_id', 'attempt', 'shard']
        self.metadata = self.metadata.drop_duplicates(subset=hashable)
        return True
        
    def load_runtime(self):
        '''convert big query runtime table into pandas dataframe
        Includes both instance_id and workflow_id.
        Use to map workflow_id to instance_ids
        '''
        log.info('Loading Runtime...')
        query_string = '''
        SELECT * FROM `''' + self.gcp_query_project + '''.cromwell_monitoring.runtime`
            WHERE DATE(start_time) >= "''' + self.run_date + '''"
            AND DATE(start_time) <= "''' + self.end_date + '''"
            AND workflow_id IN ''' + self.sub_workflow_uuid_list + '''
            LIMIT ''' + self.limit + '''
        '''
        runtime = (
            self.bqclient.query(query_string)
            .result()
            .to_dataframe(bqstorage_client=self.bqstorageclient)
        )
        self.runtime = runtime[runtime.workflow_id.isin(self.sub_workflow_uuids)].copy()
        
    def start_instance_maps(self):
        self.instance_ids = self.runtime.instance_id.unique().tolist()
        self.instance_id_map = {}
        for instance_id in self.instance_ids:
            self.instance_id_map[instance_id] = {}
            self.instance_id_map[instance_id]['workflow_ids'] = self.runtime[self.runtime.instance_id == instance_id].workflow_id.unique().tolist()[0]
            self.instance_id_map[instance_id]['task_call_name'] = self.runtime[self.runtime.instance_id == instance_id].task_call_name.tolist()[0]
            self.instance_id_map[instance_id]['instance_name'] = self.runtime[(self.runtime.instance_id == instance_id)].instance_name.tolist()[0]
        self.non_retry_instance_id_map = {}
        map = {row.instance_name: row.execution_status for index, row in self.metadata[['instance_name', 'execution_status']].drop_duplicates().iterrows()}
        execution_status_map = {row.instance_id : map[row.instance_name] for index, row in 
                                self.runtime[['instance_name', 'instance_id']].drop_duplicates().iterrows()}
        self.non_retry_instance_ids = [inst for inst in execution_status_map]
        for instance_id in self.non_retry_instance_ids:
            self.non_retry_instance_id_map[instance_id] = {}
            self.non_retry_instance_id_map[instance_id]['workflow_ids'] = self.runtime[(self.runtime.instance_id == instance_id)].workflow_id.unique().tolist()[0]
            self.non_retry_instance_id_map[instance_id]['task_call_name'] = self.runtime[(self.runtime.instance_id == instance_id)].task_call_name.tolist()[0]
            self.non_retry_instance_id_map[instance_id]['instance_name'] = self.runtime[(self.runtime.instance_id == instance_id)].instance_name.tolist()[0]
        
    def format_instance_id(self):
        instance_id_list = '(' + ', '.join([str(i) for i in self.instance_ids]) + ')'
        return instance_id_list
    
    def get_cpu_used(self, instance_id):
        metrics_instance = self.metrics[self.metrics.instance_id == instance_id].copy()
        cpu_used_percents = []
        for cpu_used_percent_list in metrics_instance.cpu_used_percent.tolist():
            for cpu_used_percent in cpu_used_percent_list:
                cpu_used_percents.append(cpu_used_percent)
        if len(cpu_used_percents) > 0:
            return max(cpu_used_percents)
        return np.nan
    
    def get_disk_used_gb(self, instance_id):
        metrics_instance = self.metrics[self.metrics.instance_id == instance_id].copy()
        disk_used_gbs = []
        for disk_used_gb_list in metrics_instance.disk_used_gb.tolist():
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
    
    def load_metrics(self):
        '''convert big query metrics table into pandas dataframe
        Includes instance_id, cpu_used_percent, mem_used_gb, disk_used_gb, 
        disk_read_iops, disk_write_iops.
        Use to find max and plot resource usage over time.
        '''
        log.info('Loading Metrics...')
        
        instance_id_list = self.format_instance_id()
        
        query_string = '''
        SELECT * FROM `''' + self.gcp_query_project + '''.cromwell_monitoring.metrics`
            WHERE DATE(timestamp) >= "''' + self.run_date + '''"
            AND DATE(timestamp) <= "''' + self.end_date + '''"
            AND instance_id IN ''' + instance_id_list 
        self.metrics = (
            self.bqclient.query(query_string)
            .result()
            .to_dataframe(bqstorage_client=self.bqstorageclient)
        )
        
    def add_to_instance_maps(self):
        for instance_id in self.instance_ids:
            max_cpu_used_percent = self.get_cpu_used(instance_id)
            disk_used_gb = self.get_disk_used_gb(instance_id)
            disk_total_gb = self.get_disk_total_gb(instance_id)
            match = self.runtime[self.runtime.instance_id == instance_id].copy()
            cpu_count = match.cpu_count.tolist()[0]
            cpu_platform = match.cpu_platform.tolist()[0]
            preemptible = match.preemptible.tolist()[0]

            self.instance_id_map[instance_id]['max_cpu_used_percent'] = max_cpu_used_percent
            self.instance_id_map[instance_id]['disk_used_gb'] = disk_used_gb
            self.instance_id_map[instance_id]['disk_total_gb'] = disk_total_gb
            self.instance_id_map[instance_id]['cpu_count'] = cpu_count
            self.instance_id_map[instance_id]['cpu_platform'] = cpu_platform
            self.instance_id_map[instance_id]['preemptible'] = preemptible
            self.instance_id_map[instance_id]['disk_total_gb'] = disk_total_gb
            self.instance_id_map[instance_id]['disk_used_gb'] = disk_used_gb
            if instance_id in self.non_retry_instance_ids:
                self.non_retry_instance_id_map[instance_id]['max_cpu_used_percent'] = max_cpu_used_percent
                self.non_retry_instance_id_map[instance_id]['disk_used_gb'] = disk_used_gb
                self.non_retry_instance_id_map[instance_id]['disk_total_gb'] = disk_total_gb
                self.non_retry_instance_id_map[instance_id]['cpu_count'] = cpu_count
                self.non_retry_instance_id_map[instance_id]['cpu_platform'] = cpu_platform
                self.non_retry_instance_id_map[instance_id]['preemptible'] = preemptible
                self.non_retry_instance_id_map[instance_id]['disk_total_gb'] = disk_total_gb
                self.non_retry_instance_id_map[instance_id]['disk_used_gb'] = disk_used_gb
                
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
                    if set(analysis_group[sampleIds]).issubset(set(all_matched_sample_ids)) \
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

    def match_id(self, inputs):
        ''' convert inputs to dictionary (problematic but current best options?)
        '''
        possible_indentifiers = []
        inputs = inputs.replace(" None}", " 'None'}").replace("'}  {'", "'},  {'").replace("'", '"')
        inputs = inputs.replace('["', '[').replace('"]', ']')
        inputs_dicts = json.loads(inputs)
        for inputs_dict in inputs_dicts:
            if inputs_dict['type'] == 'string':
                possible_indentifiers.append(inputs_dict['value'])
        if len(possible_indentifiers) > 0:
            top_match = self.search_ids(possible_indentifiers)
            return top_match
        return 'No possible identifiers'

    
def get_args():
    '''Parse input flags
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-info-file',
                        help='JSON file with outputInfo. Includes workflow uuid, '
                        'options JSON dictionary submitted to cromwell, '
                        'project pairing, sample, '
                        'genome build, library, output files, '
                        'subworkflow uuids, id sorted output files'
                        'and interval list information',
                        required=False
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
                          output_info=output_info,
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
    