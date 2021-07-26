'''
Created on Feb 8, 2021

@author: jshelton
'''
import os
import sys
import json
import argparse
import pandas as pd
import subprocess
import logging as log
import Colorer
from collections.abc import Iterable
import copy
import meta

log.basicConfig(format='%(levelname)s:  %(message)s', level=log.INFO)


class JsonModify():
    
    def __init__(self, inputs, 
                 sample_id=False, 
                 pair_id=False,
                 dir=False,
#                  pair_ids=False,
                 all_pair_ids=False,
                 all_sample_ids=False,
                 uri_list=False,
                 tool='match'):
        '''
        deep copy required to modify complex outputs structure while
        keeping an original
        tool options:
        
        skip: skip objects that have no relevant files
        convert_path: swaps path betore the task call- directory out with
        the final output directory
        match: Assumes that the sample name followed by "." or "_" are in the basename
            Sorts out which objects have names that refer to a pair or sample
        '''
        self.tool = tool
        self.sample_id = sample_id
        self.pair_id = pair_id
        self.uri_list = uri_list
        if not self.uri_list:
            self.json_obj = copy.deepcopy(inputs)
            self.final_json_obj = copy.deepcopy(inputs)
        self.dir = dir
        self.files = []
        self.sub_workflow_uuids = []
        # create list of all other pair/sample ids in project
        # that could also match current
        self.other_file_pair_ids = []
        if all_pair_ids:
            if self.pair_id:
                other_file_pair_ids = [pair_id for pair_id in all_pair_ids if
                                       (not pair_id == self.pair_id)]
                other_file_pair_ids = [pair_id for pair_id in other_file_pair_ids
                                       if len(pair_id) > len(self.pair_id)]
                other_file_pair_ids = [pair_id for pair_id in other_file_pair_ids
                                       if self.pair_id in pair_id]
            else:
                other_file_pair_ids = all_pair_ids
            self.other_file_pair_ids = [pair_id + '.' for pair_id in other_file_pair_ids]
            self.other_file_pair_ids += [pair_id + '_' for pair_id in other_file_pair_ids]
            self.other_file_pair_ids += [pair_id + '/' for pair_id in other_file_pair_ids]
        self.other_file_sample_ids = []
        if all_sample_ids:
            other_file_sample_ids = [sample_id for sample_id in all_sample_ids if
                                   (not sample_id == self.sample_id)]
            other_file_sample_ids = [sample_id for sample_id in other_file_sample_ids
                                   if len(sample_id) > len(self.sample_id)]
            other_file_sample_ids = [sample_id for sample_id in other_file_sample_ids
                                   if self.sample_id in sample_id]
            
            self.other_file_sample_ids = [sample_id + '.' for sample_id in other_file_sample_ids]
            self.other_file_sample_ids += [sample_id + '_' for sample_id in other_file_sample_ids]
            self.other_file_sample_ids += [sample_id + '/' for sample_id in other_file_sample_ids]
        # sanity checks
        if self.tool == 'convert_path':
            assert self.dir, 'dir is required when tool is convert_path'
        elif self.tool in ['match']:
            assert self.sample_id or self.pair_id, 'pair_id or sample_id is required when tool is match'
        elif self.tool in ['skip', 'sub_workflow_uuid']:
            pass
        else:
            log.error('tool options are match, skip, sub_workflow_uuid, and convert_path: ' + tool)
        # gather all uris in the json_obj
        # or process the json object with choosen tool
        if not self.uri_list:
            self.hlpr_fnc(self.json_obj, self.final_json_obj)
        else:
            # filter uris set to false
            self.results = [self.run_tool(uri) for uri in self.uri_list if uri]
        # skip objects in output with no "true" files (because they 
        # were filtered based on sample or pair)
        if self.tool == 'skip':
            self.skip_obj = len(self.files) == 0
    
    def skip(self, string):
        '''skip objects that have no relevant files'''
        if not string.startswith('gs://'):
            return False, string
        else:
            uri = string
            basename = os.path.basename(uri)
            self.files.append(uri)
            return False, uri
        
    def match(self, string):
        '''Find objects that match a given sample or pair id'''
        if not string.startswith('gs://'):
            return string
        else:
            uri = string
            basename = os.path.basename(uri)
            # pair match search
            if not self.sample_id:
                if self.pair_id + '.' in basename or self.pair_id + '_' in basename:
                    pair_in_name = any([file_pair_id in basename for file_pair_id in self.other_file_pair_ids])
                    if not pair_in_name:
                        return uri
                elif self.pair_id + '.' in uri or self.pair_id + '_' in uri:
                    pair_in_name = any([file_pair_id in uri for file_pair_id in self.other_file_pair_ids])
                    if not pair_in_name:
                        return uri
            # sample match search
            else:
                if not any([file_pair_id in basename for file_pair_id in self.other_file_pair_ids]):
                    if self.sample_id + '.' in basename or self.sample_id + '_' in basename:
                        
                        pair_in_name = any([file_pair_id in basename for file_pair_id in self.other_file_pair_ids])
                        sample_in_name = any([file_sample_id in basename for file_sample_id in self.other_file_sample_ids])
                        if not pair_in_name and not sample_in_name:
                            return uri
                    elif self.sample_id + '.' in uri or self.sample_id + '_' in uri:
                        pair_in_name = any([file_pair_id in uri for file_pair_id in self.other_file_pair_ids])
                        sample_in_name = any([file_sample_id in uri for file_sample_id in self.other_file_sample_ids])
                        if not pair_in_name and not sample_in_name:
                            return uri
        return False
        
    def convert_path(self, string):
        '''Convert from intermediate file URI to final URI of the final output'''
        if not string.startswith('gs://'):
            return string
        else:
            uri = string
            basename = '/'.join(uri.replace('/cacheCopy', '').split('/call-')[-1].split('/')[1:])
            new_uri = self.dir + '/' + basename
            new_uri = '/'.join([section for section in new_uri.split('/') 
                               if not section.startswith('attempt-') ])
            self.files.append(new_uri)
            return new_uri

    def get_sub_workflow_uuid(self, string):
        '''find task uuid for geting metrics from big query'''
        if not string.startswith('gs://'):
            pass
        else:
            uri = string
            sub_workflow_uuid = uri.replace('/cacheCopy', '').split('/call-')[-2].split('/')[-1]
            self.sub_workflow_uuids.append(sub_workflow_uuid)
        return string
    
    def run_tool(self, subnode):
        '''run the correct tool for the string (subnode)'''
        if self.tool == 'convert_path':
            result = self.convert_path(subnode)
            return result
        elif self.tool == 'match':
            result = self.match(subnode)
            return result
        elif self.tool == 'skip':
            skip, result = self.skip(subnode)
            return result
        elif self.tool == 'sub_workflow_uuid':
            self.get_sub_workflow_uuid(subnode)
            return subnode
        
    def hlpr_fnc_list(self, list1, description):
        '''replace for list
            hlpr_fnc(list_obj, final_json_obj)
        '''
        not_matched = []
        for i, subnode in enumerate(list1):
            if isinstance(subnode, dict):
                self.hlpr_fnc(subnode, description[i])
            elif not isinstance(subnode, str) and isinstance(subnode, Iterable):
                self.hlpr_fnc_list(subnode, description[i])
            else:
                if isinstance(subnode, str):
                    result = self.run_tool(subnode)
                    if result:
                        description[i] = result
                    if self.tool in ['match']:
                        if not result:
                            not_matched.append(i)
        if self.tool in ['match']:
            offset = 0
            for ix, subnodex in enumerate(list1):
                if ix in not_matched:
                    del description[ix - offset]
                    offset += 1

    def hlpr_fnc(self, dict1, description):
        '''replace for dictionary
            hlpr_fnc(json_obj, final_json_obj)
        '''
        nodes = dict1.keys()
        not_matched = []
        for node in nodes:
            subnode = dict1[node]
            if isinstance(subnode, dict):
                self.hlpr_fnc(subnode, description[node])
            elif not isinstance(subnode, str) and isinstance(subnode, Iterable):
                self.hlpr_fnc_list(subnode, description[node])
            else:
                if isinstance(subnode, str):
                    result = self.run_tool(subnode)
                    if result:
                        description[node] = result
                    if self.tool in ['match']:
                        if not result:
                            not_matched.append(node)
        if self.tool in ['match']:
            for n in not_matched:
                del description[n]
                    
                    
                    
class CloudOutput():
    def __init__(self, run_data, url):
        self.parent_dir = os.path.abspath(os.path.dirname(__file__))
        self.url = url
        self.run_data = run_data
        self.options = run_data.options
        self.workflow_uuid = self.run_data.workflow_uuid
        # fast get calls from API to search for uuids
        self.uri_list = self.list_uris()
        if len(self.uri_list) == 0:
            # slow get uris from bucket directory
            log.warning('Slow uri lookup beginning')
            self.uri_list = self.list_project()
        # get API outputs section
        self.raw_outputs = self.read_api()
        # get root dir
        self.root_dir = self.read_api(goal='get_root')
        # get files (named and unnamed) from API outputs section
        self.parse_block()
        self.unnamed_files = self.gather_unnamed()
        # Find only files with the pair or the sample in the 
        # filename with . or _ after the name
        self.divide_by_id()
        # Find only task uuids for files with the pair or the sample in the 
        # filename with . or _ after the name
        self.divide_uuid_by_id()
        self.run_data.run_info['pair_association'] = self.pair_association
        self.run_data.run_info['sample_association'] = self.sample_association
        self.run_data.run_info['outputs'] = self.named_outputs
        self.run_data.run_info['named_files'] = self.named_files
        self.run_data.run_info['unnamed_files'] = self.unnamed_files
        self.run_data.run_info['sub_workflow_uuids'] = self.sub_workflow_uuids
        
    def filter_by_id(self, association, pair=True, sample=False):
        '''skip object that have no files after filtering by id'''
        association = copy.deepcopy(association)
        new_association = {}
        for id in association:
            if not id in new_association:
                new_association[id] = {}
            for object in association[id]:
                association_object = copy.deepcopy(association[id][object])
                if pair:
                    skip = JsonModify(inputs={object: association_object},
                                      pair_id=id,
                                      tool='skip')
                if sample:
                    skip = JsonModify(inputs={object: association_object},
                                      sample_id=id,
                                      tool='skip')
                if not skip.skip_obj:
                    new_association[id][object] = skip.final_json_obj[object]
        return new_association
    
    def divide_by_id(self):
        ''' Find only files with the pair or the sample in the 
        filename with . or _ after the name'''
        pair_association = {}
        sample_association = {}
        pair_ids = list(set([pair_info["pairId"] for pair_info in self.run_data.project_info["pairInfos"]]))
        for pair_id in pair_ids:
            current_pair_association = JsonModify(inputs=self.named_outputs, 
                                                    all_pair_ids=list(pair_ids),
                                                    pair_id=pair_id,
                                                    tool='match').final_json_obj
            pair_association[pair_id] = current_pair_association
        for sample_id in self.run_data.project_info["sampleIds"]:
            current_sample_association = JsonModify(inputs=self.named_outputs,
                                                       all_sample_ids=self.run_data.project_info["sampleIds"],
                                                       sample_id=sample_id,
                                                       all_pair_ids=list(pair_ids),
                                                       tool='match').final_json_obj
            sample_association[sample_id] = current_sample_association
        self.pair_association = self.filter_by_id(pair_association, 
                                                  pair=True, 
                                                  sample=False)
        self.sample_association = self.filter_by_id(sample_association, 
                                                    pair=False, 
                                                    sample=True)
        for pair_id in self.pair_association:
            for object in self.pair_association[pair_id]:
                if isinstance(self.pair_association[pair_id][object], dict):
                    self.pair_association[pair_id][object] = {self.pair_association[pair_id][object][key] for key in self.pair_association[pair_id][object]
                                                              if not len(self.pair_association[pair_id][object][key]) == 0}
                elif not isinstance(self.pair_association[pair_id][object], str) \
                        and isinstance(self.pair_association[pair_id][object], Iterable):
                    self.pair_association[pair_id][object] = [self.pair_association[pair_id][object][i] for i in range(len(self.pair_association[pair_id][object]))
                                                              if not len(self.pair_association[pair_id][object][i]) == 0]
        
    def divide_uuid_by_id(self):
        ''' Find only task uuids for files with the pair or the sample in the 
        filename with . or _ after the name'''
        self.sub_workflow_uuids = {}
        pair_ids = list(set([pair_info["pairId"] for pair_info in self.run_data.project_info["pairInfos"]]))
        for pair_id in pair_ids:
            matches = JsonModify(inputs={},
                                 uri_list=self.uri_list,
                                 all_pair_ids=list(pair_ids),
                                 pair_id=pair_id,
                                 tool='match').results
            sub_workflow_uuids = JsonModify(inputs={},
                                    uri_list=matches,
                                    all_pair_ids=list(pair_ids),
                                    pair_id=pair_id,
                                    tool='sub_workflow_uuid').sub_workflow_uuids
            self.sub_workflow_uuids[pair_id] = list(set([uuid for uuid in sub_workflow_uuids if uuid]))
        for sample_id in self.run_data.project_info["sampleIds"]:
            matches = JsonModify(inputs={},
                                 uri_list=self.uri_list,
                                 all_pair_ids=list(pair_ids),
                                 sample_id=sample_id,
                                 tool='match').results
            sub_workflow_uuids = JsonModify(inputs={},
                                    uri_list=matches,
                                    all_pair_ids=list(pair_ids),
                                    sample_id=sample_id,
                                    tool='sub_workflow_uuid').sub_workflow_uuids
            self.sub_workflow_uuids[sample_id] = list(set([uuid for uuid in sub_workflow_uuids if uuid]))

    def out_path(self):
        if self.options["use_relative_output_paths"]:
            dir = self.options["final_workflow_outputs_dir"]
            return dir
        else:
            log.error('Not tested for non-relative paths')
            sys.exit(1)
        
    def parse_block(self):
        '''get files (named and unnamed) from API outputs section'''
        self.out_dir =  self.out_path()
        outputs = JsonModify(inputs=self.raw_outputs, 
                             dir=self.out_dir,
                             tool='convert_path')
        self.named_outputs = outputs.final_json_obj
        self.named_files = outputs.files
        
    def gather_unnamed(self):
        '''Some (generally undelivered outputs) have a list of files
        too long for the log file'''
        all_files = self.ls_gsutil()
        unnamed_files = list(set(all_files).difference(set(self.named_files)))
        return unnamed_files
    
    def list_uris(self):
        '''
        Get all top level sub workflow uuids
        '''
        calls = self.read_api(goal='get_uuids')
        self.uri_list = []
        for line in str(calls).split('\n'):
            for strings in line.split(' '):
                for string in strings.strip().split("'"):
                    if string.startswith("gs://") and 'call-' in string:
                        uri = string
                        self.uri_list.append(uri)
                    elif string.startswith("gs://"):
                        log.warning('potential skipped uri : ' + string)
        self.uri_list = list(set(self.uri_list))
        return list(set(self.uri_list))
    
    def list_project(self):
        '''slow way to get all uris
            
        '''
        skip = [':', '/',
                 '/gcs_delocalization.sh', '/rc', 
                '/gcs_localization.sh', '/gcs_transfer.sh',
                '/script', '/stderr', '/stdout']
        log.info('Retrieve uuids for run...')
        try:
            uris = subprocess.run(['gsutil', 'ls', 
                                   '-R',
                                   self.root_dir
                                   ],
                                   check=True,
                                   stdout=subprocess.PIPE).stdout.decode('utf-8')
        except subprocess.CalledProcessError as err:
            log.error(err.output.decode('utf-8'))
            log.error('Failed to recursively list directories')
            return False
        final_uris = []
        for uri in [uri for uri in uris.split('\n') if not uri == '']:
            failed = any([uri.endswith(end) for end in skip])
            if not failed:
                final_uris.append(uri)
        return final_uris
       
    def ls_gsutil(self):
        '''
            return empty list if no URLs are found
        '''
        try:
            result = subprocess.run(['gsutil', 'ls', self.out_dir + '/*'],
                                    check=True,
                                    stdout=subprocess.PIPE).stdout.decode('utf-8')
        except subprocess.CalledProcessError as err:
            log.warning(err.output.decode('utf-8'))
            log.warning('Failed to read URI: ' + self.out_dir)
            return []
        return [file for file in result.split('\n') if not file == '' 
                and not file.endswith(':')
                and not file.endswith('/')]
                    
    def read_api(self, goal='get_outputs'):
        if goal == 'get_outputs':
            script = self.parent_dir + '/get_outputs.sh'
        if goal == 'get_root':
            script = self.parent_dir + '/get_root.sh'
        if goal == 'get_uuids':
            script = self.parent_dir + '/get_uuids.sh'
        try:
            result = subprocess.run(['bash', script,
                                     '-u', self.url,
                                     '-w', self.workflow_uuid],
                                    check=True,
                                    stdout=subprocess.PIPE).stdout.decode('utf-8')
        except subprocess.CalledProcessError as err:
            log.error(err.output.decode('utf-8'))
            log.error('Failed to get outputs from api for workflow: ' + self.workflow_uuid)
            return False
        return json.loads(result)
        
        
        
class RunData():
    '''Load rundata from file and confirm project data schema '''
    def __init__(self, run_info_file):
        self.run_info = self.read(run_info_file)
        self.run_date = run_info_file.split('.')[-3:-2][0]
        self.options = self.run_info['project_data']['options']
        self.project_info = self.run_info['project_data']
        self.project_info['run_date'] = self.run_date
        self.passed = meta.test_schema(self.project_info)
        self.workflow_uuid = self.run_info['workflow_uuid']
        
    def read(self, file):
        with open(file) as project_info_file:
            project_info = json.load(project_info_file)
            return project_info
        

def main():
    args = get_args()
    run_data = RunData(run_info_file=args['run_data'])
    outputs = CloudOutput(run_data=run_data, url=args['url'])
    file_out = outputs.run_data.run_info["project_data"]['project'].replace(' ', '_') + '.' +  outputs.run_data.run_info['workflow_uuid'] + '_outputInfo.json'
    with open(file_out, 'w') as project_info_file:
        json.dump(outputs.run_data.run_info, project_info_file, indent=4)

        

            
def get_args():
    '''Parse input flags
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-data',
                        help='JSON file with RunInfo. Includes workflow uuid, '
                        'options JSON dictionary submitted to cromwell, '
                        'project pairing, sample, '
                        'genome build, library and interval list information',
                        required=True
                        )
    parser.add_argument('--url',
                        help='URL for cromwell server.',
                        required=True
                        )
    args_namespace = parser.parse_args()
    return args_namespace.__dict__

if __name__ == "__main__":
    main()       
