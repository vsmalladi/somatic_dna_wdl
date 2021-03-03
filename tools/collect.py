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

log.basicConfig(format='%(levelname)s:  %(message)s', level=log.INFO)


import meta


class JsonUpdate():
    
    def __init__(self, inputs, dir):
        '''deep copy required to modify complex outputs structure while
        keeping an original'''
        
        self.dir = dir
        self.json_obj = copy.deepcopy(inputs)
        self.final_json_obj = copy.deepcopy(inputs)
        self.files = []
        self.hlpr_fnc(self.json_obj, self.final_json_obj)
        print(self.final_json_obj)
        
    def convert_path(self, string):
        if string.endswith('...'):
            # LongObject'
            return False
        elif not string.startswith('gs://'):
            return string
        else:
            uri = string
            basename = '/'.join(uri.split('/call-')[-1].split('/')[1:])
            self.files.append(self.dir + '/' + basename)
            return self.dir + '/' + basename
        
    def hlpr_fnc_list(self, list1, description):
        for i, subnode in enumerate(list1):
            if isinstance(subnode, dict):
                self.hlpr_fnc(subnode, description[i])
            elif not isinstance(subnode, str) and isinstance(subnode, Iterable):
                self.hlpr_fnc_list(subnode, description[i])
            else:
                if isinstance(subnode, str):
                    description[i] = self.convert_path(subnode)
                
    def hlpr_fnc(self, dict1, description):
        nodes = dict1.keys()
        for node in nodes:
            subnode = dict1[node]
            if isinstance(subnode, dict):
                self.hlpr_fnc(subnode, description[node])
            elif not isinstance(subnode, str) and isinstance(subnode, Iterable):
                self.hlpr_fnc_list(subnode, description[node])
            else:
                if isinstance(subnode, str):
                    description[node] = self.convert_path(subnode)

class JsonMatch():
    
    def __init__(self, inputs, sample_id, pair_id):
        '''Assumes that the sample name followed by "." or "_" are in the basename'''
        self.sample_id = sample_id
        self.pair_id = pair_id
        self.json_obj = copy.deepcopy(inputs)
        self.final_json_obj = copy.deepcopy(inputs)
        self.files = []
        self.hlpr_fnc(self.json_obj, self.final_json_obj)
        
    def match(self, string):
        if not string.startswith('gs://'):
            return string
        else:
            uri = string
            basename = uri.split('/')[-1]
            if self.pair_id:
                if self.pair_id + '.' in basename or self.pair_id + '_' in basename:
                    return uri
            else:
                if self.sample_id+ '.' in basename or self.sample_id + '_' in basename:
                    return uri
            return False
        
    def hlpr_fnc_list(self, list1, description):
        for i, subnode in enumerate(list1):
            if isinstance(subnode, dict):
                self.hlpr_fnc(subnode, description[i])
            elif not isinstance(subnode, str) and isinstance(subnode, Iterable):
                self.hlpr_fnc_list(subnode, description[i])
            else:
                if isinstance(subnode, str):
                    description[i] = self.match(subnode)
                
    def hlpr_fnc(self, dict1, description):
        nodes = dict1.keys()
        for node in nodes:
            subnode = dict1[node]
            if isinstance(subnode, dict):
                self.hlpr_fnc(subnode, description[node])
            elif not isinstance(subnode, str) and isinstance(subnode, Iterable):
                self.hlpr_fnc_list(subnode, description[node])
            else:
                if isinstance(subnode, str):
                    description[node] = self.match(subnode)

class CloudOutput():
    def __init__(self, run_data):
        self.run_data = run_data
        self.options = run_data.options
        self.uri = run_data.options["final_workflow_log_dir"] + '/workflow.' + self.run_data.workflow_uuid + '.log'
        self.workflow_log = self.read_gsutil()
        self.gather_outputs()
        self.unnamed_files = self.gather_unnamed()
        self.divide_by_id()
#         self.run_data.run_info['pair_association'] = self.pair_association
#         self.run_data.run_info['sample_association'] = self.sample_association
        self.run_data.run_info['outputs'] = self.named_outputs
        self.run_data.run_info['named_files'] = self.named_files
        self.run_data.run_info['unnamed_files'] = self.unnamed_files
        
    def divide_by_id(self):
        self.pair_association = {}
        self.sample_association = {}
        pair_ids = list(set([pair_info["pairId"] for pair_info in self.run_data.project_info["pairInfos"]]))
        for pair_id in pair_ids:
            self.pair_association[pair_id] = JsonMatch(inputs=self.named_outputs, 
                                                       sample_id=False, 
                                                       pair_id=pair_id).final_json_obj
        for sample_id in self.run_data.project_info["sampleIds"]:
            self.sample_association[sample_id] = JsonMatch(inputs=self.named_outputs, 
                                                           sample_id=sample_id, 
                                                           pair_id=False).final_json_obj

    def out_path(self):
        if self.options["use_relative_output_paths"]:
            dir = self.options["final_workflow_outputs_dir"]
            return dir
        else:
            log.error('Not tested for non-relative paths')
            sys.exit(1)
        
    def parse_block(self):
        self.out_dir =  self.out_path()
        outputs = JsonUpdate(inputs=self.raw_outputs, 
                              dir=self.out_dir)
        self.named_outputs = outputs.final_json_obj
        self.named_files = outputs.files
        
    def gather_unnamed(self):
        '''Some (generally undelivered outputs) have a list of files
        too long for the log file'''
        all_files = self.ls_gsutil()
        unnamed_files = list(set(all_files).difference(set(self.named_files)))
        return unnamed_files
    
    def gather_outputs(self):
        '''The last block has final workflow outputs'''
        for block in self.load_outputs():
            outputs = json.loads(block)
        self.raw_outputs = outputs
        self.parse_block()
            
        
    def load_outputs(self):
        read = False
        for line in self.workflow_log.split('\n'):
            if line.startswith('{'):
                block = line
                read = True
            elif line.startswith('}'):
                block += line
                read = False
                yield block
            elif read:
                block += line
                
    def ls_gsutil(self):
        try:
            result = subprocess.run(['gsutil', 'ls', self.out_dir + '/*'],
                                    check=True,
                                    stdout=subprocess.PIPE).stdout.decode('utf-8')
        except subprocess.CalledProcessError as err:
            log.error(err.output.decode('utf-8'))
            log.error('Failed to read URI: ' + uri)
            return False
        return [file for file in result.split('\n') if not file == '' 
                and not file.endswith(':')
                and not file.endswith('/')]
                    
    def read_gsutil(self):
        try:
            result = subprocess.run(['gsutil', 'cat', self.uri],
                                    check=True,
                                    stdout=subprocess.PIPE).stdout.decode('utf-8')
        except subprocess.CalledProcessError as err:
            log.error(err.output.decode('utf-8'))
            log.error('Failed to read URI: ' + uri)
            return False
        return result
        
        
        
class RunData():
    def __init__(self, run_info_file, options_file):
        self.run_info = self.read(run_info_file)
        self.options = self.read(options_file)
        self.project_info = self.run_info['project_data']
        self.passed = meta.test_schema(self.project_info)
        self.workflow_uuid = self.run_info['workflow_uuid']
        
    def read(self, file):
        with open(file) as project_info_file:
            project_info = json.load(project_info_file)
            return project_info
        

def main():
    args = get_args()
    run_data = RunData(run_info_file=args['run_data'], 
                       options_file=args['options'])
    outputs = CloudOutput(run_data=run_data)
    file_out = outputs.run_data.run_info["project_data"]['project'].replace(' ', '_') + outputs.run_data.run_info['workflow_uuid'] + '_outputInfo.json'
    with open(file_out, 'w') as project_info_file:
        json.dump(outputs.run_data.run_info, project_info_file, indent=4)

        

            
def get_args():
    '''Parse input flags
        Need to add optional task-specific input json
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--options',
                        help='Options json file.',
                            required=True
                            )
    parser.add_argument('--run-data',
                        help='JSON file with RunIfon. Includes workflow uuid, '
                        'project pairing, sample, '
                        'genome build, library and interval list information',
                        required=True
                        )
    args_namespace = parser.parse_args()
    return args_namespace.__dict__

if __name__ == "__main__":
    main()       
