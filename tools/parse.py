import os
import sys
import pprint
import json
from collections.abc import Iterable
import subprocess
from google.cloud import storage
import logging as log
import google
import re 
import bicseq_config_prep
import make_auth

class Json_leaves():
    
    def __init__(self, inputs):
        self.files = []
        self.json_obj = inputs
        self.hlpr_fnc(self.json_obj)
        
    def hlpr_fnc_list(self, list1):
        for subnode in list1:
            if isinstance(subnode, dict):
                self.hlpr_fnc(subnode)
            elif not isinstance(subnode, str) and isinstance(subnode, Iterable):
                self.hlpr_fnc_list(subnode)
            else:
                if isinstance(subnode, str):
                    self.files.append(subnode)
                
    def hlpr_fnc(self, dict1):
        nodes = dict1.keys()
        for node in nodes:
            subnode = dict1[node]
            if isinstance(subnode, dict):
                self.hlpr_fnc(subnode)
            elif not isinstance(subnode, str) and isinstance(subnode, Iterable):
                self.hlpr_fnc_list(subnode)
            else:
                if isinstance(subnode, str):
                    self.files.append(subnode)
                

class Wdl():
    '''draft input JSON from various sources
    Need to add optional task-specific input json
    '''
    def __init__(self, wdl_file,
                 genome_input,
                 interval_input,
                 pipeline_input,
                 genome,
                 custom_inputs,
                 bicseq2_draft=False,
                 validate=True,
                 project_info_file=False):
        self.nygc_prefix = 'gs://nygc'
        self.nygc_public = 'gs://nygc-resources-public/'
        self.bicseq2_draft = bicseq2_draft
        self.genome = genome
        self.custom_inputs = custom_inputs
        self.validate = validate
        self.inputs = {}
        self.input_objects = {}
        if project_info_file:
            self.project_info = self.load_json(project_info_file)
        else:
            self.project_info = {}
        # find required input
        for type, variable, full_variable, default in self.load(wdl_file):
            self.inputs[variable] = {'variable' : variable,
                                     'default' : default,
                                     'full_variable' : full_variable,
                                     'object' : None,
                                     'type' : type}
        if self.bicseq2_draft:
            # load bicseq variables for custom config files
            if 'bicseq2ConfigFile' in self.inputs:
                assert 'readLength' in self.custom_inputs, 'readLength must be defined in the custom inputs json to use bicseq2_draft'
                read_length = self.inputs['readLength']['object']
                assert self.genome in self.load_json(genome_input), 'error genome not in interval file'
                genome_data = self.load_json(genome_input)[self.genome]
                upload_bucket = self.project_info['options']['final_workflow_log_dir'].replace('cromwell-logs', 'input')
                bicseq = bicseq_config_prep.Bicseq2Prep(list_of_chroms_full=genome_data['listOfChromsFull'],
                                                        uniq_coords=genome_data['uniqCoords'],
                                                        read_length=read_length,
                                                        upload_bucket=upload_bucket)
                self.load_custom_dict(custom_dict=bicseq.inputs)   
        # load preexisting reference variables
        self.load_genome_input(genome_input)
        # load preexisting pipeline reference variables
        self.load_pipeline_input(pipeline_input)
        # load preexisting interval list reference variables
        self.load_interval_input(interval_input)
        # load custom variables and overwrite any defaults from the reference jsons
        self.load_custom()
        # populate
        self.populate_inputs()
        self.finish_inputs()
        self.check_external()
        if self.validate:
            self.validate_inputs()
            
    def check_external_var(self, variable):
        if variable in self.inputs.keys():
            try:
                value = self.inputs[variable]['object']
            except TypeError:
                value = self.inputs[variable]
            if value == True:
                return True
        default = self.get_default(current_variable=variable)
        if default:
            if default == 'true':
                return True
            
    def check_external(self):
        self.external = False
        value = self.check_external_var(variable='production')
        if value:
            self.external = True
        value = self.check_external_var(variable='external')
        if value:
            self.external = True
            
    def skip(self, file):
        '''If "production" or "external" inputs are true this triggers the skipping of internal-only files 
        (the steps that run with these also are written to skip tasks that localize these files 
        so any inability to read these files will not negatively affect the run'''
        if file.startswith(self.nygc_prefix):
            if not file.startswith(self.nygc_public):
                if self.external:
                    return True
            return False
        return False
        
    def validate_inputs(self):
        potential_files = Json_leaves(self.inputs)
        files = [string for string in potential_files.files if string.startswith('gs://')]
        files = [file for file in files if not self.skip(file)]
        print(files)
        found = self.validate_input_gsutil(strings=files)
        if not found:
            log.error('searching for first missing/unreadible file. This may be slow...')
            self.narrow_down(strings=files)
        
    def parse_url(self, url):
        '''divide gcp bucket location into parts'''
        bucket_id = url.split('gs://')[-1].split('/')[0]
        project_id = '-'.join(bucket_id.split('-')[0:-1])
        name = '/'.join(url.split('gs://')[-1].split('/')[1:])
        return bucket_id, project_id, name
    
    def validate_input_gsutil(self, strings):
        try:
            result = subprocess.run(['gsutil', 'ls'] + strings,
                                    check=True,
                                    stdout=subprocess.PIPE).stdout.decode('utf-8')
        except subprocess.CalledProcessError as err:
            log.error(err.output.decode('utf-8'))
            log.error('Failed to locate file in bucket')
            return False
        return True
    
    def narrow_down(self, strings):
        for string in strings:
             if not self.validate_input_gsutil([string]):
                 print(string)
                 sys.exit(1)
        
    def validate_input(self, strings):
        '''validate that file exists in bucket'''
        self.credentials, self.gcp_project = make_auth.login()
        if string.startswith('gs://'):
            bucket_id, project_id, name = self.parse_url(string)
            storage_client = storage.Client(credentials=self.credentials,
                                            project=self.gcp_project)
            try:
                bucket = storage_client.get_bucket(bucket_id)
                # no preceeding slash!!!
                exists = storage.Blob(bucket=bucket, name=name).exists(storage_client)
                if not exists:
                    log.error('file not found in bucket: ' + string)
                    sys.exit(0)
            except google.api_core.exceptions.Forbidden:
                '''add listing of sweng files'''
                pass
        return True
            
    def check_gsutil(self, files):
        self.credentials, self.gcp_project = make_auth.login()
        if string.startswith('gs://'):
            bucket_id, project_id, name = self.parse_url(string)
            storage_client = storage.Client(credentials=self.credentials,
                                            project=self.gcp_project)
            try:
                bucket = storage_client.get_bucket(bucket_id)
                # no preceeding slash!!!
                exists = storage.Blob(bucket=bucket, name=name).exists(storage_client)
                if not exists:
                    log.error('file not found in bucket: ' + string)
                    sys.exit(0)
            except google.api_core.exceptions.Forbidden:
                '''add listing of sweng files'''
                pass
        return True
        
    def finish_inputs(self):
        final_inputs = {}
        for variable in self.inputs:
            if not self.inputs[variable]['default']:
                full_variable = self.inputs[variable]['full_variable']
                try:
                    final_inputs[full_variable] = self.inputs[variable]['object']
                except TypeError:
                    final_inputs[full_variable] = self.inputs[variable]
        self.final_inputs = final_inputs
        
    def add_from_ref_custom(self):
        '''Load from pre-existing reference JSONs and
           custom input JSONs
        '''
        for variable in self.inputs:
            if variable in self.input_objects:
                self.inputs[variable]['object'] = self.input_objects[variable]['object']
                self.inputs[variable]['default'] = False
                    
    def add_from_project(self):
        '''Load from user define project-specific information. 
        See "tools/project_schema.json" for examples of values stored in tools/project_schema.json
        (e.g. pairIds, pairInfos, genome, library, sample lists etc)
        '''
        for variable in self.inputs:
            if not self.inputs[variable]['object']:
                if variable in self.project_info:
                    self.inputs[variable]['object'] = self.project_info[variable]
                    self.inputs[variable]['default'] = False
                    
    def populate_inputs(self):
        self.add_from_ref_custom()
        self.add_from_project()
        
    def parse_input(self, line):
        line = line.replace(' , ', ',')
        line = line.replace(', ', ',')
        line = line.replace(' ,', ',')
        components = line.split()
        if len(components) == 0:
            return True, False, False
        if line.split()[0].startswith('#'):
            return True, False, False
        if components[0] == '}':
            return False, False, False
        if len(components) >= 2:
            return [True] + components[0:2]
        
    def load_json(self, file):
        with open(file) as input:
            inputs = json.load(input)
            return inputs
        
    def load_pipeline_input(self, file):
        '''load files specific to the pipeline'''
        data = self.load_json(file)
        assert len(set(data.keys()).intersection(self.input_objects.keys())) == 0, 'reference pipeline inputs must be not be redundant to other input variables'
        for variable in data:
            if variable in self.inputs:
                self.input_objects[variable] = data[variable]
        
    def load_interval_input(self, file):
        '''load variables according to the genome FASTA key and the interval list key'''
        assert self.genome in self.load_json(file), 'error genome not in interval file'
        genome_data = self.load_json(file)[self.genome]
        if self.project_info['library'] == 'WGS' or self.project_info['intervalList'] == 'default':
            interval_list = genome_data['default']
        else:
            interval_list = self.project_info['intervalList']
        data = genome_data[interval_list]
        assert len(set(data.keys()).intersection(self.input_objects.keys())) == 0, 'interval list inputs must be not be redundant to other input variables'
        for variable in data:
            if variable in self.inputs:
                self.input_objects[variable] = data[variable]
                  
    def load_genome_input(self, file):
        '''load variables according to the genome FASTA key'''
        assert self.genome in self.load_json(file), 'error genome not in genome file'
        data = self.load_json(file)[self.genome]
        assert len(set(data.keys()).intersection(self.input_objects.keys())) == 0, 'reference genome inputs must be not be redundant to other input variables'
        for variable in data:
            if variable in self.inputs:
                self.input_objects[variable] = data[variable]
                
    def load_custom_dict(self, custom_dict):
        '''load variables from any custom list of dictionaries'''
        for variable in custom_dict:
            if variable in self.input_objects.keys():
                log.warning(variable + 'value is being taken from custom generated dictionary.')
            self.input_objects[variable] = {}
            self.input_objects[variable]['object'] = custom_dict[variable]
                
    def load_custom(self):
        '''load variables from any custom file'''
        if self.custom_inputs:
            for custom_input in self.custom_inputs:
                data = self.load_json(custom_input)
                new_vars = [variable.split('.')[-1] for variable in set(data.keys())]
                assert len(set(new_vars)) == len(new_vars), 'input variables must have unique variable names (after removing workflow names).'
                for variable in data:
                    new_var = variable.split('.')[-1]
                    if new_var in self.input_objects.keys():
                        log.warning(new_var + 'value is being taken from custom inputs json file.')
                    self.input_objects[new_var] = {}
                    self.input_objects[new_var]['object'] = data[variable]
    
    def load_wf_name(self):
        return [self.inputs[variable]['full_variable'].split('.')[0] for variable in self.inputs][0]
                
    def load(self, file):
        '''Use womtool to list input files
        Replace with python parser later'''
        self.custom_types = {variable: type for type, variable in self.load_custom_type(file)}
        # inputs with no defaults
        try:
            result = subprocess.run(['womtool', 'inputs', 
                                     '--optional-inputs', 'false',
                                     file],
                                    check=True,
                                    stdout=subprocess.PIPE).stdout.decode('utf-8')
        except subprocess.CalledProcessError:
            log.error('Failed to read file: ' + file)
            sys.exit(1)
        required_inputs = [full_variable for full_variable in json.loads(result)]
        # all inputs
        try:
            result = subprocess.run(['womtool', 'inputs', 
                                     file],
                                    check=True,
                                    stdout=subprocess.PIPE).stdout.decode('utf-8')
        except subprocess.CalledProcessError:
            log.error('Failed to read file: ' + file)
            sys.exit(1)
        self.all_inputs_with_values = json.loads(result)
        self.all_inputs = [full_variable for full_variable in json.loads(result)]
        self.optional_inputs = list(set(self.all_inputs).difference(set(required_inputs)))
        for input in self.all_inputs:
            variable = input.split('.')[-1]
            # skip variables only defined with default
            if variable in self.custom_types:
                type = self.custom_types[variable]
                full_variable = input
                default = input in self.optional_inputs
                yield type, variable, full_variable, default
                
    def get_default(self, current_variable):
        ''' Lookup default values'''
        for input in self.all_inputs:
            variable = input.split('.')[-1]
            if current_variable == variable:
                # skip variables only defined with default
                self.simple_inputs_with_values = {key.split('.')[-1] : self.all_inputs_with_values[key] for key in self.all_inputs_with_values}
                if variable in self.custom_types:
                    value = self.simple_inputs_with_values[current_variable]
                    value = re.sub(r"\s+", '', value)
                    if 'default=' in value:
                        value = value.split('default=')[-1].replace(')', '').replace('\\', '').replace('"', '')
                        return value
        return False

    def load_custom_type(self, file):
        ''' get name of custom struct from WDL file'''
        read = False
        within_workflow = False
        with open(file) as input:
            for line in input:
                line = line.rstrip()
                if read:
                    read, type, variable = self.parse_input(line)
                    if read and type and variable:
                        yield type, variable
                if re.sub(r"\s+", '', line).startswith('workflow') and re.sub(r"\s+", '', line).endswith('{'):
                    within_workflow = True
                elif within_workflow and re.sub(r"\s+", '', line)  == 'input{':
                    read = True
                elif read and re.sub(r"\s+", '', line)  == '}':
                    read = False
                    within_workflow = False


def main():
    file = sys.argv[1]
    genome_input = sys.argv[2]
    interval_input = sys.argv[3]
    genome = 'Human_GRCh38_full_analysis_set_plus_decoy_hla'
    bicseq2_draft=False
    Wdl(file, genome_input=genome_input, 
        interval_input=interval_input, 
        genome=genome,
        bicseq2_draft=bicseq2_draft)


if __name__ == "__main__":
    main()       
                    