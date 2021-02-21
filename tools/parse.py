import os
import sys
import pprint
import json
import subprocess

class Wdl():
    '''draft input JSON from various sources
    Need to add optional task-specific input json
    '''
    def __init__(self, wdl_file,
                 genome_input,
                 interval_input,
                 genome,
                 project_info_file=False):
        self.genome = genome
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
        self.wf_name = self.load_wf_name()
        # load preexisting reference variables
        self.load_genome_input(genome_input)
        # populate
        self.populate_inputs()
        self.finish_inputs()
        
    def finish_inputs(self):
        final_inputs = {}
        for variable in self.inputs:
            if not self.inputs[variable]['default']:
                full_variable = self.inputs[variable]['full_variable']
                try:
                    final_inputs[full_variable] = self.inputs[variable]['object']
                except TypeError:
                    final_inputs[full_variable] = self.inputs[variable]
        self.inputs = final_inputs
        
    def add_from_ref(self):
        '''Load from pre-existing reference data'''
        for variable in self.inputs:
            # skip resource vars because these will be hardcoded
            if not variable.endswith('DockerImage') and not variable.endswith('threads') and not variable.endswith('Mem') and not variable.endswith('mem'):
                if variable in self.input_objects:
                    self.inputs[variable]['object'] = self.input_objects[variable]['object']
                    self.inputs[variable]['default'] = False
                    
    def add_from_project(self):
        '''Load from user define project-specific information 
        (e.g. pairing location of input files etc)
        '''
        for variable in self.inputs:
            # skip resource vars because these will be hardcoded
            if not variable.endswith('DockerImage') and not variable.endswith('threads') and not variable.endswith('Mem') and not variable.endswith('mem'):
                if not self.inputs[variable]['object']:
                    if variable in self.project_info:
                        self.inputs[variable]['object'] = self.project_info[variable]
                        self.inputs[variable]['default'] = False
                    
    def populate_inputs(self):
        self.add_from_ref()
        self.add_from_project()
        
    def parse_input(self, line):
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
        
    def load_genome_input(self, file):
        assert self.genome in self.load_json(file), 'error genome not in genome file'
        data = self.load_json(file)[self.genome]
        assert len(set(data.keys()).intersection(self.input_objects.keys())) == 0, 'inputs must be unique'
        for variable in data:
            if variable in self.inputs:
                self.input_objects[variable] = data[variable]
    
    def load_wf_name(self):
        return [self.inputs[variable]['full_variable'].split('.')[0] for variable in self.inputs][0]
                
    def load(self, file):
        '''Use womtool to list input files
        Replace with python parser later'''
        custom_types = {variable: type for type, variable in self.load_custom_type(file)}
        try:
            result = subprocess.run(['womtool', 'inputs', file],
                                    check=True,
                                    stdout=subprocess.PIPE).stdout.decode('utf-8')
        except CalledProcessError:
            log.error('Failed to read file: ' + file)
        inputs = json.loads(result)
        for input in inputs:
            variable = input.split('.')[-1]
            # skip variables only defined with default
            if variable in custom_types:
                type = custom_types[variable]
                full_variable = input
                default = 'default =' in inputs[input]
                yield type, variable, full_variable, default
        
                      
    def load_custom_type(self, file):
        ''' get name of custom struct from WDL file'''
        read = False
        with open(file) as input:
            for line in input:
                line = line.rstrip()
                if read:
                    read, type, variable = self.parse_input(line)
                    if read and type and variable:
                        yield type, variable 
                if line.replace(' ', '') == 'input{':
                    read = True


def main():
    file = sys.argv[1]
    genome_input = sys.argv[2]
    interval_input = sys.argv[3]
    genome = 'Human_GRCh38_full_analysis_set_plus_decoy_hla'
    Wdl(file, genome_input=genome_input, interval_input=interval_input, genome=genome)


if __name__ == "__main__":
    main()       
                    