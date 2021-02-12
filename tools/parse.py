import os
import sys
import pprint
import json

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
        self.wf_name = self.load_wf_name(wdl_file)
        if project_info_file:
            self.project_info = self.load_json(project_info_file)
        else:
            self.project_info = {}
        # find required input
        for type, variable in self.load(wdl_file):
            self.inputs[variable] = {'variable' : variable,
                                     'object' : None,
                                     'type' : type}
        # load preexisting reference variables
        self.load_genome_input(genome_input)
        # populate
        self.populate_inputs()
        self.finish_inputs()
        
    def finish_inputs(self):
        final_inputs = {}
        for variable in self.inputs:
            try:
                final_inputs[self.wf_name + '.' + variable] = self.inputs[variable]['object']
            except TypeError:
                final_inputs[self.wf_name + '.' + variable] = self.inputs[variable]
        self.inputs = final_inputs
        
    def add_from_ref(self):
        '''Load from pre-existing reference data'''
        for variable in self.inputs:
            # skip resource vars because these will be hardcoded
            if not variable.endswith('DockerImage') and not variable.endswith('threads') and not variable.endswith('Mem') and not variable.endswith('mem'):
                if variable in self.input_objects:
                    self.inputs[variable] = self.input_objects[variable]
                    
    def add_from_project(self):
        '''Load from user define project-specific information 
        (e.g. pairing location of input files etc)
        '''
        for variable in self.inputs:
            # skip resource vars because these will be hardcoded
            if not variable.endswith('DockerImage') and not variable.endswith('threads') and not variable.endswith('Mem') and not variable.endswith('mem'):
                if not self.inputs[variable]['object']:
                    if variable in self.project_info:
                        self.inputs[variable] = self.project_info[variable]
                    
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
        if len(components) == 2:
            return [True] + components
        if len(components) > 2:
            # skip lines where input variables are 
            # composed from other input variables with an 
            # equal sign
            return True, False, False
        
    def load_json(self, file):
        with open(file) as input:
            inputs = json.load(input)
            return inputs
        
    def load_genome_input(self, file):
        assert self.genome in self.load_json(file), 'error genome not in genome file'
        data = self.load_json(file)[self.genome]
        assert len(set(data.keys()).intersection(self.input_objects.keys())) == 0, 'inputs must be unique'
        for variable in data:
            self.input_objects[variable] = data[variable]
            
    
    def load_wf_name(self, file):
        with open(file) as input:
            for line in input:
                if line.startswith('workflow '):
                    components = line.split()
                    name = components[1]
                    return name
                      
    def load(self, file):
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
                    