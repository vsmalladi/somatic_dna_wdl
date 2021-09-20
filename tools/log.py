import json
import sys
import argparse
import datetime
import os

class Log():
    def __init__(self,
                 project_data, 
                 uuid, 
                 inputs):
        self.project_data_file = project_data
        self.uuid = uuid
        self.input_files = inputs
        self.project_data = self.load_json(self.project_data_file)
        self.run_info = {'workflow_uuid': self.uuid,
                         'project_data': self.project_data,
                         }
        self.inputs = {}
        for file in self.input_files:
            input = self.load_json(file)
            for key in input:
                self.inputs[key] = input[key]
        self.run_info['inputs'] = self.inputs
        self.write()
        
    def load_json(self, file):
        with open(file) as input:
            inputs = json.load(input)
            return inputs

    def write(self):
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
        self.run_info['project_data']['run_date'] = timestamp
        file_out = os.path.basename(self.project_data_file).replace('Info.json', '') + '.' +  timestamp + '.RunInfo.json'
        with open(file_out, 'w') as input_info_file:
            json.dump(self.run_info, input_info_file, indent=4)
        
def get_args():
    '''Parse input flags
        Need to add optional task-specific input json
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--project-data',
                        help='JSON file with project pairing, sample, '
                        'genome build, library and interval list information',
                        required=True
                        )
    parser.add_argument('--uuid',
                        help='Cromwell submission uuid',
                        required=True
                        )
    parser.add_argument('--inputs',
                        help='JSON file(s) with inputs for cromwell run',
                        required=True,
                        nargs='+'
                        )
    args_namespace = parser.parse_args()
    return args_namespace.__dict__
        
if __name__ == "__main__":
    args = get_args()
    Log(args['project_data'], args['uuid'], args['inputs'])    