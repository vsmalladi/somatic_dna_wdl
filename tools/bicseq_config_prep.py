import pandas as pd
import numpy as np
import os
from google.cloud import storage


class Bicseq2Prep():
    ''''Prep the non-sample-specific and non-filesystem-specific 
    portion of the config files.
    Also returns the Bicseq2 WDL input variables that will be required for the run.'''
    def __init__(self, pair_relationships,
                 list_of_chroms_full,
                 uniq_coords,
                 read_length,
                 upload_bucket
                 ):
        self.list_of_chroms_full = list_of_chroms_full['object']
        self.upload_bucket_uri = upload_bucket
        self.pair_relationships = pair_relationships
        self.read_length = self.match_read_length(uniq_coords, read_length)
        self.uniq_coords = uniq_coords['object'][self.read_length]
        self.inputs = {}
        self.inputs['coordReadLength'] = self.read_length
        self.inputs['readLength'] = read_length
        self.inputs['bicseq2ConfigMaps'] = {}
        for pair_relationship in self.pair_relationships:
            self.inputs['bicseq2ConfigMaps'][pair_relationship['pairId']] = {}
            self.write_sample_configs(pair_relationship)
            self.write_pair_config(pair_relationship)
        
    def write_sample_configs(self, pair_relationship):
        '''create and upload configs for the normalization steps'''
        # tumor
        data = self.prep(pair_relationship['tumor'])
        file = pair_relationship['tumor'] + '.bicseq2.config'
        config = data.to_csv(sep='\t', index=False)
        uri = self.upload(file, config)
        self.inputs['bicseq2ConfigMaps'][pair_relationship['pairId']]['tumorConfigFile'] = uri
        # normal
        data = self.prep(pair_relationship['normal'])
        file = pair_relationship['normal'] + '.bicseq2.config'
        config = data.to_csv(sep='\t', index=False)
        uri = self.upload(file, config)
        self.inputs['bicseq2ConfigMaps'][pair_relationship['pairId']]['normalConfigFile'] = uri
        
    def write_pair_config(self, pair_relationship):
        file = pair_relationship['pairId'] + '.bicseq2.seg.config'
        data = self.prep_pair(pair_relationship)
        config = data.to_csv(sep='\t', index=False)
        uri = self.upload(file, config)
        self.inputs['bicseq2ConfigMaps'][pair_relationship['pairId']]['segConfigFile'] = uri
            
    def parse_url(self, url):
        '''divide gcp bucket location into parts'''
        bucket_id = url.split('gs://')[-1].split('/')[0]
        project_id = '-'.join(bucket_id.split('-')[0:-1])
        name = '/'.join(url.split('gs://')[-1].split('/')[1:])
        return bucket_id, project_id, name
            
    def upload(self, file, config):
        '''upload file content (as a string) with a given filename'''
        bucket_id, project_id, name = self.parse_url(self.upload_bucket_uri)
        client = storage.Client()
        bucket = client.get_bucket(bucket_id)
        blob = bucket.blob(file)
        blob.upload_from_string(config)
        return os.path.join(self.upload_bucket_uri, file)
        
    def match_read_length(self, uniq_coords, read_length):
        '''
        Make config file that includes the path
        to mappability file.
        '''
        length_key = False
        for length in uniq_coords['object']:
            if np.isclose(int(read_length), int(length), atol=2):
                length_key = length
                return length_key
            
    def convert_path(self, uri):
        '''convert to gcp path (not ideal but forced by config-style running'''
        return uri.replace('gs://', '/cromwell_root/')
            
    def prep_pair(self, pair_relationship):
        ''' file: tumor--normal.bicseq2.seg.config
        prep fasta-specific but sample independent portion of config file'''
        data = pd.DataFrame({'chr' : [chrom for chrom in self.list_of_chroms_full]})
        return data
        
    def prep(self, sample_id):
        ''' sample_id.bicseq2.config
        prep fasta-specific but sample independent portion of config file'''
        data = pd.DataFrame({'chrom_name' : [chrom for chrom in self.list_of_chroms_full]})
        return data
