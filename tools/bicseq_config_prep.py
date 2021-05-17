import pandas as pd
import numpy as np
import os
from google.cloud import storage


class Bicseq2Prep():
    def __init__(self, pair_relationships,
                 chrom_fastas,
                 uniq_coords,
                 read_length,
                 upload_bucket
                 ):
        
        self.upload_bucket_uri = upload_bucket
        self.pair_relationships = pair_relationships
        self.read_length = self.match_read_length(uniq_coords, read_length)
        self.uniq_coords = uniq_coords['object'][self.read_length]
        self.chrom_fastas = chrom_fastas['object']
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
            
    def prep_pair(self, pair_relationship):
        data = pd.DataFrame({'chr' : [chrom for chrom in self.chrom_fastas]})
        data['case'] = data.apply(lambda row: pair_relationship['tumor'] + '/' + 
                                  pair_relationship['tumor'] + '_' + 
                                  row.chr + '.norm.bin.txt', axis=1)
        data['control'] = data.apply(lambda row: pair_relationship['normal'] + '/' +
                                      pair_relationship['normal'] + '_' + 
                                      row.chr + '.norm.bin.txt', axis=1)
        return data
        
    def prep(self, sample_id):
        data = pd.DataFrame({'chrom_name' : [chrom for chrom in self.chrom_fastas],
                             'fa_file' : [self.chrom_fastas[chrom] for chrom in self.chrom_fastas],
                             'mappability' : [self.uniq_coords[chrom] for chrom in self.uniq_coords]})
        data['readPosFile'] = data.apply(lambda row: sample_id + '/' + sample_id + '_' + row.chrom_name + '.seq', axis=1)
        data['bin_file_normalized'] = data.apply(lambda row: sample_id + '/' + sample_id + '_' + 
                                                 row.chrom_name + '.norm.bin.txt', axis=1)
        return data
