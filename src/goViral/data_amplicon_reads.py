from torch.utils.data import Dataset
import json
import os
import pandas as pd
import torch

class AmpliconReads(Dataset):
    def __init__(self, input_path: str, transform=None, start=None, test_mode = False):
        self.input_path = input_path
        self.transform = transform
        self.start = start
        self.test_mode = test_mode
        
        self.reference_set = pd.read_csv(os.path.join(self.input_path), sep='\t', header=0, on_bad_lines='skip')
        
        if self.test_mode: 
            if "start" in self.reference_set.columns:
                self.reference_set = self.reference_set[self.reference_set["start"] == int(self.start)]
            else: 
                self.reference_set = self.reference_set[self.reference_set["genomic_region"].str.split("_").str[0].astype(int) == int(self.start)]

        self.length = len(self.reference_set)

    def __getitem__(self, index: int):

        items = self.reference_set.iloc[index]

        if "id1" in self.reference_set.columns:
            id1 = items["id1"]
            id2 = items["id2"]
        else:
            id1 = items["read_1"]
            id2 = items["read_2"]


        reads_dir = "/".join(self.input_path.split("/")[:-1])
        fasta_1_path = os.path.join(reads_dir, "reads",  '{}.fasta'.format(id1))
        fasta_2_path = os.path.join(reads_dir, "reads", '{}.fasta'.format(id2))
        
        with open(fasta_1_path, 'r') as fasta_1_file:
            read1 = fasta_1_file.readlines()[1].strip()
        
        
        with open(fasta_2_path, 'r') as fasta_2_file:
            read2 = fasta_2_file.readlines()[1].strip()
        

        if not self.test_mode:            
            label = items["label"]
            target = 1 if label == 'positive' or label == 1 or label == "1" else 0

        data = (read1, read2)

        if self.transform:
            data = self.transform(data)

        if self.test_mode: 
            if "start" in self.reference_set.columns:
                gr = str(items["start"]) + "_" + str(items["end"])
            else: 
                gr = items["genomic_region"]
            
            return (id1, id2), data, gr
        else:
            return data, target
    
    def __len__(self):
        return self.length