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
        self.reference_set = self.reference_set[self.reference_set["start"] == int(self.start)]
        self.length = len(self.reference_set)

    def __getitem__(self, index: int):

        items = self.reference_set.iloc[index]

        id1 = items["id1"]
        id2 = items["id2"]

        reads_dir = "/".join(self.input_path.split("/")[:-1])
        fasta_1_path = os.path.join(reads_dir, "reads",  '{}.fasta'.format(id1))
        fasta_2_path = os.path.join(reads_dir, "reads", '{}.fasta'.format(id2))
        
        with open(fasta_1_path, 'r') as fasta_1_file:
            read1 = fasta_1_file.readlines()[1].strip()
        
        
        with open(fasta_2_path, 'r') as fasta_2_file:
            read2 = fasta_2_file.readlines()[1].strip()
        
        gr = str(items["start"]) + "_" + str(items["end"])
        
        if not self.test_mode:
            label = items["label"]
            target = 1 if label == 'positive' else 0

        data = (read1, read2)

        if self.transform:
            data = self.transform(data)

        if self.test_mode: 
            return (id1, id2), data, gr
        else:
            return data, target
    
    def __len__(self):
        return self.length