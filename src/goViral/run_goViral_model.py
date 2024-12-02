import argparse
import sys
import os
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, random_split
import torch.optim as optim
from pytorch_lightning.loggers import WandbLogger
import pytorch_lightning as pl
from pytorch_lightning.callbacks import ModelCheckpoint, StochasticWeightAveraging
from pytorch_lightning.callbacks.early_stopping import EarlyStopping
from transformers import TrainingArguments, AutoModelForSequenceClassification, AutoModelForMaskedLM
from peft import PeftConfig, PeftModel
from peft import LoraConfig, get_peft_model, TaskType, IA3Config

from goViral.data_amplicon_reads import AmpliconReads
from goViral.trainer import TransformerBinaryNetTrainer
from utils.utils import get_genomic_regions

def main(): 
    parser = argparse.ArgumentParser(description="Output predictions for pairs of reads")
    parser.add_argument('--outdir', dest = 'outdir', required=True, type=str, help="output directory")
    parser.add_argument('--path_to_dataset', dest = 'path_to_dataset', required=True, type=str, help="Path to input dataset")
    parser.add_argument('--gr_start', dest = 'gr_start', required=True, type=int, help="Genomic region start, for which the predictions are to be made")
    parser.add_argument('--accelerator', dest = 'accelerator', required=False, default = "gpu", type=str, help="gpu or cpu")
    args = parser.parse_args()
    
    out_file_path = os.path.join(args.outdir, "predictions.tsv")
  
    if os.path.exists(out_file_path):
        os.remove(out_file_path)
    
    f = open(out_file_path, "w")
    f.write("Genomic_region\tSequence_1_id\tSequence_1\tSequence_2_id\tSequence_2\tPredicted_label\tPredicted_probability\n")
    f.close()
    
    nt = AutoModelForSequenceClassification.from_pretrained("InstaDeepAI/nucleotide-transformer-v2-500m-multi-species", num_labels =2,  trust_remote_code=True, cache_dir = "/tudelft.net/staff-umbrella/ViralQuasispecies/inika/VQA/src/vqa/cache")
    adapter_name = "goViral/finetuned_model"

    model = PeftModel.from_pretrained(nt, adapter_name)
    model = model.merge_and_unload()

    model = TransformerBinaryNetTrainer(model = model, train_datal = None, val_datal = None, test_datal = None,optimizer = None, batch_size = 20, checkpoint_dir=None,  device="cuda:0", outdir = out_file_path)

    model.eval()

    # more than 1 GPUs not supported
    trainer = pl.Trainer(devices=1, accelerator= args.accelerator, enable_progress_bar=False)

    data = AmpliconReads(input_path = args.path_to_dataset, start = args.gr_start, test_mode = True)
    
    if data.length == 0: 
        print("No data found for the given genomic region")
        return
    
    datal = DataLoader(data, batch_size = 20, shuffle=False, pin_memory=True, num_workers=4, prefetch_factor=8)

    out =  trainer.predict(model, dataloaders = datal)
    

if __name__ == "__main__":
    sys.exit(main())