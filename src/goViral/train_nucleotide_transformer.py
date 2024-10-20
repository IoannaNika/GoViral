import os
import wandb
import argparse
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, random_split
import torch.optim as optim
from pytorch_lightning.loggers import WandbLogger
import pytorch_lightning as pl
from pytorch_lightning.callbacks import ModelCheckpoint, StochasticWeightAveraging
from pytorch_lightning.callbacks.early_stopping import EarlyStopping
from transformers import TrainingArguments, AutoModelForSequenceClassification, AutoModelForMaskedLM, AutoModel
from peft import LoraConfig, get_peft_model, TaskType, IA3Config

from goViral.data_amplicon_reads import AmpliconReads
from goViral.trainer import TransformerBinaryNetTrainer

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--devices', type=int, default=2, required=False)
    parser.add_argument('--checkpoint_dir', type=str, required=True)
    parser.add_argument('--pb', action="store_true", help="if False will assume ONT reads")
    args = parser.parse_args()

    device = "cuda:0"
    n_devices = args.devices

    check_point_dir = args.checkpoint_dir
    os.system(f"mkdir -p {check_point_dir}")

    max_length = 800 
    batch = 10

    model = AutoModelForSequenceClassification.from_pretrained("InstaDeepAI/nucleotide-transformer-v2-500m-multi-species", trust_remote_code=True, num_labels=2, cache_dir = "/tudelft.net/staff-umbrella/ViralQuasispecies/inika/VQA/src/vqa/cache")

    peft_config = IA3Config(peft_type="IA3", target_modules=["value", "key", "intermediate.dense"], feedforward_modules=["intermediate.dense"], modules_to_save=["classifier"])
    
    peft_model = get_peft_model(model, peft_config)

    # make sure the classification head will be fine-tuned
    for name, param in peft_model.named_parameters():
        if "classifier" in name:
            param.requires_grad = True
    
    if args.pb:
        print("Pac-Bio reads")
        # pacbio-hifi training set 
        data = AmpliconReads(input_path='/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/tuples_pacbio_sars_cov_2_rev_compl/dataset/samples.tsv')
    else:
        print("ONT reads")
        # ONT training set 
        data = AmpliconReads(input_path ='/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/tuples_ONT_sars_cov_2_rev_compl/dataset/samples.tsv')
    
    train_count = int(len(data)*0.9)
    val_count = len(data) - train_count
    test_count =  0 #len(data) - train_count - val_count
    
    train_data, val_data, test_data = random_split(data, [train_count, val_count, test_count])


    train_datal = DataLoader(train_data, batch_size=batch, shuffle=True, pin_memory=True, num_workers=4, prefetch_factor=1)
    val_datal = DataLoader(val_data, batch_size=batch , shuffle=False, pin_memory=True, num_workers=4, prefetch_factor=1)
    test_datal = DataLoader(test_data, batch_size=batch, shuffle=False, pin_memory=True, num_workers=4, prefetch_factor=1)
      
    optimizer = optim.Adam(peft_model.parameters(), lr=7e-3)
    
    os.environ["WANDB_DIR"] = "/tmp"
    os.environ["WANDB_START_METHOD"]="thread"
    wandb.init(project="New_nt_binary")
    wandb_logger = WandbLogger()
    # log loss per epoch
    wandb_logger.watch(peft_model) 
    
    early_stop_callback = EarlyStopping(monitor="val_acc", patience=3, verbose=False, mode="max")
    
    binary_transformer = TransformerBinaryNetTrainer(peft_model, train_datal, val_datal, test_datal, optimizer, batch, device = device, checkpoint_dir = check_point_dir, max_length=max_length, outdir = None)
    
    trainer = pl.Trainer(max_epochs=50, logger=wandb_logger,  accumulate_grad_batches=10, strategy='ddp_find_unused_parameters_true', callbacks=[early_stop_callback], devices=n_devices, accelerator="gpu", enable_progress_bar=False)
    trainer.fit(binary_transformer)
    
    wandb_logger.experiment.unwatch(peft_model)
    
    # trainer.test()

if __name__ == "__main__":
    main()