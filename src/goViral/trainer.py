import numpy as np
import pytorch_lightning as pl
import torch.optim as optim
import torch
import wandb
import torch.nn.functional as F
import torch.nn as nn
from transformers import AutoTokenizer

class TransformerBinaryNetTrainer(pl.LightningModule):
    def __init__(self, model, train_datal, val_datal, test_datal, optimizer, batch_size, checkpoint_dir, outdir, max_length=800, device="cpu"):
        super().__init__()
        self.model = model
        self.tokenizer = AutoTokenizer.from_pretrained("InstaDeepAI/nucleotide-transformer-v2-500m-multi-species", trust_remote_code=True)
        self.batch_size = batch_size
        self.max_length = max_length
        self.optimizer = optimizer
        self.train_datal = train_datal
        self.test_datal = test_datal
        self.val_datal = val_datal
        self.criterion = nn.CrossEntropyLoss()
        self.device = device
        self.checkpoint_dir = checkpoint_dir
        self.predictions_file_name = outdir 

    def forward(self, input1, input2):
        inpt = []

        for i in range(len(input1)):
            inpt.append(input1[i].strip() + "<mask><mask><mask><mask><mask><mask><mask><mask><mask><mask><mask><mask><mask><mask>" + input2[i].strip())

        inpt = self.tokenizer.batch_encode_plus(inpt, return_tensors="pt", padding="max_length", max_length = self.max_length, truncation=True)["input_ids"].to(self.device)
        
        attention_mask = inpt != self.tokenizer.pad_token_id
        
        output = self.model(inpt,
            attention_mask=attention_mask,
            output_hidden_states=True)

        return  output.logits
   
    def training_step(self, batch, batch_idx):
        (x0, x1), y = batch
        output = self.forward(x0, x1)
        loss = self.criterion(output, y)
        accuracy = torch.sum(torch.argmax(output, dim=1) == y).item() / len(y)
        wandb.log({"epoch": self.current_epoch, "train/loss": loss})
        wandb.log({"epoch": self.current_epoch, "train/accuracy": accuracy})
        return loss

    def validation_step(self, batch, batch_idx):
        self.model.save_pretrained("{}/{}".format(self.checkpoint_dir, self.current_epoch))
        (x0, x1), y = batch
        output = self.forward(x0, x1)
        loss = self.criterion(output, y)
        accuracy = torch.sum(torch.argmax(output, dim=1) == y).item() / len(y)

        self.log("val_loss", loss,  batch_size = self.batch_size) 
        self.log("val_acc", accuracy,  batch_size = self.batch_size)

        print("validation accuracy: ", accuracy)

        wandb.log({"epoch": self.current_epoch, "val/loss": loss})
        wandb.log({"epoch": self.current_epoch, "val/accuracy": accuracy})

        return loss
    
    def test_step(self, batch, batch_idx):
        (x0, x1), y = batch

        output = self.forward(x0, x1)
        loss = self.criterion(output, y)
        accuracy = torch.sum(torch.argmax(output, dim=1) == y).item() / len(y)

        print("test accuracy: ", accuracy)

        wandb.log({"epoch": self.current_epoch, "test/loss": loss})
        wandb.log({"epoch": self.current_epoch, "test/accuracy": accuracy})

        return loss
    

    def predict_step(self, batch, batch_idx):
        (read_id_1,read_id_2), (x0, x1), gr = batch
        output = self.forward(x0, x1)
        
        for i in range(len(output)):
            file = open(self.predictions_file_name, 'a')
            file.write(gr[i] + "\t" + read_id_1[i] + "\t" + x0[i] + "\t" + read_id_2[i] + "\t" + x1[i] + "\t" + str(torch.argmax(output[i]).item()) + "\t" + str(torch.max(F.softmax(output[i]), dim=0).values.item()) + "\n") 
            file.close()
            
        return output
 
    def configure_optimizers(self):
        return self.optimizer
    
    def train_dataloader(self):
        return self.train_datal
    
    def test_dataloader(self):
        return self.test_datal
    
    def val_dataloader(self):
        return self.val_datal

    def device(self): 
        return self.device
    