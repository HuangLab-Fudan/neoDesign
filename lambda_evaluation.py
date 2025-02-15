#from transformers import T5Tokenizer, T5EncoderModel
from transformers import TFT5EncoderModel, T5Tokenizer
import torch
import re
import gc
import numpy as np
import pickle
import ast
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

# Load the tokenizer
#tokenizer = T5Tokenizer.from_pretrained('prot_t5_xl_uniref50', do_lower_case=False,legacy=False)
tokenizer = T5Tokenizer.from_pretrained("/home/bobyu/G/neodesign/code/prot_t5_xl_uniref50", do_lower_case=False,legacy=False)
# Load the model
#model = T5EncoderModel.from_pretrained("prot_t5_xl_uniref50").to(device)
model = TFT5EncoderModel.from_pretrained("/home/bobyu/G/neodesign/code/prot_t5_xl_uniref50", from_pt=True)
# only GPUs support half-precision currently; if you want to run on CPU use full-precision (not recommended, much slower)
#if device==torch.device("cpu"):
#    model.to(torch.float32)

# prepare your protein sequences as a list
sequence_examples=[]
with open("/home/bobyu/G/neodesign/code/target_protein_sequence.txt","r") as f:
    for line in f.readlines():
        sequence_examples.append(line.strip())

#sequence_examples = ["P*RTEINO", "*SEQWENCE"]
# replace all rare/ambiguous amino acids by X and introduce white-space between all amino acids
sequence_examples = [" ".join(list(re.sub(r"[UZOB]", "X", sequence))) for sequence in sequence_examples]
seq_len = len(sequence_examples)
# tokenize sequences and pad up to the longest sequence in the batch
ids = tokenizer.batch_encode_plus(sequence_examples, add_special_tokens=True, padding=True, return_tensors="tf")
#print(ids)
input_ids = ids["input_ids"]
attention_mask = ids['attention_mask']

with torch.no_grad():
    embedding_repr = model(input_ids=input_ids, attention_mask=attention_mask)
emb_0 = embedding_repr.last_hidden_state[0,:]
emb = []
emb_per_protein = []
for i in range(0,len(sequence_examples)):
    emb.append(embedding_repr.last_hidden_state[i,:])
for i in emb:
    emb_per_protein.append(np.array(i).mean(0))
#sava embedding
features_array = np.array(emb_per_protein)
np.save('target_protTrans_features.npy', features_array)

emb_per_protein = np.load('target_protTrans_features.npy')
print(np.shape(emb_per_protein))
#model predict

import torch # Pytorch Library
from torch.nn import Conv1d # 1D Convolution Layer
from torch.nn import MaxPool1d # Max pooling Layer
from torch.nn import Flatten #Flatten Layer
from torch.nn import Linear # As dense FC NN
from torch.nn.functional import relu
from torch.utils.data import DataLoader, TensorDataset # DataLoader: Bach size, Shuffle
from torch.nn.functional import tanh, dropout
from torch.nn import AvgPool1d
from torch.nn import BatchNorm1d
import torch.nn as LeakyRelu
import torch
from torch.nn import Conv1d, MaxPool1d, Flatten, Linear, BatchNorm1d
import torch.nn.functional as F
from torch.utils.data import DataLoader, TensorDataset
from torch import nn

#@title Default title text
class CnnRegressor(torch.nn.Module):
    def __init__(self, batch_size, inputs, outputs):
        super(CnnRegressor, self).__init__()
        self.batch_size = batch_size
        self.inputs = inputs
        self.outputs = outputs

        self.input_layer = Conv1d(inputs, batch_size, 1, padding=0)
        #self.max_pooling_layer = MaxPool1d(1)
        self.conv_layer = Conv1d(batch_size, 128, 1)
        self.bn = BatchNorm1d(num_features=128)
        #self.avg_pooling_layer = AvgPool1d(1)
        self.conv_layer2 = Conv1d(128,128, 1)
        self.conv_layer3 = Conv1d(128,128, 1)
        self.conv_layer4 = Conv1d(128,128, 1)
        self.flatten_layer = Flatten()
        self.liner_layer = Linear(128,128,1)
        self.liner_layer2 = Linear(128,128,1)
        self.liner_layer3 = Linear(128,128)
        self.liner_layer4 = Linear(128,128)
        self.output_layer = Linear(128, outputs)

    def feed(self, input):
        input = input.reshape((self.batch_size, self.inputs, 1))
        output = relu(self.bn(self.input_layer(input)))
        output = relu(self.conv_layer(output))
        output = relu(self.conv_layer2(output))
        output = relu(self.conv_layer3(output))
        output = relu(self.conv_layer4(output))
        output = self.flatten_layer(output)
        output = self.liner_layer(output)
        output = self.liner_layer2(output)
        output = self.liner_layer3(output)
        output = self.liner_layer4(output)
        output = self.output_layer(output)
        return output

batch_size = 128

model = CnnRegressor(batch_size, 1024,1) # Parameter 1?

if torch.cuda.is_available():
    model = model.cuda()
else:
    model = model


#min model load
device = torch.device("cpu")
model = torch.load('/home/bobyu/G/neodesign/code/model_min.pth', map_location=device)

emb_num = len(emb_per_protein)
emb_per_protein = emb_per_protein.tolist()
merged_list = emb_per_protein*128
merged_list = merged_list[0:128]
print(len(merged_list))
if torch.cuda.is_available():
    inputs = torch.from_numpy(np.array(merged_list)).cuda().float()
else:
    inputs = torch.from_numpy(np.array(merged_list)).float()

model.eval()  
with torch.no_grad():  
    output = model.feed(inputs)  
tensor_list_min = output.tolist()
tensor_list_min = tensor_list_min[0:emb_num]


#max model load
model = torch.load('/home/bobyu/G/neodesign/code/model_max.pth',map_location=device)

if torch.cuda.is_available():
    inputs = torch.from_numpy(np.array(merged_list)).cuda().float()
else:
    inputs = torch.from_numpy(np.array(merged_list)).float()

model.eval()  
with torch.no_grad(): 
    output = model.feed(inputs)  
tensor_list_max = output.tolist()
tensor_list_max = tensor_list_max[0:emb_num]

mid = []
print(tensor_list_min)
print(tensor_list_max)
max_l = max(tensor_list_max)
min_l = min(tensor_list_min)
diff = max_l[0] - min_l[0]
for i in range(0,emb_num):
    if abs(tensor_list_min[i][0] + tensor_list_max[i][0])/2*0.0589 <20:
        mid.append(abs(tensor_list_min[i][0] + tensor_list_max[i][0])/2*0.0589)
    else:
        mid.append(20)

with open("/home/bobyu/G/neodesign/code/output_lambda.txt","w") as f:
    f.write("min")
    f.write("\t")
    f.write("max")
    f.write("\t")
    f.write("recommended lambda")
    f.write("\n")
    for i in range(0,emb_num):
        f.write(str(tensor_list_min[i][0]))
        f.write("\t")
        f.write(str(tensor_list_max[i][0]))
        f.write("\t")
        f.write(str(mid[i]))
        f.write("\n")
