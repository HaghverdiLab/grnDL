import argparse
import json
import os
import sys
import torch
import numpy as np
import pandas as pd

print("___________________________________________________________________________________________________________")
currentdir = os.path.dirname(os.path.abspath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
from GNNProject.classifier import Classifier
from GNNProject.dataset import Dataset
from GNNProject.hyperparameters import get_hyperparams
from GNNProject.utils import load_adj
parser = argparse.ArgumentParser()
parser.add_argument('-tag', type=str, default='Syn', help='Input dataset. Can either be the name of the synthetic dataset generation method, \
                                              or the path to a real dataset.')
parser.add_argument('-i', type=str,default='sign', help='Input dataset. Can either be the name of the synthetic dataset generation method, \
                                              or the path to a real dataset.')
parser.add_argument('--oo', type=str,default='../output', help='Output directory.')
parser.add_argument('--n_classes',type=int,default=3,help='Number of classes')
parser.add_argument('--n_features',type=int,default=100,help='Number of features (in case of synthetic data')
parser.add_argument('--n_char_features',type=int,default=10,help='Number of characteristic features for each class')
parser.add_argument('--n_obs_train',type=int,default=500,help='Number of observations for training (per class)')

parser.add_argument('--n_obs_test',type=int,default=500,help='Number of observations for testing (per class)')
parser.add_argument('--graph_model',type=str,default='BA',help='Graph generation model')
parser.add_argument('--signal_train',type=float,default=0,help='Signal for training data')
parser.add_argument('--signal_test',type=float,default=0,help='Signal for test data')
parser.add_argument('--diff_train',type=float,default=0.3,help='Diffusion coefficient for training data')
parser.add_argument('--diff_test',type=float,default=0.3,help='Diffusion coefficient for test data')
parser.add_argument('--noise_train',type=float,default=0.2,help='Noise for training data')
parser.add_argument('--noise_test',type=float,default=0.2,help='Noise for test data')
parser.add_argument('--n_iter',type=int,default=1,help='Number of diffusion iterations')
parser.add_argument('--n_communities',type=int,default=5,help='Number of communities (for SBM)')
parser.add_argument('--probs_train',type=float,default=0.9,help='SBM probabilities for training')
parser.add_argument('--probs_test',type=float,default=0.1,help='SBM probabilities for testing')

parser.add_argument('--classifier',type=str,default="MLP",help='Type of classifier')
parser.add_argument('--n_hidden_GNN',type=int,default=0,help='Number of hidden features for the GNN. If 0, do not use a GNN.')
parser.add_argument('--n_hidden_FC',type=int,default=0,help='Number of features in the fully connected hidden layer. If 0, do not use a hidden layer.')
parser.add_argument('--n_hidden_FC2',type=int,default=0,help='Number of features in the 2nd fully connected hidden layer. If 0, do not use a 2nd hidden layer.')
parser.add_argument('--n_hidden_FC3',type=int,default=0,help='Number of features in the 3rd fully connected hidden layer. If 0, do not use a 3rd hidden layer.')
parser.add_argument('--K',type=int,default=2,help='Parameter for Chebnet GNN.')
parser.add_argument('--epochs',type=int,default=30,help='Number of training epochs.')
parser.add_argument('--batch_size',type=int,default=16,help='Batch size for the neural networks.')

parser.add_argument('--infer_graph',type=str,default="False",help="Whether to infer the graph from the data (True) or directly use the true graph.")
parser.add_argument('--FPR',type=float,default=0.0,help="Percentage of false edges in the noisy graph.")
parser.add_argument('--FNR',type=float,default=0.0,help="Percentage of true edges which are removed in the noisy graph.")
parser.add_argument('--alpha',type=float,help='Parameter for graph inferrence.')
parser.add_argument('--CV_alpha',type=str,default='False',help='Whether or not to run CV for the dropout parameter.')
parser.add_argument('--dropout',type=float,default=0.2,help='Dropout rate.')
parser.add_argument('--CV_dropout',type=str,default='False',help='Whether or not to run CV for the dropout parameter.')

parser.add_argument('--seed',type=int,default=0,help='Seed for generating the synthetic dataset.')

print("___________________________________________________________________________________________________________")
args = parser.parse_args()

print("Starting...")
n_hidden_GNN = [] if args.n_hidden_GNN==0 else [args.n_hidden_GNN, 8]

if args.n_hidden_FC==0:
  n_hidden_FC=[]
else:
  if args.n_hidden_FC2==0:
    n_hidden_FC = [args.n_hidden_FC]
  else:
    if args.n_hidden_FC3==0:
      n_hidden_FC = [args.n_hidden_FC,args.n_hidden_FC2]
    else:
      n_hidden_FC = [args.n_hidden_FC,args.n_hidden_FC2,args.n_hidden_FC3]

use_true_graph = args.infer_graph=="False" and args.FPR==0.0 and args.FNR==0.0
infer_graph = args.infer_graph=="True"

if torch.cuda.is_available():
  dev = "cuda:0"
else:
  dev = "cpu"
device = torch.device(dev)

# Create a synthetic dataset or load a real dataset

if args.i in ["diffusion","activation","sign"]:
  dataset = Dataset(tag='Syn_' + args.i,random_seed=args.seed)
  dataset.create_syn(n_classes = args.n_classes,
                  n_obs_train = args.n_obs_train,
                  n_obs_test= args.n_obs_test,
                  n_features=args.n_features,
                  n_char_features = args.n_char_features,
                  signal =[args.signal_train, args.signal_test],
                  diff_coef=[args.diff_train, args.diff_test],
                  noise = [args.noise_train, args.noise_test],
                  n_communities = args.n_communities,
                  probs = [args.probs_train, args.probs_test],
                  n_iter=args.n_iter,
                  model =args.graph_model,
                  syn_method=args.i)
else:
  dataset = Dataset(tag=args.tag, input_dir=args.i)
  dataset.load()
  #dataset.subsample(n_obs_train=args.n_obs_train, n_obs_test=args.n_obs_test)
  args.n_features = dataset.X_train.shape[1]
  args.n_classes  = len(set(dataset.y_train.tolist()))
# Select the dropout rate and alpha parameter
CV_alpha = args.CV_alpha=="True"
CV_dropout = args.CV_dropout=="True"
dropout_rate,alpha = get_hyperparams(CV_dropout=CV_dropout,CV_alpha=CV_alpha,dataset=dataset,n_features=args.n_features, n_obs_train=args.n_obs_train,\
        n_classes=args.n_classes,n_hidden_GNN=n_hidden_GNN,n_hidden_FC=n_hidden_FC,K=args.K,classifier=args.classifier,\
        lr=0.001,momentum=0.9,epochs=args.epochs,device=device,batch_size=args.batch_size,use_true_graph=use_true_graph,dropout_rate=args.dropout,alpha=args.alpha)

# Create pytorch geometric dataloader
#if infer_graph: # Infer the graph with glasso
 # dataset.create_graph(alphas=alpha)
#else: # use a noisy version of the true graph
#  dataset.create_noisy_true_graph(FPR=args.FPR,FNR=args.FNR)

if args.classifier == 'MLP':
  adj = np.zeros((args.n_features, args.n_features))
else:
  adj = load_adj(args.i + args.tag)[0]

if args.classifier == 'MLP':
  adj = np.zeros((args.n_features, args.n_features))
else:
  adj = load_adj(args.i + args.tag)[0]

dataset.Ah_train = dataset.Ah_test = adj

print("Dataloader train...")
train_dataloader = dataset._dataloader('train',use_true_graph=False,batch_size=args.batch_size)

print("Data loader test...")
test_dataloader  = dataset._dataloader('test',use_true_graph=False,batch_size=args.batch_size, shuffle=False)
# Fit and evaluate a classifier
clf = Classifier(
        n_features=args.n_features,
        n_classes=args.n_classes,
        n_hidden_GNN=n_hidden_GNN,
        n_hidden_FC=n_hidden_FC,
        dropout_GNN=dropout_rate,
        dropout_FC=dropout_rate,
        K=args.K,
        classifier=args.classifier,
        lr=.001,
        momentum=.9,
        log_dir=None, # 
        device=device)

print(".......")
def save_np_txt(ndarray, path, colnames=None, rownames=None):
    df = pd.DataFrame(data=ndarray, index=rownames, columns=colnames)
    df.to_csv(path, sep='\t', index=True, header=True)

    
clf.fit(train_dataloader, epochs = args.epochs, test_dataloader=None,verbose=True)
#print("Weights?")
#print(clf.net)
#weight =clf.net.layers_FC[0].weight.data.cpu().detach().numpy()
#print(type(weight))
#print(weight.shape)
#fl_file = os.path.join(args.oo, "finallayer.txt")
#save_np_txt(weight, fl_file)

results = clf.eval(test_dataloader, verbose=True,path=args.oo)
print(results)



print("Get files...")
# Store the results in a json file.
output = {'tag': args.tag, "accuracy":results[0],"precision":results[2],"recall":results[3],"f1":results[4], "n_classes":args.n_classes,"n_features":args.n_features,\
      "n_char_features":args.n_char_features,"n_obs_train":args.n_obs_train,"n_obs_test":args.n_obs_test,\
      "noise_train":args.noise_train,"noise_test":args.noise_test,\
        "classifier":args.classifier,"n_hidden_GNN":args.n_hidden_GNN,"n_hidden_FC":args.n_hidden_FC,"n_hidden_FC2":args.n_hidden_FC2,"n_hidden_FC3":args.n_hidden_FC3,\
          "seed":args.seed,"alpha":str(args.alpha),\
          "FPR":args.FPR,"FNR":args.FNR,"infer_graph":args.infer_graph,'K': args.K}

filename = "_".join([args.tag, args.classifier,str(args.n_hidden_GNN),str(args.K),str(args.n_hidden_FC),str(args.n_hidden_FC2),str(args.n_hidden_FC3)]) + ".json"
if not os.path.exists(args.oo):
    os.makedirs(args.oo)
with open(os.path.join(args.oo,filename), 'w') as outfile:
    json.dump(output, outfile)


np.savetxt(os.path.join(args.oo, "confusionmatrix_new.txt"), results[1], delimiter=",", newline="\n")



print("....")

#save_np_txt(y_pred, preds_f)


imp_train=clf.interpret(train_dataloader,n_features=args.n_features,n_classes=args.n_classes)
imp_test =clf.interpret(test_dataloader,n_features=args.n_features,n_classes=args.n_classes)
imp =imp_train+imp_test
imp_file = os.path.join(args.oo, "attentionmap.txt")

save_np_txt(imp, imp_file)

print("Done")
