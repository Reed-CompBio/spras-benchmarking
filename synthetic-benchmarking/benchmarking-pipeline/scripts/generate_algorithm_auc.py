from pathlib import Path
from matplotlib import pyplot as plt
import pandas as pd
import networkx as nx
from sklearn.metrics import auc, precision_recall_curve, roc_curve

score_path = snakemake.input[0]
output_roc_path = snakemake.output[0]
output_pr_path = snakemake.output[1]

f = open(output_roc_path, "w+")
f = open(output_pr_path, "w+")