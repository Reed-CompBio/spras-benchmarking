import os
from pathlib import Path
from matplotlib import pyplot as plt
import pandas as pd
import networkx as nx
from sklearn.metrics import auc, precision_recall_curve, roc_curve

input_files = snakemake.input.scores
output_pr_path = snakemake.output[0]

print(output_pr_path.split("/")[-1].split("_")[-1].split(".")[0])
auc_type = output_pr_path.split("/")[-1].split("_")[-1].split(".")[0]
print("")

algorithm = output_pr_path.split("/")[-1].split("_")[0]

filtered_inputs = []
i = 0
for input in input_files:
    current_alg = ""
    if auc_type == "alg":
        current_alg = input.split("/")[-1].split("_")[0]
    else:
        current_alg = input.split("/")[-1].split("_")[-1].split(".")[0]
    if current_alg == algorithm:
        filtered_inputs.append(input)

for input in input_files:
    print(input)

print("")

for input in filtered_inputs:
    print(input)

print("")
print(output_pr_path)


plt.figure(figsize=(10, 5))
colors = ["b", "g", "r", "c", "m", "y", "k", "orange"]

for idx, file in enumerate(filtered_inputs):
    data = []
    y_true = []
    y_score = []
    with open(file, "r") as f:
        next(f)
        for line in f:
            line = line.strip()
            col = line.split("\t")
            y_true.append(int(col[2]))
            y_score.append(int(col[3]))

    precision, recall, _ = precision_recall_curve(y_true, y_score)
    fpr, tpr, _ = roc_curve(y_true, y_score)

    plt.subplot(1, 2, 1)
    plt.plot(
        recall, precision, color=colors[idx % len(colors)], label=os.path.basename(file)
    )

    plt.subplot(1, 2, 2)
    plt.plot(fpr, tpr, color=colors[idx % len(colors)], label=os.path.basename(file))

plt.subplot(1, 2, 1)
plt.xlabel("Recall")
plt.ylabel("Precision")
plt.title("Precision-Recall Curve")
plt.legend()

plt.subplot(1, 2, 2)
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC Curve")
plt.legend()

plt.tight_layout()
plt.savefig(output_pr_path, format="png", dpi=300)
