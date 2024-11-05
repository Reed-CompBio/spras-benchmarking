from pathlib import Path
from matplotlib import pyplot as plt
import pandas as pd
import networkx as nx
from sklearn.metrics import auc, precision_recall_curve, roc_curve

score_path = snakemake.input[0]
output_roc_path = snakemake.output[0]
output_pr_path = snakemake.output[1]

score_df = pd.read_csv(score_path, sep="\t")
y_score = list(score_df["y_score"])
y_true = list(score_df["y_true"])

# Compute ROC curve and PR curve
fpr, tpr, _ = roc_curve(y_true, y_score)
precision, recall, _ = precision_recall_curve(y_true, y_score)

# Compute AUC
roc_auc = auc(fpr, tpr)
pr_auc = auc(recall, precision)

# Plot ROC curve
plt.figure()
plt.plot(fpr, tpr, color='blue', lw=2, label='ROC curve (area = %0.2f)' % roc_auc)
plt.plot([0, 1], [0, 1], color='gray', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic')
plt.legend(loc="lower right")
plt.savefig(output_roc_path, format="png", dpi=300)
plt.close()

# Plot Precision-Recall curve
plt.figure()
plt.plot(recall, precision, color='blue', lw=2, label='PR curve (area = %0.2f)' % pr_auc)
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Precision-Recall Curve')
plt.legend(loc="lower left")
plt.savefig(output_pr_path, format="png", dpi=300)
plt.close()