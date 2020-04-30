import pandas as pd 
import sys
from scipy import stats 
import numpy as np
from matplotlib import pyplot as plt

f = open(sys.argv[1])
features = pd.read_csv(f)
f.close()

control_f = open(sys.argv[2])
control_features = pd.read_csv(control_f)
control_f.close()

metrics = ["LocalizationMetric", "StartToBPStemMetric", "BPToEndStemMetric", "StartProtectionMetric", 
		"EndProtectionMetric", "BPProtectionMetric"]

ii = 1
for metric in metrics:
	intron_vals = []
	control_vals = []
	for idx, row in features.iterrows():
		intron_vals += [row[metric]]
	for idx, row in control_features.iterrows():
		control_vals += [row[metric]]
	plt.subplot(1, len(metrics), ii)
	ii += 1
	plt.hist(np.array(intron_vals) - np.array(control_vals), alpha=0.5, rwidth=0.85, color="forestgreen")
	plt.axvline(x=0, linestyle='--', color='black')
	print(metric)
	print(stats.ttest_rel(np.array(intron_vals), np.array(control_vals)))
	print(stats.wilcoxon(np.array(intron_vals), np.array(control_vals)))

plt.show()
