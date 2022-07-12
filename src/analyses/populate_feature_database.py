import pandas as pd 
import sys
from scipy import stats 
import numpy as np
from matplotlib import pyplot as plt
from util import features_db

pred_class = sys.argv[1]
package = sys.argv[2]
pred_approach = sys.argv[3]

secstruct_options = {'secstruct_pkg': package,
					'secstruct_type': pred_approach, 
					'use_bpp': False,
					'verbose': True,
					'force_eval': True
					}

all_features = ["LocalizationMetric", "ZipperStemStartMetric", "ZipperStemEndMetric", \
	"LongestStemMetric", "MLDMetric"]

feature_options_all = {}
for feature in all_features:
	feature_options_all[feature] = secstruct_options

standard_feature_df = features_db.get_features(all_features, \
	pred_class, feature_options_all=feature_options_all)
