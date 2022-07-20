### **Fig. 6B, Fig. S12, Fig. S13, Fig. S14**

Violin plots comparing intron and control sequence sets on secondary structure metrics.

Run: `python analyses/plot_secstruct_features.py --make_violin_plots --intron_class_violin standard_allsize_min_50_max_600 --control_class_violin standard_allsize_min_50_max_600_shuffle --secstruct_package Vienna --secstruct_type ens`

Example output violin plot: 

<img width="326" alt="standard_shuffle_vienna_ens_4" src="https://user-images.githubusercontent.com/2606810/179938801-7d629136-8018-4641-8051-43f419a36b84.png">

Example output statistics (t-test and Wilcoxon ranked sum test comparisons): 
```
(288, 6)
(288, 6)
ZipperStemStartMetric_Vienna_ens: -21.611228-0.000000
Ttest_relResult(statistic=-10.476321303007937, pvalue=5.836862851981001e-22)
WilcoxonResult(statistic=1518.5, pvalue=4.710729677627628e-21)

ZipperStemEndMetric_Vienna_ens: -8.762750-0.000000
Ttest_relResult(statistic=-3.9562013322789515, pvalue=9.60272790339314e-05)
WilcoxonResult(statistic=4185.0, pvalue=0.00027499125398430265)

LongestStemMetric_Vienna_ens: 8.410000-35.783500
Ttest_relResult(statistic=4.306278769430331, pvalue=2.281656771266943e-05)
WilcoxonResult(statistic=16176.5, pvalue=0.0010599253627330688)

MLDMetric_Vienna_ens: 0.163169-0.391351
Ttest_relResult(statistic=5.70985901771819, pvalue=2.8212594171703126e-08)
WilcoxonResult(statistic=12699.5, pvalue=2.3123155263386365e-08)

LocalizationMetric_Vienna_ens: 7.594000-42.143500
Ttest_relResult(statistic=-6.564554894161535, pvalue=2.437312425260556e-10)
WilcoxonResult(statistic=12259.0, pvalue=1.508149693167277e-09)
```

Fig. 6B, Fig. S12, Fig. S13, and Fig. S14 were produced with this function using different intron classes, secondary structure prediction packages, and secondary structure prediction approaches. Below we list the four parameters (intron class, control class, secondary structure package, secondary structure prediction type) for each figure.

Fig. 6B: 
* Panel 1: (`standard_allsize_dms_matched_to_shift500`, `standard_allsize_dms_shift_500`, `RNAstructure_DMS`, `mfe`)
* Panel 2: (`standard_allsize_dms_matched_to_shift500`, `standard_allsize_dms_shift_500`, `RNAstructure`, `mfe`)
* Panel 3: (`standard_allsize_min_50_max_600`, `standard_allsize_min_50_max_600_shuffle`, `Vienna`, `ens`)

Fig. S12: 
* Panel 1: (`standard_allsize_dms_matched_to_shift500`, `standard_allsize_dms_shift_500`, `RNAstructure_DMS`, `mfe`)
* Panel 2: (`decoy_dms_matched_to_shift500`, `decoy_dms_shift_500`, `RNAstructure_DMS`, `mfe`)

Fig. S13:
* (`standard_allsize_min_50_max_600`, `standard_allsize_min_50_max_600_shuffle`, `Vienna`, `ens`)

Fig. S14: 
* A, Panel 1: (`standard_allsize_min_50_max_600`, `standard_allsize_min_50_max_600_shift_500`, `Vienna`, `ens`)
* A, Panel 2: (`standard_allsize_min_50_max_600`, `standard_allsize_min_50_max_600_shift_500_seq_matched`, `Vienna`, `ens`)
* A, Panel 3: (`standard_allsize_min_50_max_600`, `standard_allsize_min_50_max_600_shuffle`, `Vienna`, `ens`)
* B, Panel 1: (`standard_allsize_min_50_max_600`, `standard_allsize_min_50_max_600_shuffle`, `RNAstructure`, `ens`)
* B, Panel 2: (`standard_allsize_min_50_max_600`, `standard_allsize_min_50_max_600_shuffle`, `Vienna`, `ens`)
* C, Panel 1: (`standard_allsize_min_50_max_600_extend50`, `standard_allsize_min_50_max_600_extend50_shift_500`, `Vienna`, `ens`)


### **Fig. 6C**
Zipper stem and downstream stem dG values across species in the Saccharomyces genus, compared between intron and shuffled controls in violin plots. 

Run: `python analyses/plot_dg_zipper_stem_species.py standard_min_50_max_600 standard_min_50_max_600_shuffle --make_violin`

Example figure:

<img width="634" alt="stem_dg_start_species_noscer_3" src="https://user-images.githubusercontent.com/2606810/179938900-e9f22497-70de-4634-af8c-1457e3a99da3.png">

### **Fig. 6D**
Heatmap showing zipper stem conservation and per-species / per-intron statistics on zipper stem presence. 

Run: `python analyses/analyze_conservation.py ../database/alignments/hooks_alignments/ standard_min_50_max_600 --do_zipper_species`

Note: This code can also be used to produce heatmaps of zipper stem dG values in addition to binary presence vs absence, and can be used to produce heatmaps that separate ohnologous genes rather than combining them. 

Example figure: 

<img width="384" alt="Screen Shot 2022-07-20 at 1 34 16 AM" src="https://user-images.githubusercontent.com/2606810/179938653-faa1a613-ac58-4762-9f10-aea49da470a4.png">

### **Fig. S15A-B**

Heatmap summarizing log p-values for comparisons between intron and control sequence sets across secondary structure features and species.

Run: `python analyses/plot_secstruct_features.py --make_species_heatmap --intron_class_species standard_min_50_max_600 --control_class_species standard_min_50_max_600_phylo_control`

To compare against the shuffled control use the `control_class_species`: `standard_min_50_max_600_shuffle`. To compare against the phylogenetic control use the `control_class_species`: `standard_min_50_max_600_phylo_control`.

Example figure:

<img width="738" alt="feature_pval_heatmap_shuffle_log10_noscer_4" src="https://user-images.githubusercontent.com/2606810/179938508-d61e7c75-0592-4428-9fa2-b476ce595e9e.png">

To generate t-test and Wilcoxon ranked sum test statistics for these comparisons, run: `python analyses/stats_compare_species.py standard_min_50_max_600 standard_min_50_max_600_shuffle`

Example output statistics: 
```
Species: ndai
(218, 5)
(218, 5)
LocalizationMetric_Vienna_ens
Ttest_relResult(statistic=-3.585458673523364, pvalue=0.0004156450074144659)
WilcoxonResult(statistic=9081.0, pvalue=0.0022017715980185626)
ZipperStemStartMetric_Vienna_ens
Ttest_relResult(statistic=-6.423141467762445, pvalue=8.315522868368384e-10)
WilcoxonResult(statistic=2596.0, pvalue=1.9616570981257576e-09)
ZipperStemEndMetric_Vienna_ens
Ttest_relResult(statistic=2.3141468550248145, pvalue=0.02159515476604946)
WilcoxonResult(statistic=6277.0, pvalue=0.07561697190536668)
LongestStemMetric_Vienna_ens
Ttest_relResult(statistic=-2.060938743790534, pvalue=0.0405010858524741)
WilcoxonResult(statistic=8698.0, pvalue=0.000515909692312839)
MLDMetric_Vienna_ens
Ttest_relResult(statistic=4.215043482398518, pvalue=3.664594303046477e-05)
WilcoxonResult(statistic=7794.0, pvalue=8.915322004584187e-06)
```


### **Fig. S15C-D**

Barplots showing the number of orthologs for introns (full sequence and zipper stem regions), and the average sequence conservation for introns (full sequence and zipper stem regions). Statistics for these bar plots were generated through running the following commands: 

For Fig. S15C-D, zipper stem bars: 
Run: `python analyses/analyze_conservation.py ../database/alignments/hooks_alignments/ standard_allsize_min_50_max_600  --do_zipper`

For Fig. S15C-D, full intron bars: 
Run: `python analyses/analyze_conservation.py ../database/alignments/hooks_alignments/ standard_allsize_min_50_max_600  --do_stats`
