import mygene
import pandas as pd

def get_ensembl_names(refseq_names):
	all_empty = True
	for refseq_name in refseq_names:
		if refseq_name != "":
			all_empty = False
	if all_empty:
		return refseq_names
	mg = mygene.MyGeneInfo()
	resdf = mg.querymany(refseq_names, scopes='refseq', fields='ensembl.gene', species='all', as_dataframe=True)
	return list(resdf['ensembl.gene'].values)

def get_ensemble_name_from_symbol(symbol):
	# Symbol is like "RPL36B"
	mg = mygene.MyGeneInfo()
	resdf = mg.querymany([symbol], scopes='symbol', fields='all', species='all', as_dataframe=True)
	return list(resdf['ensembl.gene'])