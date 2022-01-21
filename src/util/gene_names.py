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

def get_symbol_from_refseq_name(refseq_names):
	# Symbol is like "RPL36B"
	mg = mygene.MyGeneInfo()
	resdf = mg.querymany(refseq_names, scopes='refseq', fields='all', species='all', as_dataframe=True)
	return list(resdf['symbol'].values)

def get_rpg_intron_mask(intron_set):
	intron_set.fill_ensembl_names()

	refseq_names = []
	for ii, intron in enumerate(intron_set.introns):
		if intron.name != "":
			refseq_names += [intron.name]
	symbols = get_symbol_from_refseq_name(refseq_names)

	is_rpg = []
	cur_idx = 0
	for intron in intron_set.introns:
		if intron.name == "":
			is_rpg += [False]
			continue
		symbol = symbols[cur_idx]
		cur_is_rpg = False
		if symbol[:2] == "RP":
			cur_is_rpg = True
		if symbol == "YML6" or symbol == "MRPL44":
			cur_is_rpg = True
		is_rpg += [cur_is_rpg]
		cur_idx += 1
	print(is_rpg)
	return is_rpg
