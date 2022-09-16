import gzip
import html_to_json #
import json
import numpy #
import pandas #
import re

# -----=====| HGNC |=====-----
# https://www.genenames.org/download/custom/
# Last download: 2022-09-16

def DatabaseLinksFunc(Line):
	if Line != Line: return list()
	List = html_to_json.convert(Line.replace('<!--,-->', '<a></a>'))['a']
	Result = list()
	for Item in List:
		if not Item: continue
		else:
			if '_value' in Item: Name = Item['_value']
			else: Name = ' '.join(Item['_values'])
			Result.append({'name': Name, 'link': Item['_attributes']['href']})
	return Result

def IndexFunc(Line):
	if Line != Line: return list()
	Result = [i.strip() for i in Line.split(',')]
	return Result

def LoadJsonFunc(Line):
	if Line != Line: return list()
	Stripped = Line.strip('"')
	List = json.loads(f'["{Stripped}"]')
	return List

def LocusSpecificFunc(Line):
	if Line != Line: return list()
	Result = list()
	List = LoadJsonFunc(Line)
	for Item in List:
		Name, Link = Item.split('|')
		Result.append({ "name": str(Name), "link": str(Link) })
	return Result

def ListFunc(Line, sep = ','):
	if Line != Line: return list()
	List = [i.strip() for i in Line.split(sep)]
	Result = [k for k in List if k != ""]
	return Result

def IntListFunc(Line, sep = ','):
	if Line != Line: return list()
	Result = [int(i) for i in ListFunc(Line, sep=sep)]
	return Result

def GeneGroupIDFunc(Line): return IntListFunc(Line, sep = '|')

def GeneGroupNameFunc(Line): return ListFunc(Line, sep = '|')

def AsStrFunc(Line): return None if Line != Line else str(Line)

def ParserHGNC(TableHGNC, DatabaseFileHGNC):
	Func = {
		'HGNC ID': IndexFunc,
		'Approved symbol': IndexFunc,
		'Previous symbols': IndexFunc,
		'Alias symbols': IndexFunc,
		'RefSeq IDs': IndexFunc,
		'Ensembl gene ID': IndexFunc,
		'CCDS IDs': IndexFunc,
		'Vega IDs': IndexFunc,
		'RefSeq(supplied by NCBI)': IndexFunc,
		'Ensembl ID(supplied by Ensembl)': IndexFunc,
		'Vega ID(supplied by Vega)': IndexFunc,
		'UCSC ID(supplied by UCSC)': IndexFunc,
		'LNCipedia ID (supplied by LNCipedia)': IndexFunc,
		'GtRNAdb ID (supplied by GtRNAdb)': IndexFunc,
		'AGR HGNC ID (supplied by Alliance of Genomic Resources)': IndexFunc,
		'Approved name': AsStrFunc,
		'Status': AsStrFunc,
		'Chromosome': AsStrFunc,
		'Accession numbers': ListFunc,
		'Locus type': AsStrFunc,
		'Locus group': AsStrFunc,
		'Previous name': LoadJsonFunc,
		'Alias names': LoadJsonFunc,
		'Date approved': AsStrFunc,
		'Date modified': AsStrFunc,
		'Date symbol changed': AsStrFunc,
		'Date name changed': AsStrFunc,
		'Enzyme IDs': ListFunc,
		'NCBI Gene ID': IntListFunc,
		'Mouse genome database ID': ListFunc,
		'Specialist database links': DatabaseLinksFunc,
		'Specialist database IDs': ListFunc,
		'Pubmed IDs': IntListFunc,
		'Gene group ID': GeneGroupIDFunc,
		'Gene group name': GeneGroupNameFunc,
		'Locus specific databases': LocusSpecificFunc,
		'NCBI Gene ID(supplied by NCBI)': IntListFunc,
		'OMIM ID(supplied by OMIM)': IntListFunc,
		'UniProt ID(supplied by UniProt)': ListFunc,
		'Mouse genome database ID(supplied by MGI)': ListFunc,
		'Rat genome database ID(supplied by RGD)': ListFunc,
		'MANE Select Ensembl transcript ID (supplied by NCBI)': ListFunc,
		'MANE Select RefSeq transcript ID (supplied by NCBI)': ListFunc
		}
	Result = { 'FullInfo': None, 'Index': {} }
	Data = pandas.read_csv(TableHGNC, sep = '\t', dtype = str, quotechar='№')
	Data['Index'] = Data['Approved symbol'].apply(str)
	for Col, ApplyFunc in Func.items():
		Data[Col] = Data[Col].apply(ApplyFunc)
		if ApplyFunc == IndexFunc:
			for _, Line in Data[['Index', Col]].iterrows():
				for Code in Line[Col]:
					if (Code in Result['Index']):
						if (Line['Index'] not in Result['Index'][Code]): Result['Index'][str(Code)].append(str(Line['Index']))
					else: Result['Index'][str(Code)] = [str(Line['Index'])]
	Data['Approved symbol'] = Data['Index'].apply(str)
	Data = Data.drop(columns=['Index'])
	Result['FullInfo'] = Data.groupby('Approved symbol').apply(lambda x: list(x.transpose().to_dict().values())).to_dict()
	json.dump(Result, gzip.open(DatabaseFileHGNC, 'wt'), separators = (',', ':'))
