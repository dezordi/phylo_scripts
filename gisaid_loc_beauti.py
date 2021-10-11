#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Filipe Zimmer Dezordi"
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Filipe Zimmer Dezordi"
__email__ = "zimmer.filipe@gmail.com"
__status__ = "Development"
__username__ = 'dezordi'

import argparse, re, csv
from Bio import SeqIO
import pandas as pd
import numpy as np

###############################>ARGUMENTS<###############################
parser = argparse.ArgumentParser(description = 'Esse script cria um arquivo _coordinates.txt com o nome da sequencia, latitude e longitude para a preparação do arquivo xml com a ferramenta BEAUti. Foi criado para trabalhar com o formato de dados do EpiCoV-gisaid.')
parser.add_argument("-fa", "--fasta_file", help="Arquivo Fasta.", required=True)
parser.add_argument("-mt", "--meta_data", help="Arquivo de metadados do gisaid.", required=True)
parser.add_argument("-ct", "--city_csv", help="Arquivo csv com informações de cidades brasileiras.", required=True)

args = parser.parse_args()
fasta_file = args.fasta_file
meta_data_file = args.meta_data
city_csv_file = args.city_csv

fasta_epi_dict = dict() 
loc_epi_dict = dict()
code_to_uf = {11:"rondonia",
12:"acre",
13:"amazonas",
14:"roraima",
15:"para",
16:"amapa",
17:"tocantins",
21:"maranhao",
22:"piaui",
23:"ceara",
24:"rio_grande_do_norte",
25:"paraiba",
26:"pernambuco",
27:"alagoas",
28:"sergipe",
29:"bahia",
31:"minas_gerais",
32:"espirito_santo",
33:"rio_de_janeiro",
35:"sao_paulo",
41:"parana",
42:"santa_catarina",
43:"rio_grande_do_sul",
50:"mato_grosso_do_sul",
51:"mato_grosso",
52:"goias",
53:"distrito_federal"}

#convertendo os headers do fasta em um dicionario epi:header
for record in SeqIO.parse(fasta_file, "fasta"):
    epi_code = re.sub(r'.*\|EPI','EPI',record.id).rstrip('\n')
    epi_code = re.sub(r'\|.*','',epi_code)
    fasta_epi_dict[epi_code] = record.id

#convertendo as colunas Acession e Localization em um dicionario com epi:localizacao
with open(meta_data_file,'r') as meta_data_file_reader:
    meta_data_file_reader_csv = csv.reader(meta_data_file_reader, delimiter='\t')
    for line in meta_data_file_reader_csv:
        if line[2] in fasta_epi_dict:
            loc_epi_dict[line[2]] = line[4]

#transformando os dicionarios em dataframes
fasta_epi_df = pd.DataFrame(fasta_epi_dict.items())
fasta_epi_df.columns = ['epi','header']
loc_epi_df = pd.DataFrame(loc_epi_dict.items())
loc_epi_df.columns = ['epi','loc']

#concatenando os dataframes baseado no codigo epi
fasta_loc_epi_df = pd.merge(fasta_epi_df, loc_epi_df, on = 'epi')

#limpando a memoria
del fasta_epi_dict, fasta_epi_df, loc_epi_dict, loc_epi_df, meta_data_file_reader, meta_data_file_reader_csv

#criando o dataframe para armazenar a localizacao normalizada para cada sequencia
fasta_loc_epi_df[['continente','pais','codigo_uf','nome']] = fasta_loc_epi_df['loc'].str.split(' / ', expand = True)
fasta_loc_epi_df = fasta_loc_epi_df.replace(r'^\s*$',np.NaN, regex=True)
fasta_loc_epi_df['codigo_uf'] = fasta_loc_epi_df['codigo_uf'].str.normalize('NFKD').str.encode('ascii', errors='ignore').str.decode('utf-8').str.replace(' ','_').str.lower()
fasta_loc_epi_df['nome'] = fasta_loc_epi_df['nome'].str.normalize('NFKD').str.encode('ascii', errors='ignore').str.decode('utf-8').str.replace(' ','_').str.lower()

#convertendo os uf_codes para o nome do estado
city_csv_file_df = pd.read_csv(city_csv_file, header=0)
for codigo, estado in code_to_uf.items():
    city_csv_file_df.loc[city_csv_file_df.codigo_uf == codigo, "codigo_uf"] = estado
city_csv_file_df['nome'] = city_csv_file_df['nome'].str.normalize('NFKD').str.encode('ascii', errors='ignore').str.decode('utf-8').str.replace(' ','_').str.lower()

#selecionando as sequencias com latitude e longitude
fasta_epi_coord = pd.merge(fasta_loc_epi_df, city_csv_file_df, on=['codigo_uf','nome'])
fasta_epi_coord = fasta_epi_coord.round({'latitude':4, 'longitude':4})
beauti_traits = fasta_epi_coord[['header','latitude','longitude']]
beauti_traits = beauti_traits.rename(columns ={'header':'traits','latitude':'lat','longitude':'long'})

#criando o arquivo de coordenadas para o beauti
beauti_traits.to_csv(fasta_file+'_coordinates.tsv',
                        index=False,
                        header=True,
                        sep='\t')

#criando um aruqivo de alinhamento apenas com as sequencias com latitude e longitude
header_with_loc_list = fasta_epi_coord['header'].tolist()
fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file,'fasta'))
seqs_com_loc = [fasta_dict[seq_id] for seq_id in header_with_loc_list]
SeqIO.write(seqs_com_loc, fasta_file+'_coordinates.fasta','fasta')

#criando um arquivo com a localizao das sequencias que nao tem cidade
epi_with_loc = fasta_epi_coord['epi'].tolist()
epi_without_loc = fasta_loc_epi_df[~fasta_loc_epi_df.epi.isin(epi_with_loc)]
epi_without_loc.to_csv(fasta_file+'_without_coordinates.tsv',
                        index=False,
                        header=True,
                        sep='\t')