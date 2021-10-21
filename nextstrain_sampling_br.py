#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Filipe Zimmer Dezordi"
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Filipe Zimmer Dezordi"
__email__ = "zimmer.filipe@gmail.com"
__status__ = "Development"
__username__ = 'dezordi'

import argparse, subprocess, shlex
from math import trunc


###############################>ARGUMENTS<###############################
parser = argparse.ArgumentParser(description = 'Esse script automatiza as analises do nextstrain, criado para a rotina de análise de SARS-CoV-2 da rede')
parser.add_argument("-fa", "--fasta_file", help="Arquivo com os genomas do gisaid", required=True)
parser.add_argument("-dr", "--ncov_dir", help="caminho para o diretório do ncov", required=True)
parser.add_argument("-mt", "--meta_data", help="Arquivo com os metadados do gisaid", required=True)
parser.add_argument("-st", "--br_state", help="Estado alvo do sampling", required=True)
parser.add_argument("-sp", "--sampling", help="Fazer um subsampling do estado?")
parser.add_argument("-li", "--lineage_data", help="Arquivo com as linhagens desejadas, uma linhagem por linha", required=True)
parser.add_argument("-n", "--number_genomes", help="Número de genomas no Brasil (exluindo estado) e no mundo", default=4000)
parser.add_argument("-le", "--genome_length", help="Filtro de tamanho dos genomas, default 28400, o que representa pelo menos 95'%' de cobertura do genoma do SARS-CoV-2", default=28400)

args = parser.parse_args()

def get_br(metadata,index,length,lineage,number,region,output):
    print("="*10,f"Step4.1: Create Brazilian without {region} state subsampling","="*10)
    sampling_br = f'''augur filter \
    --metadata {metadata} \
    --sequence-index {index} \
    --min-length {str(length)} \
    --query "(country == 'Brazil') & (pango_lineage == {lineage}) & (division != '{region}')" \
    --subsample-max-sequences {str(number)}  \
    --exclude-ambiguous-dates-by any \
    --group-by division pango_lineage year month \
    --output-strains {output}'''
    sampling_br = shlex.split(sampling_br)
    cmd_sampling_br = subprocess.Popen(sampling_br)
    cmd_sampling_br.wait()
    print("="*10,"Step4.1: Done","="*10)

def get_region(metadata,index,length,lineage,region,sampling,output):
    print("="*10,f"Step4.2: Create {region} state subsampling","="*10)
    if sampling == None:
        sampling_state = f'''augur filter \
        --metadata {metadata} \
        --sequence-index {index} \
        --min-length {str(length)} \
        --query "(country == 'Brazil') & (pango_lineage == {lineage}) & (division == '{region}')" \
        --exclude-ambiguous-dates-by any \
        --group-by location pango_lineage year month 
        --output-strains {output}'''
    else:
        sampling_state = f'''augur filter \
        --metadata {metadata} \
        --sequence-index {index} \
        --min-length {str(length)} \
        --query "(country == 'Brazil') & (pango_lineage == {lineage}) & (division == '{region}')" \
        --subsample-max-sequences {str(sampling)}
        --exclude-ambiguous-dates-by any \
        --group-by location pango_lineage year month 
        --output-strains {output}'''
    sampling_state = shlex.split(sampling_state)
    cmd_sampling_state = subprocess.Popen(sampling_state)
    cmd_sampling_state.wait()
    print("="*10,"Step4.2: Done","="*10)

def get_global(metadata,index,length,lineage,number,output):
    print("="*10,f"Step4.3: Create global subsampling","="*10)
    sampling_global = f'''augur filter \
    --metadata {metadata} \
    --sequence-index {index} \
    --min-length {str(length)} \
    --query "(country != 'Brazil') & (pango_lineage == {lineage})" \
    --subsample-max-sequences {str(number)} 
    --exclude-ambiguous-dates-by any \
    --group-by country pango_lineage year month \
    --output-strains {output}'''
    sampling_global = shlex.split(sampling_global)
    cmd_sampling_global = subprocess.Popen(sampling_global)
    cmd_sampling_global.wait()
    print("="*10,"Step4.3: Done","="*10)

def get_outgroup(metadata,index,length,output):
    print("="*10,f"Step4.4: Create A and B lineages subsampling","="*10)
    sampling_groupA = f'''augur filter \
    --metadata {metadata} \
    --sequence-index {index} \
    --min-length {str(length)} \
    --query "pango_lineage == ['A','B'] " \
    --subsample-max-sequences 30 \
    --exclude-ambiguous-dates-by any \
    --group-by country pango_lineage year month \
    --output-strains {output}'''
    sampling_groupA = shlex.split(sampling_groupA)
    cmd_sampling_groupA = subprocess.Popen(sampling_groupA)
    cmd_sampling_groupA.wait()
    print("="*10,"Step4.4: Done","="*10)

with open(args.lineage_data, 'r') as lineage_data:
    lineage_list = [line.strip('\n') for line in lineage_data.readlines()]

##sanitize sequences
print("="*40,"Step1: Sanitize sequences","="*40)
sanitize_sequences = f"python {args.ncov_dir}scripts/sanitize_sequences.py \
    --sequences {args.fasta_file} \
    --strip-prefixes 'hCoV-19/' \
    --output {args.ncov_dir}data/sequences_gisaid.fasta.gz"
sanitize_sequences = shlex.split(sanitize_sequences)
cmd_sanitize_sequences = subprocess.Popen(sanitize_sequences)
cmd_sanitize_sequences.wait()
print("="*40,"Step1: Done","="*40)

##index sequences
print("="*40,"Step2: Index sequences","="*40)
index_sequences = f"augur index \
    --sequences {args.ncov_dir}data/sequences_gisaid.fasta.gz \
    --output {args.ncov_dir}data/sequence_index_gisaid.tsv.gz"
index_sequences = shlex.split(index_sequences)
cmd_index_sequences = subprocess.Popen(index_sequences)
cmd_index_sequences.wait()
print("="*40,"Step2: Done","="*40)

##sanitize metadata
print("="*40,"Step3: Sanitize metadata","="*40)
sanitize_metadata = f"python {args.ncov_dir}scripts/sanitize_metadata.py \
    --metadata {args.meta_data} \
    --parse-location-field Location \
    --rename-fields 'Virus name=strain' 'Accession ID=gisaid_epi_isl' 'Collection date=date' 'Pango lineage=pango_lineage' \
    --strip-prefixes 'hCoV-19/' \
    --output {args.ncov_dir}data/metadata_gisaid.tsv.gz"
sanitize_metadata = shlex.split(sanitize_metadata)
cmd_sanitize_metadata = subprocess.Popen(sanitize_metadata)
cmd_sanitize_metadata.wait()
print("="*40,"Step3: Done","="*40)

###filter
print("="*40,"Step4: Filter","="*40)
get_br(f"{args.ncov_dir}data/metadata_gisaid.tsv.gz",
       f"{args.ncov_dir}data/sequence_index_gisaid.tsv.gz",
       args.genome_length,
       lineage_list,trunc(int(args.number_genomes)*0.5),
       args.br_state, 
       f"{args.ncov_dir}data/br_sub.txt")
get_region(f"{args.ncov_dir}data/metadata_gisaid.tsv.gz",
           f"{args.ncov_dir}data/sequence_index_gisaid.tsv.gz",
           args.genome_length,
           lineage_list,
           args.br_state,
           args.sampling,
           f"{args.ncov_dir}data/br_state_sub.txt")
get_global(f"{args.ncov_dir}data/metadata_gisaid.tsv.gz",
           f"{args.ncov_dir}data/sequence_index_gisaid.tsv.gz",
           args.genome_length,
           lineage_list,
           trunc(int(args.number_genomes)*0.5), 
           f"{args.ncov_dir}data/global_sub.txt")
get_outgroup(f"{args.ncov_dir}data/metadata_gisaid.tsv.gz",
           f"{args.ncov_dir}data/sequence_index_gisaid.tsv.gz",
           args.genome_length,
           f"{args.ncov_dir}data/outgroup_sub.txt")
print("="*40,"Step4: Done","="*40)

###extract metadata and sequences
print("="*40,"Step5: Extract samples","="*40)
extract_sequences = f"augur filter \
    --metadata {args.ncov_dir}data/metadata_gisaid.tsv.gz \
    --sequence-index {args.ncov_dir}data/sequence_index_gisaid.tsv.gz \
    --sequences {args.ncov_dir}data/sequences_gisaid.fasta.gz \
    --exclude-all \
    --include {args.ncov_dir}data/br_sub.txt {args.ncov_dir}data/br_state_sub.txt {args.ncov_dir}data/global_sub.txt {args.ncov_dir}data/outgroup_sub.txt \
    --output-metadata {args.ncov_dir}data/subsampled_metadata_gisaid.tsv.gz \
    --output-sequences {args.ncov_dir}data/subsampled_sequences_gisaid.fasta.gz"
extract_sequences = shlex.split(extract_sequences)
cmd_extract_sequences = subprocess.Popen(extract_sequences)
cmd_extract_sequences.wait()
print("="*40,"Step5: Done","="*40)