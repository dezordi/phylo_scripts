Phylo_scripts
=========

Esse repositório contém um conjunto de scripts para automatizar a formatação de arquivos usados em análises filogenéticas

=====
gisaid_loc_beauti.py
=====

Esse script produz um arquivo de latitude e longitude com 3 colunas: traits, lat, long. Utilizado para criação do arquivo xml na ferramenta BEAUTi para análises filogeográficas. Foi criado para seguir o padrão de nomenclatura presente no gisaid. Para uso, são necessários 3 arquivos de input, utilizado nos 3 argumentos:

* -fa: Arquivo fasta, com o alinhamento das sequências
* -mt: Arquivo de metadados, que pode ser baixado diretamente no meno download do gisaid, seguindo o padrão do gisaid com as seguintes colunas: "'Virus name','Type','Accession ID','Collection date','Location','Additional location information','Sequence length','Host','Patient age','Gender','Clade','Pango lineage','Pangolin version','Variant','AA Substitutions','Submission date','Is reference?','Is complete?','Is high coverage?','Is low coverage?','N-Content','GC-Content'"
* -ct: Arquivo de latitude e longitude por município brasileiro, seguindo o padrão `deste arquivo <https://github.com/kelvins/Municipios-Brasileiros/blob/main/csv/municipios.csv>`_

Serão retornados 3 arquivos de output:

* _coordinates.fasta: Um alinhamento apenas com as sequências que possuem informação de localização a nível de município.
* _coordinates.tsv: Um tabular .tsv no formato a ser carregado no BEAUti, com o nome do genoma, e as coordenadas de latitude e longitude do município correspondente.
* _without_coordinates.tsv: Um tabular .csv com uma série de informações sobre os genomas que não possuem informação a nível de município.

Aviso: Caso o município esteja escrito de forma incorreta no gisaid, a busca por coordenadas falhará.

Caso necessário, um abiente conda pode ser importado para execução do script: 

.. code:: bash
    conda env create -f envs/gisaid_loc_beast.yml
    conda activate gisaid_loc_beast
    python gisaid_loc_beauti.py -fa <alinhamento.fasta> -mt <metadados.tsv> -ct <municipios.csv> 

