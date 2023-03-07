# PrediXcan Workflow

## Directory Structure

```bash
PrediXcan
├── PrediXcan.py
├── PrediXcan_Workflow.md
├── convert_plink_to_dosage.py
├── genotype
│   └── chr{1..23}.txt.gz
└── results
    ├── ensemble2gene-tab.map
    └── map_genes.R
    └── {*}_predicted_expression.txt
```

## Converting Plink Files to Dosage

Modified the provided by the authors for `python3` support. Original script is available here:[GitHub](https://github.com/hakyimlab/PrediXcan/tree/master/Software).

```bash
python3 convert_plink_to_dosage.py \ 
        -b <file_basename> \
        -o <output_directory/prefix> \
        -p <path_to_plink_executable>
```

`<file_basename>` is the prefix of the `.bim`, `.bed`, `.fam` files. `<path_to_plink_executable>` can be found by typing `which plink` in shell. `output_directory` is `genotype` in the structure above.

For example,

```bash
python3 convert_plink_to_dosage.py \
        -b racer-mega-imputed-f-matched-v3-QC \ 
        -o genotype1/chr \ 
        -p /ihome/crc/install/plink/1.90_20230116/plink
```


## Running PrediXcan

```bash
python3 PrediXcan.py \
        --predict \
        --dosages <output_directory> \
        --dosages_prefix <prefix> \
        --samples <sample_file> \
        --weights <weights_file.db> \
        --output_prefix <predix_output_prefix>
```

`output_directory` is the same as above. Do not mention the full path for the sample file, only mention the file name. It should be stored in `output_directory` along with the dosage files.

The weights can be downloaded from PredictDB [here](https://predictdb.org/post/2021/07/21/gtex-v8-models-on-eqtl-and-sqtl/). 

Direct link to tar with files for PrediXcan’s support on expression is [here](https://zenodo.org/record/3518299/files/mashr_eqtl.tar?download=1)

For example,
```bash
python3 PrediXcan.py \
        --predict \ 
        --dosages ./genotype \ 
        --dosages_prefix chr \ 
        --samples racer-mega-imputed-f-matched-v3-QC.fam \ 
        --weights en_Whole_Blood.db \ 
        --output_prefix results/en_Whole_Blood
```




## Mapping Gene Names from Ensembl IDs(ENSG...) to Names

```bash
Rscript map_genes.R <predicted_expression_file>
```

For example,

```bash
Rscript map_genes.R en_Whole_Blood_predicted_expression.txt
```


## Combined Code


```bash
mkdir genotype
python3 convert_plink_to_dosage.py \
        -b racer-mega-imputed-f-matched-v3-QC \ 
        -o genotype1/chr \ 
        -p /ihome/crc/install/plink/1.90_20230116/plink


python3 PrediXcan.py \
        --predict \ 
        --dosages ./genotype \ 
        --dosages_prefix chr \ 
        --samples racer-mega-imputed-f-matched-v3-QC.fam \ 
        --weights en_Whole_Blood.db \ 
        --output_prefix results/en_Whole_Blood

cd results
Rscript map_genes.R en_Whole_Blood_predicted_expression.txt     
```