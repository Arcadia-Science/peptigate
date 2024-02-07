```
mkdir sandbox
cd sandbox
mkdir try_tools
cd try_tools
conda create -n pepsand rnasamba cpc2 cpat borf pip
conda activate pepsand
```

Download files from S3 to test software on.
* sORF prediction from small contigs: `s3://arcadia-github-repository-outputs/2023-amblyomma-americanum-txome-assembly/outputs/assembly/small_contigs.fa`
* cleavage peptide prediction from ORF predictions: `s3://arcadia-github-repository-outputs/2023-amblyomma-americanum-txome-assembly/outputs/annotation/dammit/orthofuser_final_clean.fa.transdecoder.pep`
* nucleotide sequences for eventual clustering, dn/ds estimation, etc: `s3://arcadia-github-repository-outputs/2023-amblyomma-americanum-txome-assembly/outputs/annotation/dammit/orthofuser_final_clean.fa.transdecoder.cds`
* **note** we need full length transcripts (long), these are also in the dammit folder. I don't think any tool take those as input tho, so proceeding with other files.

or, on AWS instance where amblyomma txome work was done:
```
cp ~/2023-amblyomma-americanum-txome-assembly/outputs/annotation/dammit/orthofuser_final_clean.fa.dammit.fasta .
cp ~/2023-amblyomma-americanum-txome-assembly/outputs/annotation/dammit/orthofuser_final_clean.fa.transdecoder.pep .
cp ~/2023-amblyomma-americanum-txome-assembly/outputs/annotation/dammit/orthofuser_final_clean.fa.transdecoder.cds
cp ~/2023-amblyomma-americanum-txome-assembly/outputs/annotation/dammit/orthofuser_final_clean.fa.transdecoder.cds .
cp ~/2023-amblyomma-americanum-txome-assembly/outputs/assembly/small_contigs.fa .
```

## sORF prediction

### download databases of small ORFs to see size distribution

```
mkdir smprot
cd smprot
wget http://bigdata.ibp.ac.cn/SmProt/datadownload/SmProt_LiteratureMining.fa.gz
wget http://bigdata.ibp.ac.cn/SmProt/datadownload/SmProt_KnownDatabase.fa.gz
http://bigdata.ibp.ac.cn/SmProt/datadownload/SmProt_RibosomeProfiling.fa.gz
wget http://bigdata.ibp.ac.cn/SmProt/datadownload/SmProt_MS.fa.gz
```

```
mamba install seqkit csvtk
seqkit stats *.fa.gz -T | csvtk csv2md -t
```

|file                          |format|type   |num_seqs|sum_len|min_len|avg_len|max_len|Q1  |Q2  |Q3  |sum_gap|N50|Q20(%)|Q30(%)|AvgQual|GC(%)|
|:-----------------------------|:-----|:------|:-------|:------|:------|:------|:------|:---|:---|:---|:------|:--|:-----|:-----|:------|:----|
|SmProt_KnownDatabase.fa.gz    |FASTA |Protein|2998    |232656 |4      |77.6   |100    |68.0|82.0|92.0|0      |85 |0.00  |0.00  |0.00   |10.44|
|SmProt_LiteratureMining.fa.gz |FASTA |Protein|81454   |3265961|2      |40.1   |100    |17.0|32.0|60.0|53     |59 |0.00  |0.00  |0.00   |8.93 |
|SmProt_MS.fa.gz               |FASTA |Protein|117099  |1759228|5      |15.0   |52     |12.0|14.0|18.0|0      |15 |0.00  |0.00  |0.00   |8.36 |
|SmProt_RibosomeProfiling.fa.gz|FASTA |Protein|53459   |2319004|2      |43.4   |100    |21.0|39.0|64.0|0      |62 |0.00  |0.00  |0.00   |7.69 |

### BORF

```
borf -l 15 -s -o borf/small_contigs small_contigs.fa
```

### CPAT

```
mkdir cpat
cd cpat
# download logit model
wget -O Flyl_logitModel.RData https://sourceforge.net/projects/rna-cpat/files/v1.2.2/prebuilt_model/Fly_logitModel.RData/download

# download hexamer table
wget -O fly_Hexamer.tsv https://sourceforge.net/projects/rna-cpat/files/v1.2.2/prebuilt_model/fly_Hexamer.tsv/download

# -x hexamer fequency table. prebuilt for fly
# -d logit model prebuilt for fly as closest relative to chelicerate
# --min-orf 15 ## ASK WHAT I SHOULD USE AS MINIMUM FOR sORF IDENTIFICATION
# --antisense
# -g RNA sequences
# -o out file prefix
cpat.py -x fly_Hexamer.tsv -d Fly_logitModel.RData --min-orf 15 --antisense -g ../small_contigs.fa -o small_contigs_sorfs
```

### rnasamba

All three of the below installation methods failed on Mac Arm64 architecture, even when running rosetta. 
The conda and pip methods failed with `[6]    582 illegal hardware instruction`
```
conda create -n rnasamba rnasamba
conda activate rnasamba
```

```
conda create -n rnasambapip python pip tensorflow=1.14.0
conda activate rnasambapip
pip install setuptools-rust
pip install rnasamba==0.2.5
```

The docker method failed with `qemu: uncaught target signal 6 (Aborted) - core dumped`
```
docker run -ti --rm -u $(id -u) -v "$(pwd):/app" antoniopcamargo/rnasamba classify -p predicted_proteins.fa classification.tsv ../small_contigs.fa full_length_weights.hdf5
docker run --platform linux/x86_64 -ti --rm -u $(id -u) -v "$(pwd):/app" antoniopcamargo/rnasamba classify -p predicted_proteins.fa classification.tsv ../small_contigs.fa full_length_weights.hdf5
```

Switched dev over to AWS.

```
curl -JLO https://raw.githubusercontent.com/apcamargo/RNAsamba/master/data/partial_length_weights.hdf5
curl -JLO https://raw.githubusercontent.com/apcamargo/RNAsamba/master/data/full_length_weights.hdf5
rnasamba classify -p predicted_proteins_full_length_weights.fa classification_full_length_weights.tsv ../small_contigs.fa full_length_weights.hdf5
rnasamba classify -p predicted_proteins_partial_length_weights.fa classification_partial_length_weights.tsv ../small_contigs.fa partial_length_weights.hdf5
```

### cpc2
```
mkdir cpc2
cd cpc2
CPC2.py
```

## Cleavage tools

### nlpprecursor

```
conda create -n nlpprecursor python=3.7
conda activate nlpprecursor
pip install torch==1.0.0
pip install git+https://github.com/nishanthmerwin/fastai.git@deepripp_install
pip install git+https://github.com/magarveylab/nlpprecursor
mamba install biopython
```

download  models
```
curl -JLO https://github.com/magarveylab/NLPPrecursor/releases/download/1.0/nlpprecursor_models.tar.gz
tar xf nlpprecursor_models.tar.gz
```
make a small protein set
```
head -n 32 ../orthofuser_final_clean.fa.transdecoder.pep > head.pep
sed '/^[^>]/s/\*//g' head.pep > head.noasterisk.pep

```
```
python run_nlpprecursor.py models head.noasterisk.pep out.tsv
```

### deeppetide

prepare sequences by removing transdecoder's stop codon asterisk
```
sed '/^[^>]/s/\*//g' orthofuser_final_clean.fa.transdecoder.pep > orthofuser_final_clean.fa.transdecoder.noasterisk.pep
```

Installed via git clone, will need to build a docker container and put on docker hub for use with a singularity environment for snakemake.

```
python3 predict.py -ff orthofuser_final_clean.fa.transdecoder.noasterisk.pep -od testrun_outputs2/ --batch_size 1000
```

if not sped up by GPU & increased batch size, [split file and run in separate deeppeptide runs](https://github.com/fteufel/DeepPeptide/issues/2).
