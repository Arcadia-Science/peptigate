import pandas as pd
import os
import re

metadata = pd.read_csv("inputs/models/datasets/train_data_links.tsv", sep="\t")
GENOMES = metadata["genome"].unique().tolist()
RNA_TYPES = ["cdna", "ncrna"]  # inherits names from ensembl
VALIDATION_TYPES = [
    "mRNAs",
    "ncRNAs",
]  # inherits names from https://github.com/cbl-nabi/RNAChallenge
# Note that for variables CODING_TYPES and DATASET_TYPES, the order matters;
# The snakemake expand() function maintains the order of these arguments.
# This behavior is used later in the snakefile to specify inputs and outputs of different scripts/shell commands.
# The order of CODING_TYPES must correspond to the order of the positional args that are later passed to rnasamba.
# In addition, the order of CODING_TYPES and DATASET_TYPES is used to declare outputs in the rule/R script process_sequences_into_nonoverlapping_sets.
CODING_TYPES = ["coding", "noncoding"]
DATASET_TYPES = ["train", "test", "validation"]
MODEL_TYPES = ["eukaryote", "eukaryotesmall" "human"]


rule all:
    input:
        "outputs/models/datasets/3_stats/set_summary.tsv",
        expand(
            "outputs/models/build/rnasamba/1_evaluation/{model_type}/accuracy_metrics_{dataset_type}.tsv",
            model_type=MODEL_TYPES,
            dataset_type=DATASET_TYPES,
        ),


rule download_ensembl_data:
    """
    Download ensembl cDNA and ncRNA files.
    Ensembl annotates protein coding and non-coding RNA transcripts in their files.
    This information will be used to separate protein coding from non-coding RNAs to build an RNAsamba model.
    Note this download renames genome files from their names on ensembl to make them simpler to point to.
    See example transformations below:
    - cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -> cdna/Homo_sapiens.GRCh38.cdna.fa.gz (dropped "all.")
    - ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz   -> ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz (no change)
    """
    output:
        "inputs/models/datasets/ensembl/{rna_type}/{genome}.{rna_type}.fa.gz",
    run:
        genome_df = metadata.loc[(metadata["genome"] == wildcards.genome)]
        root_url = genome_df["root_url"].values[0]
        if wildcards.rna_type == "cdna":
            suffix = genome_df["cdna_suffix"].values[0]
        else:
            suffix = genome_df["ncrna_suffix"].values[0]

        url = root_url + suffix
        shell("curl -JLo {output} {url}")


rule extract_protein_coding_orfs_from_cdna:
    """
    Ensembl cDNA files consist of transcript sequences for actual and possible genes, including pseudogenes, NMD and the like.
    Transcripts in the cDNA files have headers like: >TRANSCRIPT_ID SEQTYPE LOCATION GENE_ID GENE_BIOTYPE TRANSCRIPT_BIOTYPE, where the gene_biotype and transcript_biotype both contain information about whether the gene is coding or not.
    """
    input:
        "inputs/models/datasets/ensembl/cdna/{genome}.cdna.fa.gz",
    output:
        "outputs/models/datasets/0_coding/{genome}.fa.gz",
    conda:
        "envs/seqkit.yml"
    shell:
        """
        seqkit grep --use-regexp --by-name --pattern "transcript_biotype:protein_coding" -o {output} {input}
        """


rule download_validation_data:
    output:
        "inputs/models/datasets/validation/rnachallenge/{validation_type}.fa.gz",
    shell:
        """
        curl -JL https://raw.githubusercontent.com/cbl-nabi/RNAChallenge/main/RNAchallenge/{wildcards.validation_type}.fa | gzip > {output}
        """


rule combine_sequences:
    input:
        coding=expand("outputs/models/datasets/0_coding/{genome}.fa.gz", genome=GENOMES),
        noncoding=expand(
            "inputs/models/datasets/ensembl/ncrna/{genome}.ncrna.fa.gz", genome=GENOMES
        ),
        validation=expand(
            "inputs/models/datasets/validation/rnachallenge/{validation_type}.fa.gz",
            validation_type=VALIDATION_TYPES,
        ),
    output:
        "outputs/models/datasets/1_homology_reduction/all_sequences.fa",
    shell:
        """
        cat {input} | gunzip > {output}
        """


rule grab_all_sequence_names_and_lengths:
    input:
        "outputs/models/datasets/1_homology_reduction/all_sequences.fa",
    output:
        "outputs/models/datasets/1_homology_reduction/all_sequences.fa.seqkit.fai",
    conda:
        "envs/seqkit.yml"
    shell:
        """
        seqkit faidx -f {input}
        """


rule reduce_sequence_homology:
    """
    To reduce pollution between training and testing set, cluster sequences at 80% sequence identity.
    """
    input:
        "outputs/models/datasets/1_homology_reduction/all_sequences.fa",
    output:
        "outputs/models/datasets/1_homology_reduction/clustered_sequences_rep_seq.fasta",
        "outputs/models/datasets/1_homology_reduction/clustered_sequences_cluster.tsv",
    params:
        prefix="outputs/models/datasets/1_homology_reduction/clustered_sequences",
    conda:
        "envs/mmseqs2.yml"
    shell:
        """
        mmseqs easy-cluster {input} {params.prefix} tmp_mmseqs2 --min-seq-id 0.8 --cov-mode 1 --cluster-mode 2
        """


rule grab_validation_set_names_and_lengths:
    """
    The train/test data set sequences are identifiable by the genome information in the header, which is consistently formatted by Ensembl.
    The same is not true for the validation data.
    This rule grabs the validation sequence header names so they can be separated from the train/test sets.
    """
    input:
        "inputs/models/datasets/validation/rnachallenge/{validation_type}.fa.gz",
    output:
        validation="inputs/models/datasets/validation/rnachallenge/{validation_type}.fa",
        validation_fai="inputs/models/datasets/validation/rnachallenge/{validation_type}.fa.seqkit.fai",
    conda:
        "envs/seqkit.yml"
    shell:
        """
        cat {input} | gunzip > {output.validation}
        seqkit faidx -f {output.validation}
        """


rule process_sequences_into_nonoverlapping_sets:
    input:
        metadata="inputs/models/datasets/train_data_links.tsv",
        all_fai="outputs/models/datasets/1_homology_reduction/all_sequences.fa.seqkit.fai",
        validation_fai=expand(
            "inputs/models/datasets/validation/rnachallenge/{validation_type}.fa.seqkit.fai",
            validation_type=VALIDATION_TYPES,
        ),
        clusters="outputs/models/datasets/1_homology_reduction/clustered_sequences_cluster.tsv",
    output:
        expand(
            "outputs/models/datasets/2_sequence_sets/{coding_type}_{dataset_type}.txt",
            coding_type=CODING_TYPES,
            dataset_type=DATASET_TYPES,
        ),
    conda:
        "envs/tidyverse.yml"
    script:
        "scripts/process_sequences_into_nonoverlapping_sets.R"


rule filter_sequence_sets:
    input:
        fa="outputs/models/datasets/1_homology_reduction/clustered_sequences_rep_seq.fasta",
        names="outputs/models/datasets/2_sequence_sets/{coding_type}_{dataset_type}.txt",
    output:
        "outputs/models/datasets/2_sequence_sets/{coding_type}_{dataset_type}.fa",
    conda:
        "envs/seqtk.yml"
    shell:
        """
        seqtk subseq {input.fa} {input.names} > {output}
        """


##################################################################
## Build RNAsamba model
##################################################################


rule pip_install_rnasamba_no_deps:
    """
    To take advantage of nvidia GPU on AWS instance "Deep Learning Base OSS Nvidia Driver GPU AMI (Ubuntu 20.04) 20240122" (ami-07eb000b3340966b0),
    we need to install specific versions of tensorflow and other dependencies.
    This is accomplished in the envs/rnasamba.yml file, however rnasamba itself is not installed there because we need to use the command:
    pip install --no-deps rnasamba
    and there is no way to specify the "--no-deps" flag in a yaml file.
    This rule installs rnasamba into the conda-generated environment.
    """
    output:
        "outputs/models/build/rnasamba/rnasamba_installed.txt",
    conda:
        "envs/rnasamba.yml"
    shell:
        """
        pip install 'nvidia-tensorflow~=1.15'
        pip install --no-deps rnasamba # used version 0.2.5
        touch {output}
        """


rule build_rnasamba_model:
    """
    Build a new rnasamba model from the training data curated above.
    The --early_stopping parameter reduces training time and can help avoid overfitting.
    It is the number of epochs after lowest validation loss before stopping training.
    """
    input:
        rnasamba="outputs/models/build/rnasamba/rnasamba_installed.txt",
        fa=expand(
            "outputs/models/datasets/2_sequence_sets/{coding_type}_train.fa",
            coding_type=CODING_TYPES,
        ),
    output:
        "outputs/models/build/rnasamba/0_model/eukaryote_rnasamba.hdf5",
    conda:
        "envs/rnasamba.yml"
    benchmark:
        "benchmarks/models/build/rnasamba/0_model/eukaryote_rnasamba.tsv"
    shell:
        """
        rnasamba train --early_stopping 5 --verbose 2 {output} {input.fa[0]} {input.fa[1]}
        """


rule assess_rnasamba_model:
    input:
        model="outputs/models/build/rnasamba/0_model/{model_type}_rnasamba.hdf5",
        faa="outputs/models/datasets/2_sequence_sets/{coding_type}_{dataset_type}.fa",
    output:
        faa="outputs/models/build/rnasamba/1_evaluation/{model_type}/{coding_type}_{dataset_type}.fa",
        predictions="outputs/models/build/rnasamba/1_evaluation/{model_type}/{coding_type}_{dataset_type}.tsv",
    benchmark:
        "benchmarks/models/build/rnasamba/1_evaluation/{model_type}/{coding_type}_{dataset_type}.tsv"
    conda:
        "envs/rnasamba.yml"
    shell:
        """
        rnasamba classify --protein_fasta {output.faa} {output.predictions} {input.faa} {input.model}
        """


rule calculate_rnasamba_model_accuracy:
    input:
        expand(
            "outputs/models/build/rnasamba/1_evaluation/{{model_type}}/{coding_type}_{{dataset_type}}.tsv",
            coding_type=CODING_TYPES,
        ),
    output:
        freq="outputs/models/build/rnasamba/1_evaluation/{model_type}/confusionmatrix_{dataset_type}.tsv",
        metrics="outputs/models/build/rnasamba/1_evaluation/{model_type}/accuracy_metrics_{dataset_type}.tsv",
    conda:
        "envs/caret.yml"
    script:
        "scripts/calculate_rnasamba_model_accuracy.R"


rule download_rnasamba_human_model:
    """
    Use this model to compare whether the new model performs better or worse.
    It's saved under a new name so we can use a wildcard to run rnasamba classify and to calculate model accuracy.
    """
    output:
        "outputs/models/build/rnasamba/0_model/human_rnasamba.hdf5",
    shell:
        """
        curl -JLo {output} https://github.com/apcamargo/RNAsamba/raw/master/data/full_length_weights.hdf5
        """


##################################################################
## Get sequence statistics
##################################################################


rule get_sequence_descriptors:
    input:
        "outputs/models/datasets/2_sequence_sets/{coding_type}_{dataset_type}.fa",
    output:
        "outputs/models/datasets/2_sequence_sets/{coding_type}_{dataset_type}.fa.seqkit.fai",
    conda:
        "envs/seqkit.yml"
    shell:
        """
        seqkit faidx -f {input}
        """


rule calculate_sequence_statistics:
    input:
        expand(
            "outputs/models/datasets/2_sequence_sets/{coding_type}_{dataset_type}.fa.seqkit.fai",
            coding_type=CODING_TYPES,
            dataset_type=DATASET_TYPES,
        ),
    output:
        set_summary="outputs/models/datasets/3_stats/set_summary.tsv",
        set_length_summary="outputs/models/datasets/3_stats/set_length_summary.tsv",
        set_length_genome_summary="outputs/models/datasets/3_stats/set_length_genome_summary.tsv",
    conda:
        "envs/tidyverse.yml"
    script:
        "scripts/calculate_sequence_statistics.R"


##################################################################
## ugly repetitive code to make sORF model
##################################################################


rule filter_sequence_sets_to_small:
    """
    TER note this will probably create imbalanced classes.
    These names could be pulled by the process R scripts...
    """
    input:
        fa="outputs/models/datasets/2_sequence_sets/{coding_type}_{dataset_type}.fa",
    output:
        fa="outputs/models/datasets/2_sequence_sets/small/{coding_type}_{dataset_type}.fa",
    conda:
        "envs/seqtk.yml"
    shell:
        """
        seqtk seq --max-len 300 {input.fa}  > {output}
        """


rule build_rnasamba_model:
    """
    Build a new rnasamba model from the training data curated above.
    The --early_stopping parameter reduces training time and can help avoid overfitting.
    It is the number of epochs after lowest validation loss before stopping training.
    """
    input:
        rnasamba="outputs/models/build/rnasamba/rnasamba_installed.txt",
        fa=expand(
            "outputs/models/datasets/2_sequence_sets/small/{coding_type}_train.fa",
            coding_type=CODING_TYPES,
        ),
    output:
        "outputs/models/build/rnasamba/0_model/eukaryotesmall_rnasamba.hdf5",
    conda:
        "envs/rnasamba.yml"
    benchmark:
        "benchmarks/models/build/rnasamba/0_model/eukaryotesmall_rnasamba.tsv"
    shell:
        """
        rnasamba train --early_stopping 5 --verbose 2 {output} {input.fa[0]} {input.fa[1]}
        """
