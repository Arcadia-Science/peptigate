import pandas as pd
import os
import re

metadata = pd.read_csv("inputs/models/rnasamba/build/train_data_links.tsv", sep="\t")
GENOMES = metadata["genome"].unique().tolist()

RNA_TYPES = ["cdna", "ncrna"]  # inherits names from ensembl

VALIDATION_TYPES = [
    "mRNAs",
    "ncRNAs",
]  # inherits names from https://github.com/cbl-nabi/RNAChallenge

CODING_TYPES = ["coding", "noncoding"]

DATASET_TYPES = ["train", "test", "validation"]

# KC: 'human' is the existing model used as a reference, and 'eu' is the model trained here
MODEL_TYPES = ["eu", "human"]


wildcard_constraints:
    genome="|".join(GENOMES),
    rna_type="|".join(RNA_TYPES),
    validation_type="|".join(VALIDATION_TYPES),
    coding_type="|".join(CODING_TYPES),
    dataset_type="|".join(DATASET_TYPES),
    model_type="|".join(MODEL_TYPES),


rule all:
    input:
        "outputs/models/rnasamba/build/5_stats/set_summary.tsv",
        expand(
            "outputs/models/rnasamba/build/4_evaluation/{model_type}/accuracy_metrics_{dataset_type}.tsv",
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
        "inputs/ensembl/{rna_type}/{genome}.{rna_type}.fa.gz",
    run:
        genome_df = metadata.loc[(metadata["genome"] == wildcards.genome)]
        root_url = genome_df["root_url"].values[0]
        if wildcards.rna_type == "cdna":
            suffix = genome_df["cdna_suffix"].values[0]
        else:
            suffix = genome_df["ncrna_suffix"].values[0]

        url = root_url + suffix
        shell("curl -JLo {output} {url}")


rule subsample_emsembl_data:
    """
    """
    input:
        rules.download_ensembl_data.output,
    output:
        "inputs/ensembl/{rna_type}/{genome}.{rna_type}.subsampled.fa.gz",
    conda:
        "envs/seqkit.yml"
    shell:
        """
        seqkit sample -n 1000 -o {output} {input}
        """

rule extract_protein_coding_orfs_from_cdna:
    """
    Ensembl cDNA files consist of transcript sequences for actual and possible genes, including pseudogenes, NMD and the like.
    Transcripts in the cDNA files have headers like: >TRANSCRIPT_ID SEQTYPE LOCATION GENE_ID GENE_BIOTYPE TRANSCRIPT_BIOTYPE, 
    where the gene_biotype and transcript_biotype both contain information about whether the gene is coding or not.
    """
    input:
        "inputs/ensembl/cdna/{genome}.cdna.subsampled.fa.gz",
    output:
        "outputs/models/rnasamba/build/0_coding/{genome}.fa.gz",
    conda:
        "envs/seqkit.yml"
    shell:
        """
    seqkit grep --use-regexp --by-name --pattern "transcript_biotype:protein_coding" -o {output} {input}
    """


rule download_validation_data:
    """
    KC: download datasets of coding and nonoding RNA from the RNA challenge.
    """
    output:
        "inputs/validation/rnachallenge/{validation_type}.fa.gz",
    shell:
        """
    curl -JL https://raw.githubusercontent.com/cbl-nabi/RNAChallenge/main/RNAchallenge/{wildcards.validation_type}.fa | gzip > {output}
    """

rule subsample_validation_data:
    input:
        "inputs/validation/rnachallenge/{validation_type}.fa.gz",
    output:
        "inputs/validation/rnachallenge/{validation_type}.subsampled.fa.gz",
    conda:
        "envs/seqkit.yml"
    shell:
        """
        seqkit sample -n 1000 -o {output} {input}
        """

rule combine_sequences:
    """
    KC: this rule combines all of the coding and noncoding sequences from ensembl
    and the validation sequences from the RNA challenge.
    """
    input:
        coding=expand("outputs/models/rnasamba/build/0_coding/{genome}.fa.gz", genome=GENOMES),
        noncoding=expand("inputs/ensembl/ncrna/{genome}.ncrna.subsampled.fa.gz", genome=GENOMES),
        validation=expand(
            "inputs/validation/rnachallenge/{validation_type}.subsampled.fa.gz",
            validation_type=VALIDATION_TYPES,
        ),
    output:
        "outputs/models/rnasamba/build/1_homology_reduction/all_sequences.fa.gz",
    shell:
        """
    cat {input} > {output}
    """


rule reduce_sequence_homology:
    """
    To reduce pollution between training and testing set, cluster sequences at 80% sequence identity.
    """
    input:
        "outputs/models/rnasamba/build/1_homology_reduction/all_sequences.fa.gz",
    output:
        "outputs/models/rnasamba/build/1_homology_reduction/clustered_sequences_rep_seq.fasta",
        "outputs/models/rnasamba/build/1_homology_reduction/clustered_sequences_cluster.tsv",
    params:
        prefix="outputs/models/rnasamba/build/1_homology_reduction/clustered_sequences",
    conda:
        "envs/mmseqs2.yml"
    shell:
        """
    mmseqs easy-cluster {input} {params.prefix} tmp_mmseqs2 --min-seq-id 0.8 --cov-mode 1 --cluster-mode 2
    """


rule grab_validation_set_names_and_lengths:
    input:
        "inputs/validation/rnachallenge/{validation_type}.subsampled.fa.gz",
    output:
        validation="inputs/validation/rnachallenge/{validation_type}.fa",
        validation_fai="inputs/validation/rnachallenge/{validation_type}.fa.fai",
    conda:
        "envs/seqkit.yml"
    shell:
        """
    cat {input} | gunzip > {output.validation}
    seqkit faidx {output.validation}
    """


rule grab_traintest_coding_names_and_lengths:
    input:
        expand("outputs/models/rnasamba/build/0_coding/{genome}.fa.gz", genome=GENOMES),
    output:
        coding="outputs/models/rnasamba/build/2_sequence_sets/traintest/all_coding.fa",
        coding_fai="outputs/models/rnasamba/build/2_sequence_sets/traintest/all_coding.fa.fai",
    conda:
        "envs/seqkit.yml"
    shell:
        """
    cat {input} | gunzip > {output.coding}
    seqkit faidx {output.coding}
    """


rule grab_traintest_noncoding_names_and_lengths:
    input:
        expand("inputs/ensembl/ncrna/{genome}.ncrna.subsampled.fa.gz", genome=GENOMES),
    output:
        noncoding="outputs/models/rnasamba/build/2_sequence_sets/traintest/all_noncoding.fa",
        noncoding_fai="outputs/models/rnasamba/build/2_sequence_sets/traintest/all_noncoding.fa.fai",
    conda:
        "envs/seqkit.yml"
    shell:
        """
    cat {input} | gunzip > {output.noncoding}
    seqkit faidx {output.noncoding}
    """


rule process_sequences_into_nonoverlapping_sets:
    """
    KC: the outputs of this rule are text files of sequence_ids of a subset
    of the sequences in the combined fasta files. 

    KC: the outputs are six text files of sequence_ids for all combinations
    of coding/noncoding and train/test/validation.
    """
    input:
        traintest_fai=expand(
            "outputs/models/rnasamba/build/2_sequence_sets/traintest/all_{coding_type}.fa.fai",
            coding_type=CODING_TYPES,
        ),
        validation_fai=expand(
            "inputs/validation/rnachallenge/{validation_type}.fa.fai",
            validation_type=VALIDATION_TYPES,
        ),
        clusters="outputs/models/rnasamba/build/1_homology_reduction/clustered_sequences_cluster.tsv",
    output:
        # KC: I think the order of outputs is coding train/test/val, noncoding train/test/val
        expand(
            "outputs/models/rnasamba/build/2_sequence_sets/{coding_type}_{dataset_type}.txt",
            coding_type=CODING_TYPES,
            dataset_type=DATASET_TYPES,
        ),
    conda:
        "envs/tidyverse.yml"
    script:
        "scripts/process_sequences_into_nonoverlapping_sets.R"


rule filter_sequence_sets:
    """
    KC: this rule uses the sequence_ids in the files from `process_sequences_into_nonoverlapping_sets`
    to generate fasta files containing the corresponding sequences from the clustered sequences.
    """
    input:
        fa="outputs/models/rnasamba/build/1_homology_reduction/clustered_sequences_rep_seq.fasta",
        names="outputs/models/rnasamba/build/2_sequence_sets/{coding_type}_{dataset_type}.txt",
    output:
        "outputs/models/rnasamba/build/2_sequence_sets/{coding_type}_{dataset_type}.fa",
    conda:
        "envs/seqtk.yml"
    shell:
        """
    seqtk subseq {input.fa} {input.names} > {output}
    """


##################################################################
## Build RNAsamba model
##################################################################


rule build_rnasamba_model:
    """
    Build a new rnasamba model from the training data curated above.
    The --early_stopping parameter reduces training time and can help avoiding overfitting.
    It is the number of epochs after lowest validation loss before stopping training.
    """
    input:
        expand(
            "outputs/models/rnasamba/build/2_sequence_sets/{coding_type}_train.fa",
            coding_type=CODING_TYPES,
        ),
    output:
        "outputs/models/rnasamba/build/3_model/eu_rnasamba.hdf5",
    conda:
        "envs/rnasamba.yml"
    shell:
        # KC: input[0] should be "coding", input[1] "noncoding"
        """
        rnasamba train --early_stopping 5 --verbose 2 {output} {input[0]} {input[1]}
        """


rule assess_rnasamba_model:
    """
    KC: this rule is run for each combination of human/eu models, coding/noncoding, and train/test/validation.

    KC: note that rnasamba also outputs a protein fasta file for the predicted ORFs
    via the `--protein_fasta` CLI option. 

    KC: why run this rule with coding/noncoding and on test and validation datasets? 
    Shouldn't it only be run on test (and without a split into coding/noncoding)?
    """
    input:
        model="outputs/models/rnasamba/build/3_model/{model_type}_rnasamba.hdf5",
        faa="outputs/models/rnasamba/build/2_sequence_sets/{coding_type}_{dataset_type}.fa",
    output:
        faa="outputs/models/rnasamba/build/4_evaluation/{model_type}/{coding_type}_{dataset_type}.fa",
        predictions="outputs/models/rnasamba/build/4_evaluation/{model_type}/{coding_type}_{dataset_type}.tsv",
    benchmark:
        "benchmarks/models/rnasamba/build/4_evaluation/{model_type}/{coding_type}_{dataset_type}.tsv"
    conda:
        "envs/rnasamba.yml"
    shell:
        """
    rnasamba classify --protein_fasta {output.faa} {output.predictions} {input.faa} {input.model}
    """


rule calculate_rnasamba_model_accuracy:
    input:
        expand(
            "outputs/models/rnasamba/build/4_evaluation/{{model_type}}/{coding_type}_{{dataset_type}}.tsv",
            coding_type=CODING_TYPES,
        ),
    output:
        freq="outputs/models/rnasamba/build/4_evaluation/{model_type}/confusionmatrix_{dataset_type}.tsv",
        metrics="outputs/models/rnasamba/build/4_evaluation/{model_type}/accuracy_metrics_{dataset_type}.tsv",
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
        "outputs/models/rnasamba/build/3_model/human_rnasamba.hdf5",
    shell:
        """
    curl -JLo {output} https://github.com/apcamargo/RNAsamba/raw/master/data/full_length_weights.hdf5
    """


##################################################################
## Get sequence statistics
##################################################################


rule get_sequence_descriptors:
    """
    KC: this rule just builds the index files for the fasta files from `filter_sequence_sets`
    which are used only in `calculate_sequence_statistics`.
    """
    input:
        "outputs/models/rnasamba/build/2_sequence_sets/{coding_type}_{dataset_type}.fa",
    output:
        "outputs/models/rnasamba/build/2_sequence_sets/{coding_type}_{dataset_type}.fa.seqkit.fai",
    conda:
        "envs/seqkit.yml"
    shell:
        """
    seqkit faidx -f {input}
    """


rule calculate_sequence_statistics:
    """
    This rule does not depend on the predictions from the RNASamba models. 
    """
    input:
        expand(
            "outputs/models/rnasamba/build/2_sequence_sets/{coding_type}_{dataset_type}.fa.seqkit.fai",
            coding_type=CODING_TYPES,
            dataset_type=DATASET_TYPES,
        ),
    output:
        set_summary="outputs/models/rnasamba/build/5_stats/set_summary.tsv",
        set_length_summary="outputs/models/rnasamba/build/5_stats/set_length_summary.tsv",
        set_length_genome_summary="outputs/models/rnasamba/build/5_stats/set_length_genome_summary.tsv",
    conda:
        "envs/tidyverse.yml"
    script:
        "scripts/calculate_sequence_statistics.R"
