rule representatation_assembly:
    input: 
        trinity="results/bowtie2_index/Trinity.fasta",
        fq1="results/concatenated/{sample}.R1.fq.gz",
        fq2="results/concatenated/{sample}.R2.fq.gz"
    output:
        "results/representatation_assembly/{sample}_align_stats.txt" 
    conda:
        "../envs/assembly.yaml"
    threads: threads
    shell:
        """
        bowtie2 -p {threads} -q --no-unal -k 20 -x {input.trinity} -1 {input.fq1} -2 {input.fq2} > /dev/null 2> {output}
        """

rule representatation_assembly_build:
    input: 
        trinity="results/assembly_trinity/Trinity.fasta",
    output:
        "results/bowtie2_index/Trinity.fasta" 
    conda:
        "../envs/assembly.yaml"
    threads: threads
    shell:
        """
        ln -s `pwd`/{input.trinity} {output}
        bowtie2-build --threads {threads} {input.trinity} {output}
        """

rule full_length_transcripts_blast:
    input: 
        trinity="results/assembly_trinity/Trinity.fasta"
    output:
        blast="results/full-length_transcript/blastx.outfmt6"
    conda:
        "../envs/blast.yaml"
    params:
        **config["params"]['blast']
    threads: threads_blast
    shell: 
        """
        mkdir -p $(dirname {params.dbname}.gz ) && touch {params.dbname}.gz
        wget -O {params.dbname}.gz {params.link}
        gunzip {params.dbname}.gz
        makeblastdb -in {params.dbname} -dbtype prot
        blastx -query {input.trinity} -db {params.dbname} -out {output.blast} \
        -evalue 1e-20 -num_threads {threads} -max_target_seqs 1 -outfmt 6
        """

rule full_length_transcripts:
    input: 
        trinity="results/assembly_trinity/Trinity.fasta",
        blast="results/full-length_transcript/blastx.outfmt6"
    output:
        "results/full-length_transcript/blastx.outfmt6.txt.w_pct_hit_length"
    conda:
        "../envs/assembly.yaml"
    params:
        **config["params"]['blast']
    threads: threads_blast
    shell: 
        """
        $TRINITY_HOME/util/analyze_blastPlus_topHit_coverage.pl {input.blast} {input.trinity} {params.dbname} {output}
        """

rule full_length_transcripts_group:
    input: 
        trinity="results/assembly_trinity/Trinity.fasta",
        blast="results/full-length_transcript/blastx.outfmt6"
    output:
        blast="results/full-length_transcript/blastx.outfmt6.group"
    conda:
        "../envs/assembly.yaml"
    params:
        **config["params"]['blast']
    shell: 
        """
        $TRINITY_HOME/util/misc/blast_outfmt6_group_segments.pl {input.blast} {input.trinity} {params.dbname} > {output.blast}
        $TRINITY_HOME/util/misc/blast_outfmt6_group_segments.tophit_coverage.pl {output.blast} 
        """

rule run_busco:
    input: 
        trinity="results/assembly_trinity/Trinity.fasta"
    output: 
        "results/busco_transcriptome/trinity_busco"
    conda:
        "../envs/busco.yaml"
    shell:
        """
        busco -m transcriptome -f -i {input.trinity} --out_path $(dirname {output}) -o $(basename {output}) --auto-lineage-euk
        """ 

rule run_quantification:
    input: 
        trinity="results/assembly_trinity/Trinity.fasta",
        samples_file="results/samples_file/samples.txt"
    output:
        "results/quantification/kallisto.isoform.counts.matrix"
    conda:
        "../envs/assembly.yaml"
    threads: threads_trinity
    shell: 
        """
        $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts {input.trinity} --seqType fq \
        --samples_file {input.samples_file} --est_method kallisto --output_dir $(dirname {output}) --kallisto_add_opts '--bias -t {threads}' --thread_count {threads} \
        --gene_trans_map {input.trinity}.gene_trans_map --prep_reference

        $TRINITY_HOME/util/abundance_estimates_to_matrix.pl --est_method kallisto --gene_trans_map {input.trinity}.gene_trans_map \
        --out_prefix $(dirname {output})/kallisto --name_sample_by_basedir $(find . -name "*abundance.tsv" -print0 | sort -z | xargs -r0)
        """