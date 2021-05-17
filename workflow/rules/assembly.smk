rule create_samples_file_for_trinity:
    input: 
        left=expand("results/concatenated/{sample}.R1.fq.gz",sample=samples['sample']),
        right=expand("results/concatenated/{sample}.R2.fq.gz",sample=samples['sample']),
        samples=config["samples"]
    output: 
        samples_file="results/samples_file/samples.txt"
    run:
        import re
        samples = pd.read_csv(input.samples, sep=",").set_index(["sample"], drop=False)
        sample_old=''
        with open(output.samples_file, "w") as out:
            for file in input.left:
                sample = re.split('[/ | .]', file)[2]
                if sample==sample_old:
                    continue
                condition = samples.loc[(sample), ["condition"]].dropna()['condition'][0]
                left=file
                right=file.replace('.R1.','.R2.')
                out.write(f"{condition}\t{sample}\t{left}\t{right}\n")
                sample_old=sample



rule trinity_assembly:
    input: 
        samples_file="results/samples_file/samples.txt"
    output: 
        out_assembly="results/assembly_trinity/Trinity.fasta"
    conda:
        "../envs/assembly.yaml"
    threads: threads_trinity
    params:
        **config["params"]
    log:
        "results/logs/trinity_assembly/assembly.log"
    shell:
        """
        Trinity --seqType fq --samples_file {input.samples_file} --CPU {threads} {params.trinity}  --output $(dirname {output.out_assembly} ) 2> {log}
        """