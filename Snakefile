# Define all chromosomes to download
CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]

# Define the output directory
OUTPUT_DIR = "data/genome/human/"

rule python_path:
    shell:
        "which python"


rule download_chromosome:
    output:
        f"{OUTPUT_DIR}{{chrom}}.fa"
    params:
        url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/{chrom}.fa.gz"
    shell:
        """
        mkdir -p {OUTPUT_DIR}
        wget {params.url} -O {output}.gz
        gunzip -c {output}.gz > {output}
        rm {output}.gz
        """

rule download_gtf:
    output:
        "data/Homo_sapiens.GRCh38.109.gtf"
    shell:
        """
        wget -O {output}.gz ftp://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz
        gunzip {output}.gz
        """
rule shorten_gtf:
    output:
        "data/{filename}.gtf.pkl"
    input:
        "data/{filename}.gtf"
    run:
        from gtfparse import read_gtf
        import polars as pl
        from pprint import pprint

        gtf =  read_gtf(input[0])
        protein_coding=gtf.filter(pl.col("gene_biotype") == "protein_coding")
        transcripts = protein_coding.filter(pl.col("feature") == "exon")
        isoforms={}
        general={}
        for row in transcripts.rows(named=True):
            gene_id = row['gene_id']
            if not gene_id in isoforms:
                isoforms[gene_id] = {}
                general[gene_id] = {
                    "seqname": row["seqname"], 
                    "strand": row['strand'],
                    "name": row['gene_name'],
                    "gene_source": row['gene_source'],
                    "gene_version": row['gene_version'],
                    "gene_biotype": row['gene_biotype'],
                    "transcripts":  []
                    }
            if not row['transcript_id'] in isoforms[gene_id]:
                isoforms[gene_id][row['transcript_id']] = {}
                general[gene_id]["transcripts"].append(row['transcript_id'])
            isoforms[gene_id][row['transcript_id']][int(row['exon_number'])] = (row['start'],row['end'])

        import pickle
        with open(output[0], 'wb') as f:
            pickle.dump((general, isoforms), f)




rule set_output_folders:
    input:
        expand(f"{OUTPUT_DIR}{{chrom}}.fa", chrom=CHROMOSOMES),
        gtf="data/Homo_sapiens.GRCh38.109.gtf.pkl"
    output:
        path="output/{gene}/isoforms.fa"
    threads: 1
    run:
        import pickle
        from pprint import pprint
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        with open(input.gtf, 'rb') as f:
            general, isoforms = pickle.load(f)
        seqname=general[wildcards.gene]["seqname"]
        strand=general[wildcards.gene]["strand"]
        name=general[wildcards.gene]["name"]
        print(f"Processing {name}, {seqname} on {strand}")
        print(f"opening {seqname}")
        fname = f"data/genome/human/chr{seqname}.fa"
        with open(fname, 'r') as f:
            for record in SeqIO.parse(fname, "fasta"):
                # Find the sequence with the specified ID
                if record.id == f"chr{seqname}":
                    # Extract the desired region (start to end)
                    # Note: Python uses 0-based indexing, so subtract 1 from start if using 1-based coordinates
                    seq = record.seq
        records = []
        for transcript in general[wildcards.gene]["transcripts"]:
            print("transcribing ", transcript)
            transctipt_seq=""
            for i in sorted(isoforms[wildcards.gene][transcript].keys()):
                start,end=isoforms[wildcards.gene][transcript][i]
                if strand == "+":
                    transctipt_seq += seq[start-1:end]
                else:
                    transctipt_seq += seq[start-1:end].reverse_complement()
            transctipt_seq=transctipt_seq.transcribe()
            start_codon_pos = transctipt_seq.find('AUG')
            
            if start_codon_pos == -1:
                return "Start codon not found in the sequence."
            
            # Extract the sequence starting from the first start codon
            transctipt_seq = transctipt_seq[start_codon_pos:]            
            protein_seq = transctipt_seq.translate(to_stop=True)
            record = SeqRecord(
                Seq(protein_seq),           # The sequence itself
                id=transcript,            # The sequence ID
                description=""           # Optional description
                )
            records.append(record)
        with open(output.path, "w") as handle:
            SeqIO.write(records, handle, "fasta")

rule prepare_isoform_for_alphafold:
    input:
        "output/{gene}/isoforms.fa"
    output:
        "output/{gene}/{transcript}/seq.fa"
    run:
        from Bio import SeqIO
        
        # Extract the specific transcript ID from wildcards
        transcript_id = wildcards.transcript
        
        # Read the input fasta file and search for the transcript
        with open(input[0], "r") as fasta_in:
            for record in SeqIO.parse(fasta_in, "fasta"):
                if record.id == transcript_id:
                    # Write the specific transcript to the output file
                    with open(output[0], "w") as fasta_out:
                        record.id = "isoform"
                        record.description = "isoform"
                        SeqIO.write(record , fasta_out, "fasta")
                    break
            else:
                raise ValueError(f"Transcript {transcript_id} not found in {input[0]}")
rule run_colabfold:
    input:
        "output/{gene}/{transcript}/seq.fa"
    output:
        "output/{gene}/{transcript}/isoform_plddt.png",
        "output/{gene}/{transcript}/isoform.done.txt",
    shell:
        """
        colabfold_batch --amber "{input[0]}"  "output/{wildcards.gene}/{wildcards.transcript}/" 
        """
rule make_data_json:
    input:
        "output/{gene}/{transcript}/isoform_plddt.png",
        "output/{gene}/{transcript}/isoform.done.txt",
        pikled_data="data/Homo_sapiens.GRCh38.109.gtf.pkl",
    output:
        "output/{gene}/{transcript}/data.json"
    run:
        import json
        from pprint import pprint
        import pickle
        with open(input.pikled_data, 'rb') as f: 
            general, isoforms = pickle.load(f)
        gene = wildcards.gene
        transcript = wildcards.transcript
        data = general[gene]
        data["transcript"] = transcript
        data['coordinates'] = isoforms[gene][transcript]
        with open(output[0], 'w') as f:
            json.dump(data, f)

rule run_all_isoforms:
    input:
        "data/Homo_sapiens.GRCh38.109.gtf.pkl"
    output:
        "output/{gene}/done.txt",

    run:
        import pickle
        gene = wildcards.gene
        with open(input[0], 'rb') as f: 
            general, isoforms = pickle.load(f)
        to_do=[]
        for transcript in general[wildcards.gene]["transcripts"]:
            to_do.append( f"output/{gene}/{transcript}/isoform.done.txt")
            to_do.append( f"output/{gene}/{transcript}/data.json")
        to_do = " ".join(sorted(to_do))
        shell(
            f"snakemake {to_do} -c {threads} --rerun-incomplete -p"
            )
        shell(
            f"touch output/{gene}/done.txt"
            )


rule autofillList:
    input:
        "data/Homo_sapiens.GRCh38.109.gtf.pkl"
        
    output:
        "autofill.js"
    run:
        import pickle
        from glob import glob
        with open(input[0], 'rb') as f: 
            general, isoforms = pickle.load(f)
        with open(output[0], 'w') as f:
            f.write(f"const autofill = [\n")
            for i in ([i.split("/")[-1] for i in glob("output/*")]):
                f.write(f"    '{i}',\n")
                name = general[i]['name']
                f.write(f"'{name}',\n")
            f.write("];\n")
            f.write("const geneDict = {\n")
            for i in ([i.split("/")[-1] for i in glob("output/*")]):
                name = general[i]['name']
                f.write(f"    '{name}':'{i}',\n")
                f.write(f"    '{i}':'{i}',\n")
            f.write("};\n")
rule makeJsonForGene:
    input:
        "data/Homo_sapiens.GRCh38.109.gtf.pkl",
        "output/{gene}/done.txt",
        gtf = "data/Homo_sapiens.GRCh38.109.gtf",
        
    output:
        "output/{gene}/data.js"
    run:
        import pickle
        from glob import glob
        import json
        from gtfparse import read_gtf
        import polars as pl

        # Load the general and isoforms data from the pickle file
        with open(input[0], 'rb') as f:
            general, isoforms = pickle.load(f)

        # Get the data for the specific gene
        data = general[wildcards.gene]
        transcripts = data["transcripts"]

        # Initialize dictionaries for storing PDB content and scores
        d = {}
        s = {}

        for transcript in transcripts:
            # Read and store the content of the PDB file
            folded_path = glob(f"output/{wildcards.gene}/{transcript}/isoform_relaxed_rank_001_alphafold2_ptm_model_*_seed_*.pdb")[0]
            with open(folded_path, 'r') as pdb_file:
                pdb_content = pdb_file.read()
            d[transcript] = pdb_content  # Store the PDB content directly

            # Store the path to the scores file (or its content if needed)
            scores_path = glob(f"output/{wildcards.gene}/{transcript}/isoform_scores_rank_001_alphafold2_ptm_model_*_seed_*.json")[0]
            s[transcript] = scores_path  # Store the path to the scores file

        # Create the JavaScript object with PDB content and scores
        json_object = "const relaxed_pdb = " + json.dumps(d, indent=4) 
        json_object += ";\nconst scores = " + json.dumps(s, indent=4)
        gtf =  read_gtf(input["gtf"])
        # Filter for transcripts of the specific gene
        gtf = gtf.filter(
            (pl.col("gene_biotype") == "protein_coding") & 
            (pl.col("feature")      == "transcript") & 
            (pl.col("gene_id")      == wildcards.gene)
        )
        transcript_dict = dict(zip(gtf['transcript_id'].to_list(), gtf['transcript_name'].to_list()))
        json_object += ";\nconst transcript_names = " + json.dumps(transcript_dict, indent=4)

        print(transcript_dict)
        # Write the JavaScript object to the output file
        with open(output[0], 'w') as f:
            f.write(json_object)

genes = [
    "ENSG00000141510", #TP53
    "ENSG00000254647", #insulin
    "ENSG00000170315", #Ubiquitin UBB
    "ENSG00000197061", #Histone H4, HIST1H4A
    "ENSG00000164825", #Beta-Defensin1 DEFB1
    "ENSG00000154620", #thymosin beta 4 Y-linked 
    "ENSG00000164128", #Neuropeptide Y
    "ENSG00000101200", #Vasopressin
    ]

rule all:
    input:
        expand("output/{gene}/data.js", gene=genes),
    shell:
        "rm autofill.js"
        "snakemake autofill.js -c1 --rerun-incomplete -p"