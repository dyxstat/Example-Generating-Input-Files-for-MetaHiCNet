# Instruction to generating input files for MetaHiCNet 

Some scripts to process the intermediate data are available in the folder [Scripts](https://github.com/dyxstat/Example-Generating-Input-Files-for-MetaHiCNet/tree/main/Scripts).

## Process metagenomic Hi-C (metaHi-C) datasets
We use the sheep gut metaHi-C datasets as an example; the corresponding input files can be downloaded from the MetaHiCNet web server as the input format example.

**Version of softwares exploited in the analyses**
```
fastq_dump command from Sratoolkit: v2.10.8

bbduk.sh and clumpify.sh command from BBTools suite: v37.25

bwa command from BWA MEM: v0.7.17

samtools command from Samtools: v1.15.1

MetaCC.py command from MetaCC: v1.2.0

PPR_Meta command from PPR-Meta: v1.1

demovir.sh command from DemoVir: https://github.com/feargalr/Demovir

makeblastdb and blastn command from BLAST: v2.12.0
```

**Step 1: Preprocess the raw data**
```

fastq-dump --split-files --gzip SRR14350344.1

bbduk.sh  in1=SRR14350344.1_1.fastq.gz in2=SRR14350344.1_2.fastq.gz out1=HIC1_AQ.fastq.gz out2=HIC2_AQ.fastq.gz ref=/path_to_bbmap/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 minlen=50 tpe tbo

bbduk.sh  in1=HIC1_AQ.fastq.gz in2=HIC2_AQ.fastq.gz out1=HIC1_CL.fastq.gz out2=HIC2_CL.fastq.gz trimq=10 qtrim=r ftm=5 minlen=50

bbduk.sh in1=HIC1_CL.fastq.gz in2=HIC2_CL.fastq.gz out1=HIC1_trim.fastq.gz out2=HIC2_trim.fastq.gz ftl=10

clumpify.sh in1=HIC1_trim.fastq.gz in2=HIC2_trim.fastq.gz out1=HIC1_dedup.fastq.gz out2=HIC2_dedup.fastq.gz dedupe
```

**Step 2: Download assembled contigs and align processed Hi-C reads to contigs**
```
# Download assembled contigs (flye.v29.sheep_gut.hifi.250g.fasta) provided by original authors
wget https://zenodo.org/records/5228989/files/flye.v29.sheep_gut.hifi.250g.fasta.gz?download=1
gunzip -c flye.v29.sheep_gut.hifi.250g.fasta.gz > assembly.fa

bwa index assembly.fa
bwa mem -5SP assembly.fa HIC1.fastq.gz HIC2.fastq.gz > HIC_MAP.sam
samtools view -F 0x904 -bS HIC_MAP.sam > HIC_MAP_UNSORTED.bam
samtools sort -n HIC_MAP_UNSORTED.bam -o HIC_MAP_SORTED.bam
```

**Step 3: Generate the *Contig Information File and Raw Hi-C Contact File***
```
# Run MetaCC pipeline
python ./MetaCC.py norm -v assembly.fa HIC_MAP_SORTED.bam out_sheep_gut
python ./MetaCC.py bin --cover -v assembly.fa out_sheep_gut
```
After running MetaCC pipeline, the ***Contig Information File and Raw Hi-C Contact File*** can be obtained from 
```
* out_sheep_gut/contig_info.csv
* out_sheep_gut/Raw_contact_matrix.npz
```

**Step 4: Generate the *Binning Information File***

The MetaCC pipeline can generate contig binning results in a specified directory, where each FASTA file represents a single bin. For example:
```
out_sheep_gut/BIN
```

To convert these binning results into a standardized *Binning Information Fil*e, run the following script:
```
chmod +x generate_binning_csv.sh
./generate_binning_csv.sh out_sheep_gut/BIN fa
```
Users may also apply other contig binning tools and use the same *generate_binning_csv.sh* script to create a compatible *Binning Information File* from their results.

**Step 5: Identify viral and plasmid contigs**

We use PPR-Meta to detect both viral and plasmid contigs from the assembled metagenomic data. Users may employ alternative tools to identify viral and/or plasmid sequences. PPR-Meta assigns scores to each contig for three possible categories: phage, plasmid, and chromosome. The contig is classified based on the highest score, unless all scores are below the specified threshold, in which case it is labeled as uncertain.

Run PPR-Meta with the assembled contigs as input:
```
./PPR_Meta assembly.fa mge_results.csv
```
This command will generate a results file (mge_results.csv) containing classification scores for each contig.

To extract contigs with high confidence as phage or plasmid, use the following script. It filters contigs based on a user-defined cutoff (e.g., 0.7 for both phage and plasmid scores):
```
chmod +x extract_mge.sh
./extract_mge.sh mge_results.csv assembly.fa 0.7
```

This will generate two FASTA files:
```
* viral_contig.fa: Contigs with phage scores ≥ 0.7
* plasmid_contig.fa: Contigs with plasmid scores ≥ 0.7
```

These files can be used for downstream annotation of mobile genetic elements.

**Step 6: Generate the Taxonomy Information File**
To generate the Taxonomy Information File required by MetaHiCNet, we integrate taxonomic annotations from multiple tools applied to different types of contigs and bins.

1. Viral contigs: contigs identified as viral (in viral_contig.fa) are annotated using DemoVir (https://github.com/feargalr/Demovir) in this example, other tools such as VPF-class and geNomad can also be used here:
   ```
   ./demovir.sh viral_contig.fasta
   ```
   
2. Plasmid contigs: contigs identified as plasmid (in plasmid_contig.fa) are classified using BLAST against the NCBI RefSeq plasmid database. Follow the steps below to download and prepare the BLAST database:
   ```
   # Download RefSeq plasmid sequences (example link)
   wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/plasmid.*.genomic.fna.gz
   gunzip *.gz
  
   # Combine into a single FASTA file
   cat *.fna > refseq_plasmid.fa

   # Create a BLAST database
   makeblastdb -in refseq_plasmid.fa -dbtype nucl -out refseq_plasmid_db
   blastn -query plasmid_contig.fa -db refseq_plasmid_db -out plasmid_blast_output.tsv -outfmt 6  -perc_identity 95 -evalue 1e-5
   ```
   
3. Host genome bins: contig bins generated in Step 4 (under out_sheep_gut/BIN/) are classified using GTDB-Tk with the classify_wf workflow:
   ```
   gtdbtk classify_wf --genome_dir out_sheep_gut/BIN/ --out_dir gtdbtk_output
   ```
   
4. Compile results

   After obtaining annotations from DemoVir, BLAST, and GTDB-Tk, a custom script can be used to integrate all classifications into a single Taxonomy Information File. This file should include:
   ```
   * ID: The unique identifier, which can be either a bin ID (if contigs are grouped into bins) or a contig ID (if contigs are annotated individually). 
   * Category: Specifies whether the bin or contig represents a chromosome, virus, plasmid, or unclassified entity.
   * Taxonomic classification columns: Users may include hierarchical taxonomic ranks (e.g., phylum, class, order, family, genus, species) to provide structured biological context.
   ```



