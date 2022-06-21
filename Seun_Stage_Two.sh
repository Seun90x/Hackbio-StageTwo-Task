!usr/bin/bash

#downloading sample datasets for genomic analysis

echo -e "\n Downloading data... \n"

mkdir -p raw_data 
cd raw_data

wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz     

cd ..

#downloading reference sequence

echo -e "\nSLGFSK-N_231335\nSLGFSK-T_231336\n" > list.txt
echo -e "\n Downloading reference sequence... \n"

wget https://zenodo.org/record/2582555/files/hg19.chr5_12_17.fa.gz

#unzip reference
gunzip hg19.chr5_12_17.fa.gz

#installing fastqc

conda install -c bioconda fastqc multiqc --yes

echo -e "\n Data Preprocessing... \n"

mkdir -p Fastqc_Reports  #create directory for the fastqc output

#running fastqc on reads

for sample in `cat list.txt`
do
        fastqc raw_data/${sample}*.fastq.gz -o Fastqc_Reports
done

multiqc Fastqc_Reports -o Fastqc_Reports        

#installing trimmomatic for removing low-quality sequences

conda install -c bioconda trimmomatic --yes
multiqc Fastqc_Reports -o Fastqc_Reports        

#installing trimmomatic for removing low-quality sequences

conda install -c bioconda trimmomatic --yes

mkdir -p trimmed_reads
mkdir -p trimmed_reads/Fastqc_results

for sample in `cat list.txt`
do
       trimmomatic PE -threads 8 raw_data/${sample}_r1_chr5_12_17.fastq.gz raw_data/${sample}_r2_chr5_12_17.fastq.gz \
               trimmed_reads/${sample}_r1_paired.fq.gz trimmed_reads/${sample}_r1_unpaired.fq.gz \
               trimmed_reads/${sample}_r2_paired.fq.gz trimmed_reads/${sample}_r2_unpaired.fq.gz \
               ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:keepBothReads \
               LEADING:3 TRAILING:10 MINLEN:25

       fastqc  trimmed_reads/${sample}_r1_paired.fq.gz  trimmed_reads/${sample}_r2_paired.fq.gz \
                 -o trimmed_reads/Fastqc_results
done 

multiqc  trimmed_reads/Fastqc_results  -o trimmed_reads/Fastqc_results
'/home/micholufemi/Pictures/Screenshot from 2022-06-18 03-18-39.png' 

#conda install -y -c bioconda bwa
#conda install -c bioconda samtools
#conda install -c bioconda bamtools

#Indexing reference file        
bwa index hg19.chr5_12_17.fa 

mkdir Mapping
   
#Perform alignment
bwa mem -R '@RG\tID:231335\tSM:Normal' hg19.chr5_12_17.fa trimmed_reads/SLGFSK-N_231335_r1_paired.fq.gz \
      trimmed_reads/SLGFSK-N_231335_r2_paired.fq.gz > Mapping/SLGFSK-N_231335.sam

bwa mem -R '@RG\tID:231336\tSM:Tumor' hg19.chr5_12_17.fa trimmed_reads/SLGFSK-T_231336_r1_paired.fq.gz \
       trimmed_reads/SLGFSK-T_231336_r2_paired.fq.gz > Mapping/SLGFSK-T_231336.sam      
#converting SAM file to BAM file

for sample in `cat list.txt`
do
        Convert SAM to BAM and sort it 
        samtools view -@ 20 -S -b Mapping/${sample}.sam | samtools sort -@ 32 > Mapping/${sample}.sorted.bam

        Index BAM file
        samtools index Mapping/${sample}.sorted.bam
done

#filtering mapped reads

for sample in `cat list.txt`
do
        #Filter BAM files
        samtools view -q 1 -f 0x2 -F 0x8 -b Mapping/${sample}.sorted.bam > Mapping/${sample}.filtered1.bam
done

#Viewing results of mapped reads filtering
samtools flagstat Mapping/${sample}.filtered1.bam

#removing  duplicates

#using the command markdup
for sample in `cat list.txt`
do
        samtools collate -o Mapping/${sample}.namecollate.bam Mapping/${sample}.filtered1.bam
        samtools fixmate -m Mapping/${sample}.namecollate.bam Mapping/${sample}.fixmate.bam
        samtools sort -@ 32 -o Mapping/${sample}.positionsort.bam Mapping/${sample}.fixmate.bam
        samtools markdup -@32 -r Mapping/${sample}.positionsort.bam Mapping/${sample}.clean.bam
done


#aligning BAM files
for sample in `cat list.txt`
do      
        cat Mapping/${sample}.clean.bam  | bamleftalign -f hg19.chr5_12_17.fa -m 5 -c > Mapping/${sample}.leftAlign.bam

#-c - compressed, -m - max-iterations

#recalibrating BAM files

for sample in `cat list.txt`
do
        samtools calmd -@ 32 -b Mapping/${sample}.leftAlign.bam hg19.chr5_12_17.fa > Mapping/${sample}.recalibrate.bam
done

#refiltering BAM files

for sample in `cat list.txt`
do
        bamtools filter -in Mapping/${sample}.recalibrate.bam -mapQuality <=254 > Mapping/${sample}.refilter.bam
done

echo -e "\n Variant Calling and Classification \n"

#installing VarScan somatic tool for identifying variants from mapped samples

wget https://sourceforge.net/projects/varscan/files/VarScan.v2.3.9.jar          

mkdir Variants

for sample in `cat list.txt`
do
        samtools mpileup -f hg19.chr5_12_17.fa Mapping/${sample}.refilter.bam --min-MQ 1 --min-BQ 28 \
                > Variants/${sample}.pileup
done

#calling variants
java -jar VarScan.v2.3.9.jar somatic Variants/SLGFSK-N_231335.pileup \
        Variants/SLGFSK-T_231336.pileup Variants/SLGFSK \
        --normal-purity 1  --tumor-purity 0.5 --output-vcf 1 

#merging vcf files using bcf tools
bgzip Variants/SLGFSK.snp.vcf > Variants/SLGFSK.snp.vcf.gz
bgzip Variants/SLGFSK.indel.vcf > Variants/SLGFSK.indel.vcf.gz
tabix Variants/SLGFSK.snp.vcf.gz
tabix Variants/SLGFSK.indel.vcf.gz
bcftools merge Variants/SLGFSK.snp.vcf.gz Variants/SLGFSK.indel.vcf.gz > Variants/SLGFSK.vcf

echo -e '\n Variant Annotation \n'
#functional annotation using snpEff

#downloading SnpEff database
#download jar file
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip

# Unzip file
unzip snpEff_latest_core.zip

#download snpEff database
java -jar snpEff.jar download hg19

#annotate variants
java -Xmx8g -jar snpEff/snpEff.jar hg19 Variants/SLGFSK.vcf > Variants/SLGFSK.ann.vcf

#clinical annotation using gemini

wget https://raw.github.com/arq5x/gemini/master/gemini/scripts/gemini_install.py
#python gemini_install.py /usr/local /usr/local/share/gemini

gemini load -v Variants/SLGFSK.ann.vcf -t snpEff Annotation/gemini.db
