/opt/busco/BUSCO.py -i ../assemblies/shasta/Diaz26_shastaV2/Assembly.fasta -l /data/references/BUSCO/bacteria_odb9 -m geno -t ../tmp -o shastaV2
cd Diaz26_shastaV2/
ls
cat Assembly.fasta 
cd ..
shasta --command=assemble --memoryBacking 2M --threads 8 --input ../../data/Diaz26.ONT.fasta --assemblyDirectory Diaz26_shasta > shasta.log
cd ..
porechop -i ../../data/Diaz26.ONT.fastq -o porechopped/Diaz26.ONT.chopped.fastq --discard_middle > porechopped/Diaz26_porechop.log
cd canu+pilon/
mkdir porechopped
porechop -i ../../data/Diaz26.all.fastq -o porechopped/Diaz26.all.chopped.fastq --discard_middle > porechopped/Diaz26.all_porechop.log
porechop -i ../../data/Diaz26.ONT.fastq -o porechopped/Diaz26.all.chopped.fastq --discard_middle > porechopped/Diaz26.all_porechop.log
screen -r
screen 9
ls -l
screen -r 9663
screen -r 9663.pts-23.khaleesi
screen -r 9
y
screen -d -r
screen -d -r 9663.pts-23.khaleesi
cd /projects/classes/bioinformatics/Spring2020/vangarcia/
screen -r
cd /projects/classes/bioinformatics/Spring2020/vangarcia/
ls
cd assemblies/
cd shasta/
ls
cd ..
cd canu+pilon/
ls
screen -r
cd /projects/classes/bioinformatics/Spring2020/vangarcia/assemblies/canu+pilon/
ls
ls nanofilt/
cd nanofilt/
cat Diaz26.ONT.porechopped.nanofilted.fastq.gz 
gunzip Diaz26.all.porechopped.nanofilted.fastq.gz 
cat Diaz26.all.porechopped.nanofilted.fastq
cd ..
mkdir canu_assembly
zcat porechopped/Diaz26.ONT.chopped.fastq.gz | NanoFilt -q 9 -l 500 | gzip > nanofilt/Diaz26.ONT.porechopped.nanofilted.fastq.gz
cd porechopped/
ls
zcat porechop/Diaz26.ONT.chopped.fastq.gz | NanoFilt -q 9 -l 500 | gzip > nanofilt/Diaz26.ONT.porechopped.nanofilted.fastq.gz
cd ..
zcat porechop/Diaz26.ONT.chopped.fastq.gz | NanoFilt -q 9 -l 500 | gzip > nanofilt/Diaz26.ONT.porechopped.nanofilted.fastq.gz
cd porechop
ls
cd ..
ls
porechop -i ../../data/Diaz26.ONT.fastq -o porechopped/Diaz26.ONT.chopped.fastq --discard_middle > porechopped/Diaz26_porechop.log
cd /projects/clas
cd /projects/classes/bioinformatics/Spring2020/vangarcia/assemblies/
ls shasta/
screen -r
cd /projects/classes/bioinformatics/Spring2020/vangarcia/assemblies/
cd ..
ls
cd assemblies/
cd unicycler/
cd unicycler_conservative/
ls
/opt/busco/BUSCO.py assembly.fasta -l /data/references/BUSCO/bacteria_odb9 -m geno -t ../tmp -o unicyclerV1
/opt/busco/BUSCO.py -i assembly.fasta -l /data/references/BUSCO/bacteria_odb9 -m geno -t ../tmp -o unicyclerV1
cd ..
cd flye/
ls
cd out_Diaz26_ont/
ls
cd 00-assembly/
ls
cd ..
ls
/opt/busco/BUSCO.py -i assembly.fasta -l /data/references/BUSCO/bacteria_odb9 -o Diaz26_flye -m genome
cd tmp/
ls
cd ..
ls
cd ..
ls
cd ..
ls
cd tmp/
ls
cd ..
ls
cd assemblies/
cd flye/
ls
cd out_Diaz26_ont/
ls
cat flye.log 
ls
/opt/busco/BUSCO.py -i assembly.fasta -l /data/references/BUSCO/bacteria_odb9 -o Diaz26_flye -m genome
cd run_Diaz26_flye/
ls
cd ../../
cd ..
screen -r
screen -r -d 9663
exit
cd /projects/classes/bioinformatics/Spring2020/vangarcia/assemblies/
ls
cd canu+pilon/
ls
cd porechopped/
ls
cd ..
ls porechop
rm -rf porechopped
porechop -i ../../data/Diaz26.ONT.fastq -o porechop/Diaz26.ONT.chopped.fastq --discard_middle > porechop/Diaz26_porechop.log
cd /projects/classes/bioinformatics/Spring2020/vangarcia/assemblies/
ls
ls unicycler/
cd unicycler/
cd unicycler_conservative/
ls
cd ..
/opt/busco/BUSCO.py -i ../assemblies/unicycler/unicycler_conservative/assembly.fasta -l /data/references/BUSCO/bacteria_odb9 -m geno -t ../tmp -o unicyclerV1
/opt/busco/BUSCO.py -i /unicycler_conservative/assembly.fasta -l /data/references/BUSCO/bacteria_odb9 -m geno -t ../tmp -o unicyclerV1
cd unicycler_conservative/
ls
/projects/classes/bioinformatics/Spring2020/reference_genomes
cd /projects/classes/bioinformatics/Spring2020/reference_genomes
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/099/645/GCF_900099645.1_IMG-taxon_2671180067_annotated_assembly/GCF_900099645.1_IMG-taxon_2671180067_annotated_assembly_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/099/645/GCF_900099645.1_IMG-taxon_2671180067_annotated_assembly/GCF_900099645.1_IMG-taxon_2671180067_annotated_assembly_genomic.gff.gz
ls
gunzip GCF_900099645.1_IMG-taxon_2671180067_annotated_assembly_genomic.fna.gz 
gunzip GCF_900099645.1_IMG-taxon_2671180067_annotated_assembly_genomic.gff.gz 
ls
mv GCF_900099645.1_IMG-taxon_2671180067_annotated_assembly_genomic.fna > GCF_900099645.1_IMG-taxon_2671180067_annotated_assembly_genomic.Pseudomonas_peli.fna
/projects/classes/bioinformatics/Spring2020/vangarcia/assemblies/
cd ..
/projects/classes/bioinformatics/Spring2020/vangarcia/assemblies/
ls
cd vangarcia/assemblies/unicycler/unicycler_conservative/
ls
less Diaz26_unicycler_to_16S.txt 
screen -r
cd /projects/classes/bioinformatics/Spring2020/vangarcia/assemblies/shasta/
ls
cat shasta.log 
screen -r
cd /projects/classes/bioinformatics/Spring2020/vangarcia/assemblies/unicycler/unicycler_conservative/
ls
blastn -db /data/references/BLAST/16SMicrobial -query assembly.fasta -out Diaz11_unicyc_to_16S.txt -evalue 1e-5 -outfmt 0 -num_descriptions 5 -
num_threads 4
blastn -db /data/references/BLAST/16SMicrobial -query assembly.fasta -out Diaz26_unicyc_to_16S.txt -evalue 1e-5 -outfmt 0 -num_descriptions 20 -num_threads 4
ls
less Diaz26_unicyc_to_16S.txt 
cd ../../../../../../
cd ..
ls
cd ..
cd /projects/classes/bioinformatics/Spring2020/reference_genomes
ls
rm -rf GCF_900099645.1_IMG-taxon_2671180067_annotated_assembly_genomic.fna 
rm -rf GCF_900099645.1_IMG-taxon_2671180067_annotated_assembly_genomic.gff 
ls
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/040/435/GCF_011040435.1_ASM1104043v1/GCF_011040435.1_ASM1104043v1_genomic.fna.gz
ls
gunzip GCF_011040435.1_ASM1104043v1_genomic.fna.gz 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/040/435/GCF_011040435.1_ASM1104043v1/GCF_011040435.1_ASM1104043v1_genomic.gff.gz
gunzip GCF_011040435.1_ASM1104043v1_genomic.gff.gz 
mv GCF_011040435.1_ASM1104043v1_genomic.gff > GCF_011040435.1_ASM1104043v1_genomic.Pseudomonas_psychrophila.gff
mv GCF_011040435.1_ASM1104043v1_genomic.gff>GCF_011040435.1_ASM1104043v1_genomic.Pseudomonas_psychrophila.gff
mv GCF_011040435.1_ASM1104043v1_genomic.gff > GCF_011040435.1_ASM1104043v1_genomic.Pseudomonas_psychrophila.gff
pwd
cd ..
cd vangarcia/
cd assemblies/
ls
cd unicycler/
ls
screen -r
screen -r -d 9663
cd /projects/classes/bioinformatics/Spring2020/vangarcia/assemblies/unicycler/
ls
screen -r
cd /proje
cd /projects/classes/bioinformatics/Spring2020/vangarcia/assemblies/
ls
cd 
cd /projects/classes/bioinformatics/Spring2020/vangarcia/assemblies/
cd unicycler/
ls
ls unicycler_bold/
cd unicycler_conservative/
ls
cd ../unicycler_normal/
screen -r 
screen -r
cd /projects/classes/bioinformatics/Spring2020/vangarcia/assemblies/
cd unicycler/
cd unicycler_conservative/
ls
cd unicycler.log
less unicycler.log
ls
cd run_unicyclerV1/
ls
less short_summary_unicyclerV1.txt 
cd ../../
cd unicycler_normal/
ls
cd run_unicyclerV2/
ls
less short_summary_unicyclerV2.txt 
cd ..
cd unicycler_bold/
ls
screen -r
/opt/busco/BUSCO.py -i assembly.fasta -l /data/references/BUSCO/bacteria_odb9 -m geno -t ../tmp -o unicyclerV3
export PATH=$PATH:/opt/SPAdes-3.14.0-Linux/bin
echo $PATH
/opt/busco/BUSCO.py -i assembly.fasta -l /data/references/BUSCO/bacteria_odb9 -m geno -t ../tmp -o unicyclerV3
cd ..
export PATH=$PATH:/opt/SPAdes-3.14.0-Linux/bin
echo $PATH
/opt/busco/BUSCO.py -i assembly.fasta -l /data/references/BUSCO/bacteria_odb9 -m geno -t ../tmp -o unicyclerV3
 cd unicycler_bold/
ls
/opt/busco/BUSCO.py -i assembly.fasta -l /data/references/BUSCO/bacteria_odb9 -m geno -t ../tmp -o unicyclerV3
export PATH=$PATH:/opt/SPAdes-3.14.0-Linux/bin
echo $PATH
/opt/busco/BUSCO.py -i assembly.fasta -l /data/references/BUSCO/bacteria_odb9 -m geno -t ../tmp -o unicyclerV3
cd /projects/classes/bioinformatics/Spring2020/vangarcia/assemblies/
ls
cd unicycler/
ls
cd unicycler_conservative/
ls
cd r
cd run_unicyclerV1/
ls
cd ..
ls
cd Diaz26_flye_report/
ls
cd /projects/classes/bioinformatics/Spring2020/vangarcia/assemblies/
cd unicycler/
cd unicycler_conservative/Diaz26_flye_report/
ls
less report.txt
cd ..
cd unicycler_normal/
cd Diaz26_flye_report/
ls
less report.txt
cd ../..
cd unicycler_bold/
cd Diaz26_flye_report/
ls
less report.txt
cd ../../../../..
cd vangarcia/
cd assemblies/canu+pilon/
screen -r
screen -r
screen -r -d 9663
screen -r -d 9663
screen -r
screen -r
screen -r -d 9663
cd /projects/classes/bioinformatics/Spring2020/
ls
cd diazinon_genomes/
ls
cd Diaz26
ls
cd ragoo
ls
cd ..
cd unicycler/
ls
cd ..
cd genemark/
ls
cd ..
cd /projects/classes/bioinformatics/Spring2020/
ls
cd diazinon_genomes/
ls
cd Diaz26/
ls
cd ragoo
ls
cd ..
cd ragoo_output/
ls
cd ..
cd unicycler/
ls
cd read_alignment/
ls
cd ../../
ls
cd data/
ls
cd ..
cd data/
ls
/opt/minimap2/minimap2 -a -x map-ont -o final_assembly.fasta ../../../../data/Diaz26.ONT.greater7000bp_1.fasta
cd ..
cd unicycler/
ls
cd ..
ls
cd flye/
ls
cd ..
cd BUSCO/
ls
cd .
cd ..
cd genemark/
ls
cd ..
ls
cd ..
cd vangarcia/
ls
cd data/
ls
cd ../../
cd diazinon_genomes/
cd Diaz26/
ls
cd data/
ls
/opt/minimap2/minimap2 -a -x map-ont -o final_assembly.fasta ../../../vangarcia/data/Diaz26.ONT.greater7000bp_1.fasta
/opt/minimap2/minimap2 -a -x map-ont -o final_assembly.fasta ragoo.fasta ../../../vangarcia/data/Diaz26.ONT.greater7000bp_1.fasta
minimap2 -a -x map-ont -o ragoo.mapped.sam data/final_assembly.fasta data Diaz26.all.fasta && samtools view -@ 20 -Sb ragoo.mapped.sam > ragoo.mapped.bam && samtools sort -@ 20 ragoo.mapped.bam -o ragoo.sorted.bam && samtools index ragoo.sorted.bam 
bwa index data/final_assembly.fasta && bwa mem data/final_assembly.fasta  data/Diaz-26_R1.fq.gz_paired.fq.gz data/Diaz-26_R2.fq.gz_paired.fq.gz  >  PE_maped_to_final.sam && samtools view -Sb PE_maped_to_final.sam >  PE_maped_to_final.bam && samtools sort PE_maped_to_final.bam -o PE_maped_to_final.sorted.bam && samtools index PE_maped_to_final.sorted.bam &
cd ..
minimap2 -a -x map-ont -o ragoo.mapped.sam data/final_assembly.fasta data Diaz26.all.fasta && samtools view -@ 20 -Sb ragoo.mapped.sam > ragoo.mapped.bam && samtools sort -@ 20 ragoo.mapped.bam -o ragoo.sorted.bam && samtools index ragoo.sorted.bam 
bwa index data/final_assembly.fasta && bwa mem data/final_assembly.fasta  data/Diaz-26_R1.fq.gz_paired.fq.gz data/Diaz-26_R2.fq.gz_paired.fq.gz  >  PE_maped_to_final.sam && samtools view -Sb PE_maped_to_final.sam >  PE_maped_to_final.bam && samtools sort PE_maped_to_final.bam -o PE_maped_to_final.sorted.bam && samtools index PE_maped_to_final.sorted.bam &
minimap2 -a -x map-ont -o ragoo.mapped.sam data/final_assembly.fasta $nanopore && samtools view -@ 20 -Sb ragoo.mapped.sam > ragoo.mapped.bam && samtools sort -@ 20 ragoo.mapped.bam -o ragoo.sorted.bam && samtools index ragoo.sorted.bam
cd data/
ls
cd ..
cd ragoo_output/
ls
cd ..
ls
cd ragoo
ls
cd ..
ls
cd data/
ls
cd ..
minimap2 -a -x map-ont -o ragoo.mapped.sam data/final_assembly.fasta $nanopore && samtools view -@ 20 -Sb ragoo.mapped.sam > ragoo.mapped.bam && samtools sort -@ 20 ragoo.mapped.bam -o ragoo.sorted.bam && samtools index ragoo.sorted.bam
/opt/minimap2/minimap2
/opt/minimap2/minimap2 -a -x map-ont -o ragoo.mapped.sam data/final_assembly.fasta data Diaz26.all.fasta && samtools view -@ 20 -Sb ragoo.mapped.sam > ragoo.mapped.bam && samtools sort -@ 20 ragoo.mapped.bam -o ragoo.sorted.bam && samtools index ragoo.sorted.bam 
bwa index data/final_assembly.fasta && bwa mem data/final_assembly.fasta  data/Diaz-26_R1.fq.gz_paired.fq.gz data/Diaz-26_R2.fq.gz_paired.fq.gz  >  PE_maped_to_final.sam && samtools view -Sb PE_maped_to_final.sam >  PE_maped_to_final.bam && samtools sort PE_maped_to_final.bam -o PE_maped_to_final.sorted.bam && samtools index PE_maped_to_final.sorted.bam &
ls
cd ragoo_output/
ls
cd ..
cd ragoo
ls
cd ..
cd data/
ls
cd ..
/opt/minimap2/minimap2 -a -x map-ont -o ragoo.mapped.sam data/final_assembly.fasta data Diaz26.all.fasta && samtools view -@ 20 -Sb ragoo.mapped.sam > ragoo.mapped.bam && samtools sort -@ 20 ragoo.mapped.bam -o ragoo.sorted.bam && samtools index ragoo.sorted.bam && bwa index data/final_assembly.fasta && bwa mem data/final_assembly.fasta  data/Diaz-26_R1.fq.gz_paired.fq.gz data/Diaz-26_R2.fq.gz_paired.fq.gz  >  PE_maped_to_final.sam && samtools view -Sb PE_maped_to_final.sam >  PE_maped_to_final.bam && samtools sort PE_maped_to_final.bam -o PE_maped_to_final.sorted.bam && samtools index PE_maped_to_final.sorted.bam & minimap2 -a -x map-ont -o ragoo.mapped.sam data/final_assembly.fasta data  Diaz26.all.fasta && samtools view -@ 20 -Sb ragoo.mapped.sam 
ragoo.mapped.bam && samtools sort -@ 20 ragoo.mapped.bam -o ragoo.sorted.bam && samtools index ragoo.sorted.bam 
cd data/
ls
cd ..
ls
/opt/minimap2/minimap2 -a -x map-ont -o ragoo.mapped.sam data/final_assembly.fasta data Diaz26.all.fasta && samtools view -@ 20 -Sb ragoo.mapped.sam > ragoo.mapped.bam && samtools sort -@ 20 ragoo.mapped.bam -o ragoo.sorted.bam && samtools index ragoo.sorted.bam && bwa index data/final_assembly.fasta && bwa mem data/final_assembly.fasta  data/Diaz-26_R1.fq.gz_paired.fq.gz data/Diaz-26_R2.fq.gz_paired.fq.gz  >  PE_maped_to_final.sam && samtools view -Sb PE_maped_to_final.sam >  PE_maped_to_final.bam && samtools sort PE_maped_to_final.bam -o PE_maped_to_final.sorted.bam && samtools index PE_maped_to_final.sorted.bam & minimap2 -a -x map-ont -o ragoo.mapped.sam data/final_assembly.fasta data  Diaz26.all.fasta && samtools view -@ 20 -Sb ragoo.mapped.sam > ragoo.mapped.bam && samtools sort -@ 20 ragoo.mapped.bam -o ragoo.sorted.bam && samtools index ragoo.sorted.bam 
ls
bash make_genome.sh
rm PE_reads_onto_flye.sam
ls
rm ragoo.mapped.sam 
cd /projects/classes/bioinformatics/Spring2020/
ls
cd diazinon_genomes/
cd Diaz26
ls
cd unicycler/
ls
cd ..
cd flye/
ls
cd ..
/opt/minimap2/minimap2 -a -x map-ont -o ragoo.mapped.sam data/final_assembly.fasta data Diaz26.all.fasta && samtools view -@ 20 -Sb ragoo.mapped.sam > ragoo.mapped.bam && samtools sort -@ 20 ragoo.mapped.bam -o ragoo.sorted.bam && samtools index ragoo.sorted.bam && bwa index data/final_assembly.fasta && bwa mem data/final_assembly.fasta  data/Diaz-26_R1.fq.gz_paired.fq.gz data/Diaz-26_R2.fq.gz_paired.fq.gz  >  PE_maped_to_final.sam && samtools view -Sb PE_maped_to_final.sam >  PE_maped_to_final.bam && samtools sort PE_maped_to_final.bam -o PE_maped_to_final.sorted.bam && samtools index PE_maped_to_final.sorted.bam & minimap2 -a -x map-ont -o ragoo.mapped.sam data/final_assembly.fasta data  Diaz26.all.fasta && samtools view -@ 20 -Sb ragoo.mapped.sam > ragoo.mapped.bam && samtools sort -@ 20 ragoo.mapped.bam -o ragoo.sorted.bam && samtools index ragoo.sorted.bam 
ls
ls ragoo_output/
ls data/
ls ragoo
ls
cd /projects/classes/bioinformatics/Spring2020/
ls
cd diazinon_genomes/
ls
cd Diaz26
ls
ls genemark/
/opt/minimap2/minimap2 -a -x map-ont -o final_assembly.fasta ragoo.fasta ../../../vangarcia/data/Diaz26.ONT.greater7000bp_1.fasta
cd unicycler/
ls
cd ..
cd data/
ls
cd ..
cd ragoo_output
ls
/opt/minimap2/minimap2 -a -x map-ont -o final_assembly.fasta ragoo.fasta ../../../vangarcia/data/Diaz26.ONT.greater7000bp_1.fasta
ls
samtools view -Sb ragoo.mapped.sam > ragoo.mapped.bam
cd ..
ls
/opt/minimap2/minimap2 -a -x map-ont -o final_assembly.fasta ragoo_output/ragoo.fasta ../../../vangarcia/data/Diaz26.ONT.greater7000bp_1.fasta
ls
samtools view -Sb ragoo.mapped.sam > ragoo.mapped.bam
ls
rm ragoo.mapped.sam
ls
samtools sort -o ragoo.mapped.sorted ragoo.mapped.bam
ls
samtools sort -o ragoo.mapped.sorted.bam ragoo.mapped.bam
samtools index ragoo.mapped.sorted.bam
ls ragoo_output
bwa mem -t 4 ragoo_output/ragoo.fasta ../../../../data/Diaz-26_F_paired.fq.gz ../../../../data/Diaz-26_R_paired.fq.gz | samtools view -Sb > ragoo.MiSeq.mapped.bam
bwa mem -t 4 ragoo_output/ragoo.fasta ../../../data/Diaz-26_F_paired.fq.gz ../../../data/Diaz-26_R_paired.fq.gz | samtools view -Sb > ragoo.MiSeq.mapped.bam
cd ../../../
cd Spring2020/
ls
cd vangarcia/data/
ls
cd ../../
cd diazinon_genomes/
cd Diaz26/
bwa mem -t 4 ragoo_output/ragoo.fasta ../../vangarcia/data/Diaz-26_F_paired.fq.gz ../../vangarcia/data/Diaz-26_R_paired.fq.gz | samtools view -Sb > ragoo.MiSeq.mapped.bam
ls
bwa index ragoo_output/ragoo.fasta
ls
samtools sort -o ragoo.MiSeq.mapped.sorted.bam ragoo.MiSeq.mapped.bam
samtools index ragoo.MiSeq.mapped.sorted.bam
samtools depth -a ragoo.mapped.sorted.bam | awk '{if($3<=20) print $0}'|less
ls
samtools depth -a ragoo.mapped.sorted.bam | awk '{if($3<=20) print $0}'|less
/opt/gms2_linux_64/gms2.pl --seq ragoo_output/ragoo.fasta --genome-type bacteria --output Diaz-26.ragoo.gff --format gff --fnn Diaz-26.ragoo.fnn --faa Diaz-26.ragoo.faa
ls
head Diaz-26.ragoo.faa
ls
rm PE_maped_to_final.sam
cd /projects/classes/bioinformatics/Spring2020/diazinon_genomes/
cd Diaz26/
mkdir mapping
mv ragoo/ragoo_output/ragoo.mapped.sorted.bam mapping/
mv ragoo_output/ragoo.mapped.sorted.bam mapping/
ls
mv ragoo.mapped.sorted.bam mapping/
History | grep minimap2
history | grep minimap2
cd /projects/classes/bioinformatics/Spring2020/diazinon_genomes/Diaz26/
ls
samtools depth -a ragoo.sorted.bam | tee  depth_out | awk '{sum+=$3; counter+=1} END{print sum/counter}' 
341.495/5
cat depth_out |  awk '{if ($3< 68){ print $1,$2,$2,$3;};}' | tr " " "\t" > less_then_68.bed 
cat less_then_(68.bed | bedtools merge  | awk '{if($3-$2 > 50){print $0}}'
cat less_then_(68.bed | bedtools merge  | awk '{if($3-$2 > 50){print $0}}')
cat less_then_68.bed | bedtools merge  | awk '{if($3-$2 > 50){print $0}}'
cd /projects/classes/bioinformatics/Spring2020/
ls
cd vangarcia/
ls
cd ..
ls
cd diazinon_genomes/
ls
cd Diaz26/
ls
cd flye/
ls
chmod 775
chmod 775*
chmod+ * 775
chmod 775 scaffolds.fasta
cd ..
ls
cd flye/
ls
cd ..
cd quast/
ls
cd ..
cd unicycler/
ls
chmod + 775
cd ../../
ls
cd Diaz26/
screen -r
cd /projects/classes/bioinformatics/Spring2020/vangarcia/assemblies/
ls
cd flye/
ls
less unicyc_to_flye.delta 
screen -r
cd /projects/classes/bioinformatics/Spring2020/diazinon_genomes/Diaz26/
ls
cd ../../
ls
cd vangarcia/
ls
cd assemblies/
ls
cd flye/
ls
ls out_Diaz26_ont/
/opt/MUMmer3.23/nucmer --mumreference -c 100 -p unicyc_to_flye ../../assemblies/flye/out_Diaz26_ont/assembly.fasta ../../assemblies/unicycler//unicycler_conservative/assembly.fasta 
pwd
cd /
cd ..
cd /projects/classes/bioinformatics/Spring2020/vangarcia/
cd assemblies/
cd shasta/
shasta --command=assemble --memoryBacking 2M --threads 8 --input ../../data/Diaz26.ONT.fasta --assemblyDirectory Diaz26_shasta > shasta.log
shasta --command=assemble --memoryBacking 2M --threads 8 --input ../../data/Diaz26.ONT.fasta --assemblyDirectory Diaz26_shasta
ls
rm -rf Diaz26_shasta/
c
cd shas
cd shasta
cat shasta.log
ls
rm -rf shasta.log
shasta --command=assemble --memoryBacking 2M --threads 8 --input ../../data/Diaz26.ONT.fasta --assemblyDirectory Diaz26_shasta > shasta.log
ls
cat shasta.log
cd Diaz26_shasta/
ls
cat Assembly.fasta 
cat AssemblySummary.
cat AssemblySummary.csv 
cd ..
cd canu+pilon/
porechop -i ../../data/Diaz26.ONT.fastq -o porechopped/Diaz26.ONT.chopped.fastq --discard_middle > porechopped/Diaz26_porechop.log
cat porechop/Diaz26.ONT.chopped.fastq.gz | NanoFilt -q 9 -l 500 | gzip > nanofilt/Diaz26.ONT.porechopped.nanofilted.fastq.gz
NanoFilt -h
mkdir nanofilt
zcat porechop/Diaz26.all.chopped.fastq.gz | NanoFilt -q 9 -l 500 | gzip > nanofilt/Diaz26.all.porechopped.nanofilted.fastq.gz 
ls
cd porechop
ls
cd ..
ls
cd porechopped
ls
cd ..
zcat porechopped/Diaz26.all.chopped.fastq.gz | NanoFilt -q 9 -l 500 | gzip > nanofilt/Diaz26.all.porechopped.nanofilted.fastq.gz 
cd porechopped/
ls
cd ..
cd porechop
ls
cd ..
zcat porechopped/Diaz26.all.chopped.fastq | NanoFilt -q 9 -l 500 | gzip > nanofilt/Diaz26.all.porechopped.nanofilted.fastq.gz 
zcat porechopped/Diaz26.all.chopped.fastq | NanoFilt -q 9 -l 500 | gzip > nanofilt/Diaz26.all.porechopped.nanofilted.fastq
zcat porechopped/Diaz26.ONT.chopped.fastq.gz | NanoFilt -q 9 -l 500 | gzip > nanofilt/Diaz26.ONT.porechopped.nanofilted.fastq.gz
zcat porechopped/Diaz26.ONT.chopped.fastq | NanoFilt -q 9 -l 500 | gzip > nanofilt/Diaz26.ONT.porechopped.nanofilted.fastq.gz
cd porechop
ls
cd Diaz26.ONT.porechop.log
cat Diaz26.ONT.porechop.log 
cd ..
cd porechopped
ls
cd ..
NanoFilt -h
ls
cd nanofilt/
ls
zcat porechop/Diaz26.all.chopped.fastq.gz | NanoFilt -q 9 -l 500 | gzip > nanofilt/Diaz26.all.porechopped.nanofilted.fastq.gz
cd ..
zcat porechop/Diaz26.all.chopped.fastq.gz | NanoFilt -q 9 -l 500 | gzip > nanofilt/Diaz26.all.porechopped.nanofilted.fastq.gz
cd nanofilt/
ls
cd ..
/opt/canu/Linux-amd64/bin/canu -p Diaz26_canu_nanofilt_porechopped -d canu_assembly -genomeSize=6m -nanopore-raw nanofilt/Diaz26.ONT.porechopped.nanofilted.fastq.gz useGrid=false > canu_assembly/Diaz26.canu.log
cd ..
cd shasta
shasta --command=assemble --memoryBacking 2M --threads 8 --input ../../data/Diaz26.ONT.fasta --assemblyDirectory Diaz26_shasta > shasta.log
ls
cd ..
ls
cd flye/
ls
cd ..
ls
cd canu+pilon/
ls
porechop -i ../../data/Diaz26.ONT.fastq -o porechop/Diaz26.ONT.chopped.fastq --discard_middle > porechop/Diaz26_porechop.log
fg
screen -r
ctrl -r -d 9663
screen -r -r 9663
screen -r -d 9663
fg
screen -r
cd ..
cd unicycler/
ls
cd unicycler_conservative/
cd ..
ls unicycler.log 
cat unicycler.log
ls
cd unicycler_conservative/
ls
cd run_unicyclerV1/
ls
cd ..
ls
cd ../../../../
cd ../../
pwd
cp /opt/genemark_suite_linux_64/gmsuite/gm_key ~/.gm_key
cp /opt/gms2_linux_64/gmhmmp2_key ~/.gmhmmp2_key
cd ..
cd /projects/classes/bioinformatics/Spring2020/vangarcia/assemblies/unicycler/unicycler_conservative/
ls
blastn -h
blastn -db /data/references/
blastn -db /data/references/BLAST/16SMicrobial -query assembly.fasta -out Diaz26_unicycler_to_16S.txt -evalue 1e-5 -outfmt 0 -num_descriptions 5 -num_threads 4
ls
cat Diaz26_unicycler_to_16S.txt 
less Diaz26_unicycler_to_16S.txt 
cd ..
unicycler -1 ../../data/Diaz-26_F_paired.fq.gz -2 ../../data/Diaz-26_R_paired.fq.gz -s ../../data/Diaz-11_unpaired.fq.gz -l ../../data/Diaz26.ONT.fastq -o unicycler_normal --vcf --mode normal --spades_tmp_dir /projects/temp > unicycler.log
unicycler -1 ../../data/Diaz-26_F_paired.fq.gz -2 ../../data/Diaz-26_R_paired.fq.gz -s ../../data/Diaz-26_unpaired.fq.gz -l ../../data/Diaz26.ONT.fastq -o unicycler_normal --vcf --mode normal --spades_tmp_dir /projects/temp > unicycler.log
export PATH=$PATH:/opt/SPAdes-3.14.0-Linux/bin
echo $PATH
unicycler -1 ../../data/Diaz-26_F_paired.fq.gz -2 ../../data/Diaz-26_R_paired.fq.gz -s ../../data/Diaz-26_unpaired.fq.gz -l ../../data/Diaz26.ONT.fastq -o unicycler_normal --vcf --mode normal --spades_tmp_dir /projects/temp > unicycler.log
python /opt/quast-4.6.3/quast.py -F_0000012.15.4_Release_6_plus_ISO1_MT_genomic.fna
python /opt/quast-4.6.3/quast.py -o quast_reports -t 4 BUSCO/GCF_0000012.15.4_Release_6_plus_ISO1_MT_genomic.fna
ls
python /opt/quast-4.6.3/quast.py -o _000151625.1_ASM15162v1_genomic.fna
cd q
cd quast_reports
ls
cd ..
pwd
cd ..
pwd
ls
cd quast_report/
ls
cd ..
cd vangarcia/
cd BUSCO/
ls
shasta
ls
Q
ls
cd ..
ls
cd assemblies/
ls
cd ..
cd data/
ls
cd ../../vangarcia/
cd assemblies/
cd flye/
cd out_Diaz26_ont/
/o
p/opt/busco/BUSCO.py
python /opt/busco/BUSCO.py
pwd
/opt/busco/BUSCO.py -i assembly.fasta -l /data/references/BUSCO/bacteria_odb9 -o ../Diaz26_flye -m genome
/opt/busco/BUSCO.py -i assembly.fasta -l /data/references/BUSCO/bacteria_odb9 -o ../Diaz11_flye -m genome
/opt/busco/BUSCO.py -i assembly.fasta -l /data/references/BUSCO/bacteria_odb9 -o Diaz11_flye -m genome
Q
shasta
/opt/busco/BUSCO.py -i assembly.fasta -l /data/references/BUSCO/bacteria_odb9 -o Diaz26_flye -m genome 
pwd
cd ../
ls
cd ..
ls
cd ..
ls
cd assemblies/
cd flye/
cd out_Diaz26_ont/
ls
quast.py
quast.py -h
/usr/bin/quast.py
/opt/quast4.6.3/quast.py -o ../../../reports/ -t 4 assembly.fasta
/opt/quast-4.6.3/quast.py -o ../../../reports/ -t 4 assembly.fasta
pwd
cd ..
cd flye/out_Diaz26_ont/
shasta --comand=assemble --memoryBacking 2M --threads 8 --input ../../data/Diaz26.ONT.fasta --assemblyDirectory Diaz26_shasta > shasta.log
ls
cd .
cd ..
cd unicycler/
export PATH=$PATH:/opt/SPAdes-3.14.0-Linux/bin
echo $PATH
unicycler -1 ../../data/Diaz-26_F_paired.fq.gz -2 ../../data/Diaz-26_R_paired.fq.gz -s ../../data/Diaz-26_unpaired.fq.gz -l ../../data/Diaz26.ONT.fastq -o unicycler_conservative --vcf --mode conservative --spades_tmp_dir /projects/temp > unicycler.log
pwd
cd ..
ls
cd ..
ls
cd data/
ls
cd ..
ls
unicycler -1 ../../data/Diaz-26_F_paired.fq.gz -2 ../../data/Diaz-26_R_paired.fq.gz -s ../../data/Diaz-26_unpaired.fq.gz -l ../../data/Diaz26.ONT.fastq -o unicycler_conservative --vcf --mode conservative --spades_tmp_dir /projects/tmp > unicycler.log
cd data/
ln -s
ln -s /projects/bacterial_genomics/assembly/Rahil/Analysis_Diaz26/Flye_Only_Diaz26/Diaz-26_F_paired.fq.gz
ln -s /projects/bacterial_genomics/assembly/Rahil/Analysis_Diaz26/Flye_Only_Diaz26/Diaz-26_R_paired.fq.gz
cd ..
cd assemblies/
cd unicycler/
unicycler -1 ../../data/Diaz-26_F_paired.fq.gz -2 ../../data/Diaz-26_R_paired.fq.gz -s ../../data/Diaz-26_unpaired.fq.gz -l ../../data/Diaz26.ONT.fastq -o unicycler_conservative --vcf --mode conservative --spades_tmp_dir /projects/tmp > unicycler.log
cd ..
cd data/
ls
pwd
cd ..
pwd
mkdir Data
ls
cd Data/
ln -s /projects/bacterial_genomics/assembly/Rahil/Analysis_Diaz26/Flye_Only_Diaz26/Diaz-26_R_paired.fq.gz
ln -s /projects/bacterial_genomics/assembly/Rahil/Analysis_Diaz26/Flye_Only_Diaz26/Diaz-26_F_paired.fq.gz
ln -s /projects/bacterial_genomics/assembly/Rahil/Analysis_Diaz26/Flye_Only_Diaz26/Diaz-26.all.fastq
cd ..
cd assemblies/
cd shasta/
ls
shasta --comand=assemble --memoryBacking 2M --threads 8 --input ../../Data/Diaz26.ONT.fasta --assemblyDirectory Diaz26_shasta > shasta.log
ls
cd ..
cd unicycler/
unicycler -1 ../../Data/Diaz-26_F_paired.fq.gz -2 ../../Data/Diaz-26_R_paired.fq.gz -s ../../Data/Diaz-26_unpaired.fq.gz -l ../../Data/Diaz26.ONT.fastq -o unicycler_conservative --vcf --mode conservative --spades_tmp_dir /projects/tmp > unicycler.log
cd ..
cd Data
ls
zcat Diaz-26_
zcat Diaz-26_F_unpaired.fq.gz Diaz-26_R_unpaired.fq.gz | gzip > Diaz-26_unpaired.fq.gz
ls
cd ..
cd data/
cd ..
cd assemblies/
cd unicycler/
export PATH=$PATH:/opt/SPAdes-3.14.0-Linux/bin
echo $PATH
unicycler -1 ../../data/Diaz-26_F_paired.fq.gz -2 ../../data/Diaz-26_R_paired.fq.gz -s ../../data/Diaz-11_unpaired.fq.gz -l ../../data/Diaz26.ONT.fastq -o unicycler_conservative --vcf --mode conservative --spades_tmp_dir /projects/temp > unicycler.log
unicycler -1 ../../data/Diaz-26_F_paired.fq.gz -2 ../../data/Diaz-26_R_paired.fq.gz -s ../../data/Diaz-26_unpaired.fq.gz -l ../../data/Diaz26.ONT.fastq -o unicycler_conservative --vcf --mode conservative --spades_tmp_dir /projects/temp > unicycler.log
ls
cd unicycler_conservative/
ls
cd ../../
ls
cd shasta/
ls
cd Diaz26_shastaV2/
ls
cat Assembly.fasta 
cd ..
rm -rf Diaz26_shastaV2/
ls
cd Diaz26_shasta/
ls
cd Ass
cd Assembly.fasta
cat Assembly.fasta 
cd ..
screen
cd migrate_directory/
ls
cd mytilus_1
ls
cd Stepping_Stone_Model/
screen -S vangarcia
screen -r vangarcia
ls
cd ..
ls
nano base_parm.txt 
ls
cp -r Stepping_Stone_Model/ Stepping_Stone_Model2/
ls
cd Stepping_Stone_Model2/
ls
nano parmfile 
screen -S vangarcia2
 cd ..
cp -r Stepping_Stone_Model2/ Stepping_Stone_Model3/
cd Stepping_Stone_Model3/
nano parmfile 
screen -r vangarcia3

cd ..
cp -r Stepping_Stone_Model2/ Stepping_Stone_Model4/
cd Stepping_Stone_Model4/
ls
nano parmfile 
screen -S vangarcia4
screen -r vangarcia
cd ..
ls
cd ..
ls
cd migrate_directory/
ls
cp -r mytilus_1/ mytilus_2/
ls
cd mytilus_2/
ls
rm -rf mytilus_1/
cd ..
cp -r mytilus_1/ mytilus_2/
cd mytilus_2/
ls
mv mytilus_1/ test_replicates_two
ls
cd test_replicates_two/
ls
cd ../..
ls
cp -r mytilus_2/ mytilus_3/
cd mytilus_3/
ls
mv mytilus_2/ test_replicates_three
ls
cd test_replicates_three/
ls
cd ..
ls
cd ..
rm -rf mytilus_3
ls
cd mytilus_2
ls
cd test_replicates_two/
ls
cd ../..
mkdir mytilus_3
cd mytilus_3
ls
cd ..
cp -rf mytilus_2 mytilus_3
ls
cd mytilus_3
ls
mv -rf mytilus_2 ./mytilus_3/
mv -r mytilus_2 ./mytilus_3/
mv mytilus_2 ./mytilus_3/
ls
cd ..
ls
rm -rf mytilus_3
mkdir mytilus_3
cd mytilus_2
ls
cd test_replicates_two/
ls
cd ..
ls
rm -rf mytilus_3
cp -r mytilus_2 mytilus_3
ls
cd mytilus_3
ls
mv test_replicates_two/ test_replicates_three
ls
cd test_replicates_three/
ls
cd ..
cd mytilus_2
ls
cd test_replicates_two/
ls
cd Stepping_Stone_Model
screen -R vangarcia5
cd ../Stepping_Stone_Model2
screen -R vangarcia6
cd ../Stepping_Stone_Model3
screen -R vangarcia7
cd ../Stepping_Stone_Model4
screen -R vangarcia8
cd ../..
cd ..
cd mytilus_3
ls
cd test_replicates_three/
ls
cd Stepping_Stone_Model
screen -R vangarcia9
cd ../Stepping_Stone_Model2
screen -R vangarcia10
cd ../Stepping_Stone_Model3
screen -R vangarcia11
cd ../Stepping_Stone_Model4
screen -R vangarcia12
ls
screen -r
screen -r 21307.vangarcia7
pwd
ls
ls
cd migrate_directory/
ls
cd ..
mkdir migrate_directory_210531
ls
pwd
cd migrate_directory_210531/
ls
cd migrate_for_evolution/
ls
screen ls
htop
./run_migrate.sh 
ls
less run_migrate.sh 
bash
./run_migrate.sh 
ls */
migrate-n
dos2unix
dos2unix -n run_migrate.sh run_migrate2.sh 
./run_migrate2.sh 
screen ls
screen -ls
screen -r 21156
ls
mv parmfile.txt parmfile
ls
migrate-n
nano parmfile 
migrate-n
nano parmfile 
migrate-n
cd
ls
cd migrate_directory/
ls
mv mytilus_1/mcalifornianus_210528.mig ./mcalifornianus_210528.mig
ls
cd mytilus_1
cd Stepping_Stone_Model/
migrate-n
nano
ls
nano parmfile 
migrate-n
nano parmfile 
migrate-n
migrate-n
migrate-n
migrate-n
migrate-n
migrate-n
migrate-n
migrate-n
migrate-n
migrate-n
migrate-n
migrate-n
ls
mv migrate_directory_210531 bad_migrate_directory_210531
ls
pwd
ls
cd migrate_for_evolution/
ls
dos2unix -h run_migrate.sh run_migrate2.sh 
screen -ls
pkill screen
screen -ls
./run_migrate2.sh 
ls
./run_migrate.sh 
chmod +x
chmod +x run_migrate.sh 
./run_migrate.sh 
screen -ls
screen -r 1922
screen -r 1900
screen -r 1844
screen -r 1791
screen -r 1812
