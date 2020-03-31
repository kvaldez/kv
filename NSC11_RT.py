
def create_pipeline(base_path,names,local):

    with open(names) as f:
        fileNames = f.read()
        fileNames = fileNames.split('\n')

    even = fileNames[0:][::2]
    odd = fileNames[1:][::2]
    fastQC = open(local + 'fastQC.swarm','w')
    bwa = open(local + 'bwa.swarm','w')
    bwa_mouse = open(local + 'bwa_mouse.swarm','w')
    # trim = open(local + 'trim.swarm','w')
    fastq_screen = open(local + 'fastq_screen.swarm','w')
    query_name_sort = open(local + 'query_name_sort.swarm','w')
    query_mouse = open(local + 'query_name_sort_MOUSE.swarm','w')         #still needs to be created
    unpaired = open(local + 'unpaired.swarm','w')
    humanOnly_nameSort = open(local + 'humanOnly_nameSort.swarm','w')
    bam2fastq = open(local + 'bam2fastq.swarm','w')
    # fastq_screen_humanOnly = open(local + 'fastq_screen_humanOnly.swarm','w')
    vcf2maf = open(local + 'vcf2maf.swarm','w')
    ascat = open(local + 'ASCAT.swarm','w')
    alleleCounter = open(local + 'alleleCount_REMAINING.swarm','w')    # already did control 1
    expands = open(local + 'expands.swarm','w')

    # mark_dups = open(local + 'mark_dups.swarm','w')
    # coord_sort = open(local + 'coord_sort.swarm','w')
    # addRG = open(local + 'addRG.swarm','w')
    bamcmp = open(local + 'bamcmp.swarm','w')



    for even,odd in zip(even,odd):
        even = even.strip()
        odd = odd.strip()
        short = even.split('.')
        short = short[0]
        fastQC.write('fastqc -f fastq ' + base_path + 'fastq_files/' + even + ' \\\n'
                     + base_path + 'fastq_files/' + odd + ' \\\n'
                     + '-o ' + base_path + 'QC/' '\n')
        bwa.write('cd ' + base_path + 'fastq_files; \\\n'
                  + 'bwa mem -M -t 12 /data/valdezkm/bwa_human_index_hg19_ucsc/ucsc.hg19.fasta \\\n'
                  + even + ' ' + odd +  ' | samtools view -bS - > ' + base_path + 'Aligned/' + short + '.bam \n')
        bwa_mouse.write('cd ' + base_path + 'fastq_files; \\\n'
                        + 'bwa mem -M -t 12 /fdb/bwa/indexes/mm9.fa \\\n'
                        + even + ' ' + odd + ' | samtools view -bS - > ' + base_path + 'Aligned/' + short + '_MOUSE.bam \n')
        fastq_screen.write('cd ' + base_path + 'fastq_files; \\\n'
                           + '/data/CCBR_Pipeliner/db/PipeDB/bin/fastq_screen_v0.9.3/fastq_screen \\\n'
                           + '--conf /data/CCBR_Pipeliner/db/PipeDB/lib/fastq_screen.conf --subset 0 --aligner bowtie2 \\\n'
                           + base_path + 'fastq_files/' + even + ' ' + base_path + 'fastq_files/' + odd + ' \\\n'
                           + '--outdir ' + base_path + 'QC/ \n')
        query_name_sort.write('java -Xmx4g -XX:ParallelGCThreads=5 -jar $PICARDJARPATH/picard.jar \\\n'
                              + 'SortSam \\\n'
                              + 'I=' + base_path + 'Aligned/' + short + '.bam \\\n'
                              + 'O=' + base_path + 'Aligned/' + short + '_queryNameSorted.bam \\\n'
                              + 'SORT_ORDER=queryname \n')
        query_mouse.write('java -Xmx4g -XX:ParallelGCThreads=5 -jar $PICARDJARPATH/picard.jar \\\n'
                              + 'SortSam \\\n'
                              + 'I=' + base_path + 'Aligned/' + short + '_MOUSE.bam \\\n'
                              + 'O=' + base_path + 'Aligned/' + short + '_queryNameSorted_MOUSE.bam \\\n'
                              + 'SORT_ORDER=queryname \n')
        bamcmp.write('cd ' + base_path + 'Aligned; \\\n'
                     + 'bamcmp -1 ' + short + '_queryNameSorted.bam \\\n'
                     + '-2 ' + short + '_queryNameSorted_MOUSE.bam \\\n'
                     + '-a ' + short + '_human_only.bam \\\n'
                     + '-b ' + short + '_mouse_only.bam \\\n'
                     + '-A ' + short + '_human_better.bam \\\n'
                     + '-B ' + short + '_mouse_better.bam \\\n'
                     + '-C ' + short + '_human_worse.bam \\\n'
                     + '-D ' + short + '_mouse_worse.bam \\\n'
                     + '-t 32 -N -s as \n')
        unpaired.write('cd ' + base_path + 'Aligned; samtools view -b -f 1 \\\n'
                       + short + '_human_better.bam \\\n'
                       + '> ' + short + '_human_better_filterPaired.bam \n')
        humanOnly_nameSort.write('java -Xmx4g -XX:ParallelGCThreads=5 -jar $PICARDJARPATH/picard.jar \\\n'
                                 + 'SortSam \\\n'
                                 + 'I=' + base_path + 'Aligned/' + short + '_human_better_filterPaired.bam \\\n'
                                 + 'O=' + base_path + 'Aligned/' + short + '_human_better_filterPaired_nameSorted.bam \\\n'
                                 + 'SORT_ORDER=queryname \n')
        bam2fastq.write('java -Xmx4g -XX:ParallelGCThreads=5 -jar $PICARDJARPATH/picard.jar \\\n'
                        + 'SamToFastq \\\n'
                        + 'I=' + base_path + 'Aligned/' + short + '_human_better_filterPaired_nameSorted.bam \\\n'
                        + 'FASTQ=' + base_path + 'fastq_files/Picard_sam2fastq/' + short + '_human_better.R1.fastq.gz \\\n'
                        + 'SECOND_END_FASTQ=' + base_path + 'fastq_files/Picard_sam2fastq/' + short + '_human_better.R2.fastq.gz \n')
        vcf2maf.write('vcf2maf.pl --input-vcf /data/valdezkm/Real_Tofilon/output/mutect2_out/' + short + '.FINALmutect2.vcf \\\n'
                      + '--output-maf /data/valdezkm/Real_Tofilon/mafs_handmade/' + short + '.maf \\\n'
                      + '--vep-path $VEP_HOME --vep-data $VEPCACHEDIR --ref-fasta /data/CCBR_Pipeliner/db/PipeDB/lib/hs37d5.fa \\\n'
                      + '--filter-vcf /fdb/VEP/88/cache/ExAC.r0.3.sites.vep.vcf.gz --vep-forks 2 --vcf-tumor-id ' + short + ' \\\n'
                      + '--tumor-id ' + short + ' --ncbi-build GRCh37 --species homo_sapiens \n')
        if ('control' in str(short) or 'RT' in str(short)):
            ascat.write('mkdir /data/valdezkm/Real_Tofilon/ASCAT/' + short + '; \\\n'
                        'ascat.pl \\\n'
                        + '-outdir /data/valdezkm/Real_Tofilon/ASCAT/' + short +  ' \\\n'
                        + '-tumour /data/valdezkm/Real_Tofilon/output_IVmatchedNormal/' + short + '_human_better.recal.bam \\\n'
                        + '-normal /data/valdezkm/Real_Tofilon/output_IVmatchedNormal/in_vitro_1_S21_human_better.recal.bam \\\n'
                        + '-reference /data/CCBR_Pipeliner/db/PipeDB/lib/hs37d5.fa \\\n'
                        + '-snp_gc /data/valdezkm/Real_Tofilon/ASCAT/SnpGcCorrections.tsv \\\n'
                        + '-protocol WXS \\\n'
                        + '-platform ILLUMINA \\\n'
                        + '-gender L \\\n'
                        + '-locus /data/valdezkm/Real_Tofilon/ASCAT/GRCh37d5_Y.loci \\\n'
                        + '-species Human \\\n'
                        + '-assembly GRCh37d5 \\\n'
                        + '-cpus 8\n')
        # already did alleleCounter with control 1, remove from swarm
        alleleCounter.write('alleleCounter.pl -b /data/valdezkm/Real_Tofilon/output_IVmatchedNormal/' + short + '_human_better.recal.bam \\\n'
                            + '-o /data/valdezkm/Real_Tofilon/alleleCounter/' + short + '.counts \\\n'
                            + '-l /data/valdezkm/Real_Tofilon/alleleCounter/loci.txt \n')
        expands.write('Rscript /data/valdezkm/Real_Tofilon/expands_biowulf.R ' + '/data/valdezkm/Real_Tofilon/output_noM/mutect2_out/oncotator_out/' + short + '_human_better.maf \\\n'
                      + '/data/valdezkm/Real_Tofilon/output_noM/cnvkit_out/' + short + '_human_better_calls.cns ' + '/data/valdezkm/Real_Tofilon/noMouse_expands/' + short + '\n')



#First argument, path is the base path where the files will be located in biowulf, add a trailing /
#Within parent directory:
#  $ mkdir fastq_files; mkdir Aligned; mkdir QC; mkdir SNP_calls

#Second argument is the file listing names of fastq files, in order S1.R1 S1.R2 S2.R1 S2.R2 etc

#Third argument is local base path for swarm output files, add a trailing /

#ORIGINAL DATA - change short name assignments above (base_path,names,local)
create_pipeline('/data/valdezkm/Real_Tofilon_RT/','/Users/valdezkm/Documents/Real_Tofilon_RT/fileNames.txt','/Users/valdezkm/Documents/Real_Tofilon_RT/swarm_files/')


#HUMAN ONLY DATA - change short name assignments above
# create_pipeline('/data/valdezkm/PDX_HumanOnly/','/Users/valdezkm/Documents/PDX_HumanOnly/filenames_humanOnly.txt','/data/valdezkm/PDX_Tofilon/S04380219_V5-UTR_covered.bed','/data/valdezkm/PDX_Tofilon/PoN/siteonly.vcf','/Users/valdezkm/Documents/PDX_HumanOnly/')


# <------------> Run expands on filtered mafs <------------>
expands_filtered = open('/Users/valdezkm/Documents/Real_Tofilon_RT/swarm_files/expands_filtered.swarm', 'w')
with open('/Users/valdezkm/Documents/Real_Tofilon_RT/noMouse_expands_filteredMAFs/fileNames_filtered.txt') as f:
    snvs = f.read().strip()
    snvs = snvs.split('\n')

with open('/Users/valdezkm/Documents/Real_Tofilon_RT/noMouse_expands_filteredMAFs/fileNames_cns.txt') as g:
    cnvs = g.read().strip()
    cnvs = cnvs.split('\n')

    for s,c in zip(snvs,cnvs):
        expands_filtered.write('Rscript /data/valdezkm/Real_Tofilon/expands_biowulf.R ' + '/data/valdezkm/Real_Tofilon/noMouse_expands_filtered/maf_cnv/' + s + ' \\\n'
                               + '/data/valdezkm/Real_Tofilon/noMouse_expands_filtered/maf_cnv/' + c + ' /data/valdezkm/Real_Tofilon/noMouse_expands_filtered/' + s.split('.')[0] + '\n')









