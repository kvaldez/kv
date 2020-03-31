import re
fileNames = '/Volumes/BioDiscovery/Lab/Tofilon/P16_OB_RH/fileNames.txt'
base_path = '/data/valdezkm/P16_OB_RH/'

with open(fileNames) as f:
    fileNames = f.read()
    fileNames = fileNames.split('\n')

    even = fileNames[0:][::2]
    odd = fileNames[1:][::2]
    star_1pass = open('swarm_files/star_1pass.swarm','w')
    star_2pass = open('swarm_files/star_2pass.swarm','w')
    mouse_star_1pass = open('swarm_files/mouse_star_1pass.swarm','w')
    mouse_star_2pass = open('swarm_files/mouse_star_2pass.swarm','w')
    # htseq = open('swarm_files/htseq.swarm','w')
    fastQC = open('swarm_files/fastQC.swarm','w')
    fastq_screen = open('swarm_files/fastq_screen.swarm', 'w')
    trim = open('swarm_files/trim.swarm','w')
    fastQC_trimmed = open('swarm_files/fastQC_trimmed.swarm','w')
    sort = open('swarm_files/sort.swarm','w')
    bamcmp = open('swarm_files/bamcmp.swarm','w')
    mergeBAMs = open('swarm_files/mergeBAMs.swarm','w')
    sortCoord = open('swarm_files/sortCoord.swarm','w')
    filterUnpaired = open('swarm_files/filterUnpaired.swarm','w')
    querySort = open('swarm_files/querySort.swarm','w')
    bam2fastq = open('swarm_files/bam2fastq.swarm','w')
    screen_again = open('swarm_files/screenAgain.swarm','w')
    # DEXseq = open('swarm_files/DEXseq.swarm','w')

    for even,odd in zip(even,odd):
        short = re.split('_|\.',even)
        short = short[0] + '_' + short[1]
        fastQC.write('fastqc -f fastq ' + base_path + 'fastq_files/' + even + ' \\\n'
                     + base_path + 'fastq_files/' + odd + ' \\\n'
                     + '-o ' + base_path + 'QC/' '\n')
        trim.write(
            'java -Xmx2g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar /usr/local/apps/trimmomatic/Trimmomatic-0.33/trimmomatic-0.33.jar PE -threads 32 -phred33 \\\n'
            + base_path + 'fastq_files/' + even + ' ' + base_path + 'fastq_files/' + odd + ' \\\n'
            + base_path + 'fastq_files/' + short + '_R1_trimmed.fastq.gz ' + base_path + 'fastq_files/' + short + '_R1_unpair.fastq.gz \\\n'
            + base_path + 'fastq_files/' + short + '_R2_trimmed.fastq.gz ' + base_path + 'fastq_files/' + short + '_R2_unpair.fastq.gz \\\n'
            + 'ILLUMINACLIP:/data/CCBR_Pipeliner/db/PipeDB/dev/TruSeq_and_nextera_adapters_new.fa:2:36:10 LEADING:10 TRAILING:10 MINLEN:35 \n')
        fastq_screen.write('cd ' + base_path + 'fastq_files; \\\n'
                           + '/data/CCBR_Pipeliner/db/PipeDB/bin/fastq_screen_v0.9.3/fastq_screen \\\n'
                           + '--conf /data/CCBR_Pipeliner/db/PipeDB/lib/fastq_screen.conf --subset 0 --aligner bowtie2 \\\n'
                           + base_path + 'fastq_files/' + short + '_R1_trimmed.fastq.gz ' + base_path + 'fastq_files/' + short + '_R2_trimmed.fastq.gz \\\n'
                           + '--outdir ' + base_path + 'QC/ \n')
        fastQC_trimmed.write('fastqc -f fastq ' + base_path + 'fastq_files/' + short + '_R1_trimmed.fastq.gz ' + ' \\\n'
                     + base_path + 'fastq_files/' + short + '_R2_trimmed.fastq.gz ' + ' \\\n'
                     + '-o ' + base_path + 'QC/' '\n')
        mouse_star_1pass.write('cd ' + base_path + 'fastq_files \\\n&& mkdir ../Mouse_aligned_1pass/' + short
                         + ' \\\n&& STAR \\\n--runThreadN $SLURM_CPUS_PER_TASK \\\n'
                         + '--genomeDir /fdb/STAR_indices/2.7.0f/GENCODE/Gencode_mouse/release_M16/genes-150 \\\n'
                         + '--outSAMtype BAM Unsorted \\\n'
                         + '--outSAMattributes NM AS \\\n'
                         + '--outFileNamePrefix ../Mouse_aligned_1pass/' + short + '/ \\\n'
                         + '--readFilesIn ' +  short + '_R1_trimmed.fastq.gz ' + short + '_R2_trimmed.fastq.gz ' + ' ' + '\\\n'
                         + '--readFilesCommand zcat\n')
        mouse_star_2pass.write('cd ' + base_path + 'fastq_files \\\n&& mkdir ../Mouse_aligned_2pass/' + short
                         + ' \\\n&& STAR \\\n--runThreadN $SLURM_CPUS_PER_TASK \\\n'
                         + '--genomeDir /fdb/STAR_indices/2.7.0f/GENCODE/Gencode_mouse/release_M16/genes-150 \\\n'
                         + '--outSAMtype BAM SortedByCoordinate \\\n'
                         + '--outFileNamePrefix ../Mouse_aligned_2pass/' + short + '/ \\\n'
                         + '--outSAMattributes NM AS \\\n'
                         + '--sjdbGTFfile /fdb/STAR_indices/2.7.0f/GENCODE/Gencode_mouse/release_M16/genes.gtf --sjdbOverhang 150 \\\n'
                         + '--readFilesIn ' + short + '_R1_trimmed.fastq.gz ' + short + '_R2_trimmed.fastq.gz ' + ' ' + '\\\n'
                         + '--sjdbFileChrStartEnd ' + base_path + 'Mouse_aligned_1pass/*/SJ.out.tab \\\n'
                         + '--quantMode TranscriptomeSAM GeneCounts \\\n'
                         + '--readFilesCommand zcat\n')
        star_1pass.write('cd ' + base_path + 'fastq_files \\\n&& mkdir ../Aligned_1pass/' + short
                  + ' \\\n&& STAR \\\n--runThreadN $SLURM_CPUS_PER_TASK \\\n'
                  + '--genomeDir /fdb/STAR_indices/2.7.0f/GENCODE/Gencode_human/release_27/genes-150 \\\n'
                  + '--outSAMtype BAM Unsorted \\\n'
                  + '--outSAMattributes NM AS \\\n'
                  + '--outFileNamePrefix ../Aligned_1pass/' + short + '/ \\\n'
                  + '--readFilesIn ' + short + '_R1_trimmed.fastq.gz ' + short + '_R2_trimmed.fastq.gz ' + ' ' + '\\\n'
                  + '--readFilesCommand zcat\n')
        star_2pass.write('cd ' + base_path + 'fastq_files \\\n&& mkdir ../Aligned_2pass/' + short
                        + ' \\\n&& STAR \\\n--runThreadN $SLURM_CPUS_PER_TASK \\\n'
                        + '--genomeDir /fdb/STAR_indices/2.7.0f/GENCODE/Gencode_human/release_27/genes-150 \\\n'
                        + '--outSAMtype BAM SortedByCoordinate \\\n'
                        + '--outFileNamePrefix ../Aligned_2pass/' + short + '/ \\\n'
                        + '--outSAMattributes NM AS \\\n'
                        + '--sjdbGTFfile /fdb/STAR_indices/2.7.0f/GENCODE/Gencode_human/release_27/genes.gtf --sjdbOverhang 150 \\\n'
                        + '--readFilesIn ' + short + '_R1_trimmed.fastq.gz ' + short + '_R2_trimmed.fastq.gz ' + ' ' + '\\\n'
                        + '--sjdbFileChrStartEnd ' + base_path + 'Aligned_1pass/*/SJ.out.tab \\\n'
                        + '--quantMode TranscriptomeSAM GeneCounts \\\n'
                        + '--readFilesCommand zcat\n')
        # human
        sort.write('java -Xmx4g -XX:ParallelGCThreads=5 -jar $PICARDJARPATH/picard.jar \\\n'
                          + 'SortSam \\\n'
                          + 'I=' + base_path + 'BAMs/' + short + '.bam \\\n'
                          + 'O=' + base_path + 'BAMs/' + short + '_queryNameSorted.bam \\\n'
                          + 'SORT_ORDER=queryname \n')
        # mouse
        sort.write('java -Xmx4g -XX:ParallelGCThreads=5 -jar $PICARDJARPATH/picard.jar \\\n'
                   + 'SortSam \\\n'
                   + 'I=' + base_path + 'BAMs/mouse_' + short + '.bam \\\n'
                   + 'O=' + base_path + 'BAMs/mouse_' + short + '_queryNameSorted.bam \\\n'
                   + 'SORT_ORDER=queryname \n')
        bamcmp.write('cd ' + base_path + 'BAMs; \\\n'
                     + 'bamcmp -1 ' + short + '_queryNameSorted.bam \\\n'
                     + '-2 ' + 'mouse_' + short + '_queryNameSorted.bam \\\n'
                     + '-a ' + short + '_human_only.bam \\\n'
                     + '-b ' + short + '_mouse_only.bam \\\n'
                     + '-A ' + short + '_human_better.bam \\\n'
                     + '-B ' + short + '_mouse_better.bam \\\n'
                     + '-C ' + short + '_human_worse.bam \\\n'
                     + '-D ' + short + '_mouse_worse.bam \\\n'
                     + '-t 32 -N -s as \n')
        mergeBAMs.write('cd ' + base_path + 'BAMs/human_BAMs; \\\n'
                        + 'samtools cat -o ' + short + '_merged_HUMAN.bam ' + short + '* \n')
        sortCoord.write('cd ' + base_path + 'BAMs/human_BAMs; \\\n'
                        + 'samtools sort -o ' + short + '_merged_HUMAN_sorted.bam ' + short + '_merged_HUMAN.bam \n')
        filterUnpaired.write('cd ' + base_path + 'BAMs/human_BAMs; \\\n'
                             + 'samtools view -b -f 1 \\\n'
                             + short + '_merged_HUMAN_sorted.bam \\\n'
                             + '> ' + short + '_human_filterUnpaired.bam \n')
        querySort.write('java -Xmx4g -XX:ParallelGCThreads=5 -jar $PICARDJARPATH/picard.jar \\\n'
                   + 'SortSam \\\n'
                   + 'I=' + base_path + 'BAMs/human_BAMs/' + short + '_human_filterUnpaired.bam \\\n'
                   + 'O=' + base_path + 'BAMs/human_BAMs/' + short + '_sortedAgain.bam \\\n'
                   + 'SORT_ORDER=queryname \n')
        bam2fastq.write('java -Xmx4g -XX:ParallelGCThreads=5 -jar $PICARDJARPATH/picard.jar \\\n'
                        + 'SamToFastq \\\n'
                        + 'I=' + base_path + 'BAMs/human_BAMs/' + short + '_sortedAgain.bam \\\n'
                        + 'FASTQ=' + base_path + 'fastq_files/bam2fastq/' + short + '_human.R1.fastq.gz \\\n'
                        + 'SECOND_END_FASTQ=' + base_path + 'fastq_files/bam2fastq/' + short + '_human.R2.fastq.gz \n')
        screen_again.write('cd ' + base_path + 'fastq_files/bam2fastq; \\\n'
                           + '/data/CCBR_Pipeliner/db/PipeDB/bin/fastq_screen_v0.9.3/fastq_screen \\\n'
                           + '--conf /data/CCBR_Pipeliner/db/PipeDB/lib/fastq_screen.conf --subset 0 --aligner bowtie2 \\\n'
                           + base_path + 'fastq_files/bam2fastq/' + short + '_human.R1.fastq.gz ' + base_path + 'fastq_files/bam2fastq/' + short + '_human.R2.fastq.gz \\\n'
                           + '--outdir ' + base_path + 'QC/filtered_QC/ \n')



        aligned = short + 'Aligned.sortedByCoord.out.bam'
        # htseq.write('cd ' + base_path + ' Aligned; htseq-count -m intersection-nonempty -i gene_name -s no -f bam \\\n'
        #             +  aligned + ' /fdb/STAR_indices/2.5.4a/GENCODE/Gencode_human/release_26/genes.gtf \\\n > ' + base_path + 'Counted/'
        #             + short + '.txt\n')
        # DEXseq.write('python /usr/local/apps/R/3.6/site-library_3.6.0/DEXSeq/python_scripts/dexseq_count.py -p yes -r pos -f bam ../exon_DEXseq.gff ' + '../Aligned_2pass_BAMSonly/' + aligned + ' ../Counted/' + short + '_DEXseq.count\n')


#****** Manually take out extra new line character at bottom of files

# mkdir fastq_files Aligned_1pass Aligned_2pass Counted QC
