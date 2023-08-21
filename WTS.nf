

// TODO:
// - ADD MULTIQC TO STAR ALIGNMENT 
// - MAKE THE CODE EASY TO ACCESS (universal path variable, length in fastp, etc.)
// - COMMENT EVERYTHING

/*
 * Define output folders
 */ 

fastq_in = '/archive/projects/MOBIT/WTS/fasta'
//fastq_in = '/home/ajan/mobit_wts/subset_two_large'
//fastq_in = '/home/ajan/mobit_wts/subset_rnaseq'
fastqc_raw_out = '/home/ajan/mobit_wts/fastqc_raw'
fastp_out = '/home/ajan/mobit_wts/fastp'
fastqc_fastp_out = '/home/ajan/mobit_wts/fastp_fastqc'
sortmerna_out = '/home/ajan/mobit_wts/sortmerna'

rRNA_1 = '/home/ajan/mobit_wts/rRNA_databases/rfam-5.8s-database-id98.fasta'
rRNA_2 = '/home/ajan/mobit_wts/rRNA_databases/rfam-5s-database-id98.fasta'
rRNA_3 = '/home/ajan/mobit_wts/rRNA_databases/silva-arc-16s-id95.fasta'
rRNA_4 = '/home/ajan/mobit_wts/rRNA_databases/silva-arc-23s-id98.fasta'
rRNA_5 = '/home/ajan/mobit_wts/rRNA_databases/silva-bac-16s-id90.fasta'
rRNA_6 = '/home/ajan/mobit_wts/rRNA_databases/silva-bac-23s-id98.fasta'
rRNA_7 = '/home/ajan/mobit_wts/rRNA_databases/silva-euk-18s-id95.fasta'
rRNA_8 = '/home/ajan/mobit_wts/rRNA_databases/silva-euk-28s-id98.fasta'

sortmerna_fastqc = '/home/ajan/mobit_wts/sortmerna_fastqc'

repair_bbmap_out = '/home/ajan/mobit_wts/repair_bbmap'

star_index_in = '/home/ajan/mobit_wts/reference_files/genomeDir'
star_alignments_out = '/home/ajan/mobit_wts/star_alignments'
gtf_in = '/home/ajan/mobit_wts/reference_files/Homo_sapiens.GRCh38.104.gtf'
qualimap_out = '/home/ajan/mobit_wts/qualimap'
picard_out = '/home/ajan/mobit_wts/picard'
multiqc_fastqc_raw_out = '/home/ajan/mobit_wts/multiqc/fastqc_raw'
multiqc_fastqc_fastp_out = '/home/ajan/mobit_wts/multiqc/fastqc_fastp'
multiqc_fastp_out = '/home/ajan/mobit_wts/multiqc/fastp'
repair_fastqc_out = '/home/ajan/mobit_wts/fastqc_repaired_sortmerna'
multiqc_star_alignments_out = '/home/ajan/mobit_wts/multiqc/star_alignments'
multiqc_qualimap_out = '/home/ajan/mobit_wts/multiqc/qualimap'
multiqc_picard_out = '/home/ajan/mobit_wts/multiqc/picard'
ref_flat_in = '/home/ajan/mobit_wts/reference_files/noChr_refFlat.txt'
samtools_index_out = star_alignments_out
samtools_flagstat_out_path = '/home/ajan/mobit_wts/samtools_flagstat'
sortmerna_multiqc_out = '/home/ajan/mobit_wts/multiqc/sortmerna'
multiqc_samtools_flagstat_out = '/home/ajan/mobit_wts/multiqc/samtools_flagstat'
multiqc_repair_fastqc_out = '/home/ajan/mobit_wts/multiqc/repaired_sortmerna'
multiqc_featureCounts_out = '/home/ajan/mobit_wts/multiqc/featureCounts'
featureCounts_out = '/home/ajan/mobit_wts/featureCounts'

log.info """\


====================================================================================================================

 Transcriptomics 0.2a    

=====================================================================================================================
|                                                                                                                   
| .fastq files                           : $fastq_in                                                                
| STAR index                             : $star_index_in                                                           
| GTF                                    : $gtf_in
| refFlat                                : $ref_flat_in 
| Head output dir                        : /home/ajan/mobit_wts                              
|                                                                                                                     
=====================================================================================================================


"""

/*
 *  Parse the input parameters
 */ 
Channel
.fromFilePairs("$fastq_in/*_{R1,R2}_001.fastq.gz", flat: true)
.into{ raw_fastq_fastqc; raw_fastq_fastp}

/*
 * Process 1: Fastqc on raw fastq
 */

process fastqc_raw {

    label 'fastqc_raw'

    module = 'MOD/fastqc/0.11.9'

    publishDir "$fastqc_raw_out", mode:'copy'
   
    input:
        set val(id), file(fastq_1), file(fastq_2) from raw_fastq_fastqc
    
    output:
        file ("*.zip") into fastqc_multiqc
        file "*.html"

    script:
        """
        fastqc -t ${task.cpus} ${fastq_1} ${fastq_2}
	"""
}

/*
 * Process 2: Filtering
 */
 
process fastp {

    label 'fastp'

    module 'MOD/fastp/0.20.1'
    
    publishDir "$fastp_out", mode:'copy', pattern: '*.fastq.gz'
    publishDir "$fastp_out", mode:'copy', pattern: '*fastp.json'
    publishDir "$fastp_out", mode:'copy', pattern: '*fastp.html'

    input:
        set val(id), file(fastq_1), file(fastq_2) from raw_fastq_fastp
    
    output:
        set val(id), file("${id}_1_fastp.fastq.gz"), file("${id}_2_fastp.fastq.gz") into (fastp_fastqc, fastp_sortmerna)
        file("${id}_fastp.json") into fastp_multiqc
        file "${id}_fastp.html"

    script:
        """
	fastp -i ${fastq_1} -I ${fastq_2} -o ${id}_1_fastp.fastq.gz -O ${id}_2_fastp.fastq.gz -j ${id}_fastp.json -h ${id}_fastp.html -q 30 -e 25 -n 3 -l 65 -c -x -p -w ${task.cpus}

	"""
}  

 /*
  * Process 3: Fastqc on filtered
  */

 process fastqc_fastp {

     label 'fastqc_fastp'

     module = 'MOD/fastqc/0.11.9'

     publishDir "$fastqc_fastp_out", mode:'copy'
 
     input:
         set val(id), file(fastq_1), file(fastq_2) from fastp_fastqc
  
     output:
         file("*.zip") into fastqc_fastp_multiqc
         file "*.html"

     script:
         """
         fastqc -t ${task.cpus} ${fastq_1} ${fastq_2}
 	"""
 }

 /*
  * Process 3: Remove rRNA
  */

 process sortmerna {

     label 'sortmerna'

     conda '/home/ajan/.conda/envs/sortmerna'

     publishDir "$sortmerna_out", mode:'copy'
 
     input:
         set val(id), file(fastq_1), file(fastq_2) from fastp_sortmerna
  
     output:
         set val(id), file("clean_${id}_fwd.fq.gz"), file("clean_${id}_rev.fq.gz") into (sortmerna_repair)
         file("*.log") into sortmerna_out_path

     script:
         """

         mkdir /home/ajan/mobit_wts/sortmerna/${id}

         sortmerna \
         --ref ${rRNA_1} \
         --ref ${rRNA_2} \
         --ref ${rRNA_3} \
         --ref ${rRNA_4} \
         --ref ${rRNA_5} \
         --ref ${rRNA_6} \
         --ref ${rRNA_7} \
         --ref ${rRNA_8} \
         --reads ${fastq_1} \
         --reads ${fastq_2} \
         --workdir /home/ajan/mobit_wts/sortmerna/${id} \
         -a ${task.cpus} \
         --fastx \
         -v \
         -m 31744 \
         --out2 \
         --other clean_${id} \
         --aligned aligned_rRNA_${id}
 	"""
 }

 /*
  * Process 4: Run STAR 
  */

 process repair {

     label 'repair'

     conda '/home/ajan/.conda/envs/bbmap'

     publishDir "$repair_bbmap_out", mode:'copy'
   
     input:
         set val(id), file(fastq_1), file(fastq_2) from sortmerna_repair
  
     output:
         set val(id), file("repaired_clean_${id}_fwd.fq.gz"), file("repaired_clean_${id}_rev.fq.gz") into (repair_star, repair_fastqc)
         file("*")

     script:
         """
         repair.sh in=${fastq_1} in2=${fastq_2} out=repaired_clean_${id}_fwd.fq.gz out2=repaired_clean_${id}_rev.fq.gz 
 	"""
 }

  /*
  * Process 3: Fastqc on filtered
  */

 process repair_fastqc {

     label 'repair_fastqc'

     module = 'MOD/fastqc/0.11.9'

     publishDir "$repair_fastqc_out", mode:'copy'
 
     input:
         set val(id), file(fastq_1), file(fastq_2) from repair_fastqc
  
     output:
         file("*.zip") into repair_fastqc_multiqc
         file "*.html"

     script:
         """
         fastqc -t ${task.cpus} ${fastq_1} ${fastq_2}
 	"""
 }

 /*
  * Process 4: Run STAR 
  */

 process star_align {

     label 'star_align'

     module 'MOD/star/2.7.9a'

     publishDir "$star_alignments_out", mode:'copy'
   
     input:
         set val(id), file(fastq_1), file(fastq_2) from repair_star
    
     output:
         set val(id), file ("${id}*.bam") into (star_qualimap, star_samtools_index, star_samtools_flagstat, star_picard, star_subread, star_featureCounts)
         path ("${id}*") into star_alignments_path

     script:
         """
         STAR --runThreadN ${task.cpus} --runMode alignReads --genomeDir $star_index_in --readFilesIn ${fastq_1} ${fastq_2} --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --sjdbGTFfile $gtf_in --sjdbOverhang 75 --outFileNamePrefix ${id}

 	"""
 }

 /*
  * Process 5: Run Qualimap 
  *
  * Process *5*: Run samtools sort by name for qualimap to omit // this is possible but to tidy the code up option:
  *                                                           // -Djava.io.tmpdir=<my/temp/path>
  *                                                           // https://www.biostars.org/p/42613/
  *                                                           // will be used to set temp directory for bam files for the time of sorting
  */

 process qualimap {

     label 'qualimap'

     conda '/home/ajan/.conda/envs/qualimap'

     publishDir "$qualimap_out", mode:'copy'
   
     input:
         set val(id), file(bam) from star_qualimap
    
     output:
         path id into qualimap_out_path

     script:
         """
         export JAVA_OPTS="-Djava.io.tmpdir=/home/ajan/mobit_wts/trash"

         qualimap rnaseq -bam ${bam} -gtf ${gtf_in} -outdir ${id}/ -pe -p strand-specific-reverse -outformat PDF:HTML --java-mem-size=16G
	"""
 }

 process picard_matrix {

     label 'picard_matrix'

     conda '/home/ajan/.conda/envs/picard'

     publishDir "$picard_out", mode:'copy'
   
     input:
         set val(id), file(bam) from star_picard
    
     output:
         path id into picard_out_path

     script:
         """
         picard CollectRnaSeqMetrics -I ${bam} -O ${id} --REF_FLAT ${ref_flat_in} --STRAND SECOND_READ_TRANSCRIPTION_STRAND
	"""
 }

 process samtools_index {

     label 'samtools_index'

     conda '/home/ajan/.conda/envs/samtools'

     publishDir "$samtools_index_out", mode:'copy'
   
     input:
         set val(id), file(bam) from star_samtools_index
    
     output:
         file("*")

     script:
         """
         samtools index -@ 4 ${bam}
	"""
 }

 process samtools_flagstat {

     label 'samtools_flagstat'

     conda '/home/ajan/.conda/envs/samtools'

     publishDir "$samtools_flagstat_out_path", mode:'copy'
   
     input:
         set val(id), file(bam) from star_samtools_flagstat
    
     output:
         path ("${id}.txt") into samtools_flagstat_out

     script:
         """
         samtools flagstat ${bam} -@ 4 > ${id}.txt
	"""
 }

 process featureCounts {

     label 'featureCounts'

     conda '/home/ajan/.conda/envs/subread'

     publishDir "$featureCounts_out", mode:'copy'
   
     input:
         set val(id), file(bam) from star_featureCounts
    
     output:
         path ("*.summary") into featureCounts_multiqc
         file ("*")

     script:
         """
         featureCounts -t gene -s 2 -T ${task.cpus} --verbose -p -a ${gtf_in} -o ${id}_featureCounts_matrix.txt ${bam}
	"""
 }

/*
 * Process 5: Run multiqc on qualimap
 */

process multiqc_featureCounts {

    label 'multiqc_featureCounts'

    conda '/home/ajan/.conda/envs/multiqc'

    publishDir "$multiqc_featureCounts_out", mode:'copy'
   
    input:
        file ("*") from featureCounts_multiqc.collect()
    
    output:
        file("multiqc_featureCounts.html")

    script:
        """
        multiqc . -n multiqc_featureCounts
	"""
}

process multiqc_samtools_flagstat {

    label 'multiqc_samtools_flagstat'

    conda '/home/ajan/.conda/envs/multiqc'

    publishDir "$multiqc_samtools_flagstat_out", mode:'copy'
   
    input:
        file("*") from samtools_flagstat_out.collect()
    
    output:
        file("multiqc_samtools_flagstat.html")

    script:
        """
        multiqc . -n multiqc_samtools_flagstat
	"""
}

/*
 * Process 5: Run multiqc on qualimap
 */

process multiqc_qualimap {

    label 'multiqc_qualimap'

    conda '/home/ajan/.conda/envs/multiqc'

    publishDir "$multiqc_qualimap_out", mode:'copy'
   
    input:
        file("*") from qualimap_out_path.collect()
    
    output:
        file("multiqc_qualimap.html")

    script:
        """
        multiqc . -n multiqc_qualimap
	"""
}

/*
 * Process 5: Run multiqc on picard
 */

process multiqc_picard {

    label 'multiqc_picard'

    conda '/home/ajan/.conda/envs/multiqc'

    publishDir "$multiqc_picard_out", mode:'copy'
   
    input:
        file("*") from picard_out_path.collect()
    
    output:
        file("multiqc_picard.html")

    script:
        """
        multiqc . -n multiqc_picard
	"""
}

/*
 * Process 5: Run multiqc on picard
 */

process sortmerna_multiqc {

    label 'sortmerna_multiqc'

    conda '/home/ajan/.conda/envs/multiqc'

    publishDir "$sortmerna_multiqc_out", mode:'copy'
   
    input:
        file("*") from sortmerna_out_path.collect()
    
    output:
        file("multiqc_sortmerna.html")

    script:
        """
        multiqc . -n multiqc_sortmerna
	"""
}

/* 
 * Process 5: Run multiqc on fastp filtered fastq
 */

process multiqc_fastp {

    label 'multiqc_fastp'

    conda '/home/ajan/.conda/envs/multiqc'

    publishDir "$multiqc_fastp_out", mode:'copy'
   
    input:
        file("*") from fastp_multiqc.collect()
    
    output:
        file("multiqc_fastp.html")

    script:
        """
        multiqc *.json -n multiqc_fastp
	"""
}

/*
 * Process 6: Run multiqc on raw fastq
 */

process multiqc_fastqc_raw {

    label 'multiqc_fastqc_raw'

    conda '/home/ajan/.conda/envs/multiqc'

    publishDir "$multiqc_fastqc_raw_out", mode:'copy'
   
    input:
        file("*") from fastqc_multiqc.collect()
    
    output:
        file("multiqc_fastqc_raw.html")

    script:
        """
        multiqc *.zip

        mv multiqc_report.html multiqc_fastqc_raw.html
	"""
}

/*
 * Process 7: Run multiqc on fastp filtered FastQC
 */

process multiqc_fastqc_fastp {

    label 'multiqc_fastqc_fastp'

    conda '/home/ajan/.conda/envs/multiqc'

    publishDir "$multiqc_fastqc_fastp_out", mode:'copy'
   
    input:
        file("*") from fastqc_fastp_multiqc.collect()
    
    output:
        file("multiqc_fastqc_fastp.html")

    script:
        """
        multiqc *.zip

        mv multiqc_report.html multiqc_fastqc_fastp.html
	"""
}

/*
 * Process 6: Run multiqc on raw fastq
 */

process multiqc_repair_fastqc {

    label 'multiqc_repair_fastqc'

    conda '/home/ajan/.conda/envs/multiqc'

    publishDir "$multiqc_repair_fastqc_out", mode:'copy'
   
    input:
        file("*") from repair_fastqc_multiqc.collect()
    
    output:
        file("multiqc_repair_fastqc.html")

    script:
        """
        multiqc *.zip

        mv multiqc_report.html multiqc_repair_fastqc.html
	"""
}

workflow.onComplete {
log.info ( workflow.success ? "\n The workflow was complete!" : "Oops .. something went wrong" )
}
