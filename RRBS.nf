

// TODO:
// - ADD MULTIQC TO STAR ALIGNMENT 
// - MAKE THE CODE EASY TO ACCESS (universal path variable, length in fastp, etc.)
// - COMMENT EVERYTHING

/*
 * Define output folders
 */ 

head_dir = '/archive/users/ajan/MOBIT/mobit_rrbs'

//fastq_in = '/archive/users/ajan/MOBIT/mobit_rrbs/test_samples/10k'
fastq_in = '/archive/projects/MOBIT/RRBS/fasta/'
bismark_genome = '/archive/users/ajan/MOBIT/mobit_rrbs/reference_files'

fastqc_raw_out = '/archive/users/ajan/MOBIT/mobit_rrbs/fastqc_raw'
fastp_out = '/archive/users/ajan/MOBIT/mobit_rrbs/fastp'
fastqc_fastp_out = '/archive/users/ajan/MOBIT/mobit_rrbs/fastp_fastqc'
bismark_alignments_out = '/archive/users/ajan/MOBIT/mobit_rrbs/bismark'
bismark_extractions_out = '/archive/users/ajan/MOBIT/mobit_rrbs/bismark'
qualimap_out = '/archive/users/ajan/MOBIT/mobit_rrbs/qualimap'
preseq_out = '/archive/users/ajan/MOBIT/mobit_rrbs/preseq'

multiqc_fastp_out = '/archive/users/ajan/MOBIT/mobit_rrbs/multiqc/fastp'
multiqc_fastqc_raw_out = '/archive/users/ajan/MOBIT/mobit_rrbs/multiqc/fastqc_raw'
multiqc_fastqc_fastp_out = '/archive/users/ajan/MOBIT/mobit_rrbs/multiqc/fastqc_fastp'
multiqc_bismark_out = '/archive/users/ajan/MOBIT/mobit_rrbs/multiqc/bismark_align'
multiqc_bismark_extractor_out = '/archive/users/ajan/MOBIT/mobit_rrbs/multiqc/multiqc_bismark_extractor_out'
multiqc_qualimap_out = '/archive/users/ajan/MOBIT/mobit_rrbs/multiqc/qualimap'
multiqc_preseq_out = '/archive/users/ajan/MOBIT/mobit_rrbs/multiqc/preseq'

bismark2report_alignment_out = '/archive/users/ajan/MOBIT/mobit_rrbs/bismark_reports'

log.info """\


====================================================================================================================

 Methyl-seq - Development v.0.1a 

=====================================================================================================================
|                                                                                                                   
| .fastq files                           : ${fastq_in}                                                                                                                          
| Head output dir                        : ${head_dir}                             
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
 * Process 2: fastp preprocessing
 * For RRBS and WGS fastp instead of Trimmomatic
 * https://robertslab.github.io/sams-notebook/2020/03/06/TrimmingMultiQC-Methcompare-Bisulfite-FastQs-with-fastp-on-Mox.html
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
        set val(id), file("${id}_1_fastp.fastq.gz"), file("${id}_2_fastp.fastq.gz") into (fastp_fastqc, fastp_bismark_align)
        file("${id}_fastp.json") into fastp_multiqc
        file "${id}_fastp.html"

    script:
        """
	fastp -i ${fastq_1} -I ${fastq_2} -o ${id}_1_fastp.fastq.gz -O ${id}_2_fastp.fastq.gz -j ${id}_fastp.json -h ${id}_fastp.html \
    -q 30 \
    -e 25 \
    -n 3 \
    -l 61 \
    -c \
    -x \
    -p \
    --detect_adapter_for_pe \
    --trim_tail1 2 \
    --trim_front2 2 \
    -w ${task.cpus}

	"""
}  

 /*
  * Process 3: Fastqc on preprocessed reads
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
  * Process 4: Align with bismark
  */

 process bismark_align {

     label 'bismark_align'

     module 'MOD/bismark/0.23.1'

     publishDir "$bismark_alignments_out", mode:'copy'
   
     input:
         set val(id), file(fastq_1), file(fastq_2) from fastp_bismark_align
    
     output:
         set val(id), file ("${id}*.bam") into (bismark_align_bismark_extract)
         path ("${id}*") into bismark_alignment_paths
         set val(id), file ("${id}*report.txt") into (bismark_aligned_path_report)
         set val(id), file ("${id}*nucleotide_stats.txt") into (bismark_aligned_path_nucleotide)
         set val(id), file ("${id}*bam") into (bismark_aligned_path_bam)


     script:
         """
         bismark --parallel ${task.cpus} --gzip -n 1 --nucleotide_coverage --genome ${bismark_genome} -1 ${fastq_1} -2 ${fastq_2}
 	"""
 }

   /*
  * Process 5: Extract methylation calls 
  */

 process bismark_methylation_extractor {

     label 'bismark_extract'

     module 'MOD/bismark/0.23.1'

     publishDir "$bismark_extractions_out", mode:'copy'
   
     input:
         set val(id), file(bam) from bismark_align_bismark_extract
    
     output:
         set val(id), file ("${id}*txt") into (bismark_meth_extracted)
         path ("${id}*") into bismark_meth_extracted_path
         set val(id), file ("${id}*report.txt") into (bismark_meth_extracted_path_report)
         set val(id), file ("${id}*M-bias.txt") into (bismark_mbias_path_report)

     script:
         """
         bismark_methylation_extractor --gzip -p --no_overlap --bedGraph ${bam}
 	"""
 }

 /* 
 * Process 6: Create bismark2reports
 */
process bismark2report {

    label 'bismark2report_alignment'

    module 'MOD/bismark/0.23.1'

    publishDir "$bismark2report_alignment_out", mode:'copy'
   
    input:
        set val(id), file(alignment) from bismark_aligned_path_report
        set val(id), file(methylation) from bismark_meth_extracted_path_report
        set val(id), file(mbias) from bismark_mbias_path_report
        set val(id), file(nuc) from bismark_aligned_path_nucleotide

    
    output:
        file("*") into bismark2report_path_out
        //path id into bismark2report_path_out

    script:
        """
        bismark2report --alignment_report ${alignment} --splitting_report ${methylation} --mbias_report ${mbias} --nucleotide_report ${nuc}
	"""
}

/* 
 * Process 7: Summarize bismark alignemnt and methylation calls
 */
process bismark2summary {

    label 'bismark2summary'

    module 'MOD/bismark/0.23.1'
   
    input:
        set val(id), file(alignment) from bismark2report_path_out.collect()

    script:
        """
        cd ${bismark_extractions_out}
        bismark2summary
        mv bismark_summary_report.html ${head_dir}/bismark_summary_report.html
        mv bismark_summary_report.txt ${head_dir}/bismark_summary_report.txt
	"""
}

/* 
 * Process 8: Sort bam file for qualimap QC
 */
 process samtools_sort {

     label 'samtools_sort'

     conda '/home/ajan/.conda/envs/samtools'

     // publishDir "$samtools_out", mode:'copy'
   
     input:
         set val(id), file(bam) from bismark_aligned_path_bam
    
     output:
         set val(id), file("${id}_sorted.bam") into (sorted_bam, sorted_bam2)

     script:
     """
     samtools sort --threads ${task.cpus} -o ${id}_sorted.bam ${bam}
	"""
 }

/* 
 * Process 9: Qualimap QC
 */
 process qualimap {

     label 'qualimap'

     conda '/home/ajan/.conda/envs/qualimap'

     publishDir "$qualimap_out", mode:'copy'
   
     input:
         set val(id), file(bam) from sorted_bam
    
     output:
         path id into qualimap_out_path

     script:
         """
         export JAVA_OPTS="-Djava.io.tmpdir=/archive/users/ajan/MOBIT/mobit_rrbs/trash"

         qualimap bamqc -bam ${bam} -outdir ${id}/ -c -p strand-specific-reverse -outformat PDF:HTML --java-mem-size=16G
	"""
 }


/* 
 * Process 10: preseq QC
 */
 process preseq {

     label 'preseq'

     conda '/home/ajan/.conda/envs/preseq'

     publishDir "$preseq_out", mode:'copy'
   
     input:
         set val(id), file(bam) from sorted_bam2
    
     output:
         set val(id), file("${id}_ccurve.txt") into preseq_out_file

     script:
         """
         preseq c_curve -output ${id}_ccurve.txt -verbose -pe -bam -seed 42 ${bam}
	"""
 }

  /*
 * Process 11: Run multiqc on preseq
 */

process multiqc_preseq {

    label 'multiqc_preseq'

    conda '/home/ajan/.conda/envs/multiqc'

    publishDir "$multiqc_preseq_out", mode:'copy'
   
    input:
        file("*") from preseq_out_file.collect()
    
    output:
        file("multiqc_preseq.html")

    script:
        """
        multiqc *txt -n multiqc_preseq
	"""
}



 /*
 * Process 12: Run multiqc on qualimap
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
 * Process 13: Run multiqc on bismark alignments
 */
process multiqc_bismark_extractor {

    label 'bismark_extractor'

    conda '/home/ajan/.conda/envs/multiqc'

    publishDir "$multiqc_bismark_extractor_out", mode:'copy'
   
    input:
        file("*") from bismark_meth_extracted_path.collect()
    
    output:
        file("bismark_extractor.html")

    script:
        """
        multiqc * -n bismark_extractor --interactive
	"""
}

/* 
 * Process 14: Run multiqc on bismark alignments
 */
process multiqc_bismark {

    label 'multiqc_bismark'

    conda '/home/ajan/.conda/envs/multiqc'

    publishDir "$multiqc_bismark_out", mode:'copy'
   
    input:
        file("*") from bismark_alignment_paths.collect()
    
    output:
        file("multiqc_bismark.html")

    script:
        """
        multiqc * -n multiqc_bismark --interactive
	"""
}


/* 
 * Process 15: Run multiqc on fastp filtered fastq
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
        multiqc *.json -n multiqc_fastp --interactive
	"""
}

/*
 * Process 16: Run multiqc on raw fastq
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

        mv multiqc_report.html multiqc_fastqc_raw.html --interactive
	"""
}

/*
 * Process 17: Run multiqc on fastp filtered FastQC
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

        mv multiqc_report.html multiqc_fastqc_fastp.html --interactive
	"""
}

workflow.onComplete {
log.info ( workflow.success ? "\n The workflow was complete!" : "Oops .. something went wrong" )
}
