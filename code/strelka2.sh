
#!/bin/bash

module load strelka/2.9.10
source activate strelka-2.9.10

mkdir strelka_out/${NORMAL}-$TUMOR

configureStrelkaSomaticWorkflow.py \
--normalBam ${NORMAL}.fix.markdup.bam \
--tumorBam ${TUMOR}.fix.markdup.bam \
--referenceFastaTAIR10_chr_all.fasta \
--runDir strelka_out/${NORMAL}-$PREFIX

strelka_out/${NORMAL}-$TUMOR/runWorkflow.py -m local -j 8

conda deactivate
