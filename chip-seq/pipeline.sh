#!/usr/bin/env bash

GSE=$1
INPUT_SRR=$2
GENOME=$3
SPECIES=$4
INDEXES=$5
SRRS=${@:6}

INDEX_FILES=$(find ${INDEXES} -name "*.bt2*")
for F in ${INDEX_FILES[@]}; do TAG=${F##*/}; ln -s $F $TAG; done

mkdir ${GSE}

for SRR in ${SRRS}; do :
    echo "Downloading: ${SRR} to ${GSE}"
    rsync -azvvit --partial-dir=.rsync-partial --human-readable --progress rsync://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/${SRR:0:6}/${SRR}/${SRR}.sra ${GSE}
    echo "Fastq-dump: ${GSE}/${SRR}.sra"
    fastq-dump --log-level err --dumpbase --outdir ${GSE} ${GSE}/${SRR}.sra
    echo "Bowtie2: ${GSE}/${SRR}_${GENOME}.sam"
    bowtie2 -p 30 --trim5 5 -S ${GSE}/${SRR}_${GENOME}.sam -x ${GENOME} -U ${GSE}/${SRR}.fastq
    echo "Samtools view: ${GSE}/${SRR}_${GENOME}_not_sorted.bam"
    samtools view -@ 30 -bS -q 10 ${GSE}/${SRR}_${GENOME}.sam -o ${GSE}/${SRR}_${GENOME}_not_sorted.bam
    echo "Samtools sort: ${GSE}/${SRR}_${GENOME}.bam"
    samtools sort -@ 30 ${GSE}/${SRR}_${GENOME}_not_sorted.bam -o ${GSE}/${SRR}_${GENOME}.bam
    echo "Clearing temporary files"
    rm -f ${GSE}/${SRR}.sra ${GSE}/${SRR}.fastq ${GSE}/${SRR}_${GENOME}.sam ${GSE}/${SRR}_${GENOME}_not_sorted.bam
done

echo "Moving ${GSE}/${INPUT_SRR} to ${GSE}/input.bam"
mv ${GSE}/${INPUT_SRR}_${GENOME}.bam ${GSE}/input.bam

cd ${GSE}
INPUT=
echo "Running macs2 narrow fdr=5E-2 in ${GSE} "
for FILE in $(find . -name '*.bam' | sed 's#\./##g' | grep -v 'input')
    do :
        macs2 callpeak --tempdir /tmp -t ${FILE} -c input.bam -f BAM -g ${SPECIES} -n ${FILE%%.bam}_5E-2 -q 5E-2
    done

rm -rf pileup
rm *.r *.summits.bed *.xls *.bw
cd ..