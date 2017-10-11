#!/bin/bash

# Arg 1: tag to add to analysis
# Arg 2: Read trim length argument
# Arg 3: Expected-errors value for -fastq_filter; settled on 1.5

set -x

TAG=${1:?Must provide analysis tag as first argument}
TrimLength=200
ExpectedErrors=1.5
MinSize=2

DIR=$PWD

SED=sed; which gsed >& /dev/null && SED=gsed

# These are just 3 of the IonTorrent adapters, to see if we can find
# overrepresented sequences there.  FastQC dies if we use a full list
# of adapters.
#
# http://mendel.iontorrent.com/ion-docs/Technical-Note---Transition-from-SFF-to-BAM-format_37421247.html
#
# It seems we don't have to worry about adapter sequences in the BAM
# files produced by TorrentSuite.
#
Adapters=$DIR/candidate-adapter-sequences.txt

# FastQC should be available in the PATH as fastqc
function run_fastqc() {
    for F in ugit_101_*.fastq ; do
        echo "$F ..."
        D=${F%.fastq}
        mkdir -p $D
        fastqc --adapters $Adapters --extract --outdir ${D} $F
    done
}

# run_fastqc



# None of the adapters are overrepresented
# open ugit_101_1.IonXpress_001.2014-07-03/*.html

# get usearch 8.1 from http://www.drive5.com/usearch/download.html
USEARCH=$DIR/usearch8.1.1861_i86osx32
[[ -x $USEARCH ]] || { echo "$USEARCH executable must be present"; exit 1; }

######################################################################
#
# following pipeline at 
#
#      http://drive5.com/usearch/manual/upp_454.html
#
function edgar_pipeline() {
    L=$1
    EE=$2
    tag="usearch_${TAG}_${L}_${EE}_minsize${MinSize}"
    LOG=${tag}.log
    echo "*** Edgar pipeline" > $LOG
    echo "*** http://drive5.com/usearch/manual/upp_454.html" >> $LOG
    echo "***" >> $LOG
    echo >> $LOG
    data_prefix="ugit_101_"
    Sample_prefix="edg_S"
    ALLREADS=allreads.${tag}.fa
    echo "*** Creating $ALLREADS" >> $LOG
    echo "***" >> $LOG
    rm -f $ALLREADS
    for F in ugit_101_*.fastq ; do
        D=${F%.fastq}
        S=${D%%.*}
        S=${Sample_prefix}${S#$data_prefix}
        O=${D}.${tag}.fa
        echo "Trimming $F ..." >> $LOG
        $USEARCH -fastq_filter $F -fastq_maxee ${EE} -fastq_trunclen ${L} -fastaout $O 2>>$LOG
        S=";barcodelabel=$S"
        echo "Adding barcodelabel '$S' ..." >> $LOG
        $SED -e "s/^\(>.\+\)$/\1$S/" $O > ${O}.tmp
        mv -f ${O}.tmp $O
        cat $O >> $ALLREADS
        rm -f $O
    done

    # dereplication
    F=$ALLREADS
    D=${F%.fa}.derep.fa
    echo "*** Dereplication $F to $D" >> $LOG
    echo "***" >> $LOG
    $USEARCH -derep_fulllength $F -fastaout ${D} -sizeout >> $LOG 2>&1

    # OTU clustering
    F=$D
    OUT=otus1.${tag}.fa
    echo "*** OTU clustering $F to $OUT" >> $LOG
    echo "***" >> $LOG
    $USEARCH -cluster_otus $F -minsize $MinSize -otus $OUT -relabel OTU_${tag}_ >> $LOG 2>&1
    OTUS=$OUT

    # map reads to OTUs for quantification
    echo "*** Map reads $ALLREADS to OTUs $OTUS, results to $UC" >> $LOG
    echo "***" >> $LOG
    $USEARCH -usearch_global $ALLREADS -db $OTUS -strand plus -id 0.97 -otutabout otutab.${tag}.txt -biomout otutab.${tag}.json >> $LOG 2>&1

}

edgar_pipeline $TrimLength $ExpectedErrors

