#! /bin/bash
#$ -b y
#$ -q long
#$ -t 1-22
#$ -N dash
#$ -l m_mem_free=96G
#$ -w e
#$ -e /home/unix/armartin/finrisk/haplotypes/germline/logs
#$ -o /home/unix/armartin/finrisk/haplotypes/germline/logs

MATCH=$1
FAM=$2
OUT=$3

THIS_MATCH=`echo $MATCH | perl -pi -e "s/chr\d+/chr${SGE_TASK_ID}/g"`
THIS_OUT=`echo $OUT | perl -pi -e "s/chr\d+/chr${SGE_TASK_ID}/g"`

if [[ $MATCH == *gz ]];
then
	zcat $THIS_MATCH | cut -f 1,2,4,10,11 | dash_cc $FAM $THIS_OUT
else
	cat $THIS_MATCH | cut -f 1,2,4,10,11 | dash_cc $FAM $THIS_OUT
fi

gzip ${THIS_OUT}.clst
