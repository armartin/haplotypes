#! /bin/bash
#$ -b y
#$ -q long
#$ -t 1-22
#$ -N dash
#$ -l m_mem_free=64G
#$ -w e
#$ -e /home/unix/armartin/finrisk/haplotypes/germline/logs
#$ -o /home/unix/armartin/finrisk/haplotypes/germline/logs

MATCH=$1
FAM=$2
OUT=$3

if [[ $MATCH == *gz ]];
then
	zcat $MATCH | \
	cut -f 1,2,4,10,11 | \
	dash_cc $FAM
	$OUT
else
	cat $MATCH | \
	cut -f 1,2,4,10,11 | \
	dash_cc $FAM \
	$OUT
fi

gzip ${OUT}.clst
