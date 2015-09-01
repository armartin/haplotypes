Haplotypes
=================
Overall description: Perform haplotype calling, run QC across various metrics, convert files, perform demographic inference, perform IBD mapping

## Pipeline Map ##
#### 1.) Phase #####
##### 2.) Infer shared haplotypes #####
##### 3.) Perform haplotype QC #####
##### 4.) Infer demographic history from haplotype sharing  #####
##### 5.) Perform IBD mapping for traits of interest #####

## 1.) Phase ###
A number of phasing options are available, for example BEAGLE, SHAPEIT, and HAPI-UR. Beagle has been around the longest and was considered the gold standard for improvement. SHAPEIT made substantial computational and accuracy gains for intermediate sample sizes (~100s to 1000s). HAPI-UR is more accurate than SHAPEIT on larger datasets (5000s-10000s). After phasing the data, convert it so that haplotypes can be inferred.

#### HAPI-UR running notes ####
* Indels break it. Ref/alt alleles must have a length no longer than 1 character
* Add genetic position for more reliable inference. If first genetic position is 0.0, it will try to set genetic from physical positions. Instead, set it to a very small number (e.g. 1e-5).

#### Convert phased haplotypes to format for inferring shared IBD ####
Example from HAPI-UR to RefinedIBD for haplotype inference:
```python hapiur2vcf.py \
--phgeno dbs_1-10_pchip-qc_hardintersect_chr22_3x.phgeno \
--phind dbs_1-10_pchip-qc_hardintersect_chr22_3x.phind \
--phsnp dbs_1-10_pchip-qc_hardintersect_chr22_3x.phsnp \
--out dbs_1-10_pchip-qc_hardintersect_chr22_3x.vcf.gz```


## 2.) Infer shared haplotypes ###
There are several software packages available for calling shared haplotypes, including e.g. GERMLINE and BEAGLE (i.e. RefinedIBD)
