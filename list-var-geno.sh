#!/bin/bash

PROJECT="/nfs/projects/dbGap/2022-05_All_WGS_WES"
SAMPLE_FILE="$PROJECT/2022-05-12CasesAll.txt"

INCLUDES="--include-gnomad-genome --include-gnomad-exome --include-exac --include-mtr --include-gerp --include-rvis --include-sub-rvis --include-limbr --include-revel --include-trap --include-discovehr --include-known-var --include-primate-ai --include-ccr --include-loftee --include-gnomad-gene-metrics --include-pext --include-mpc --flag-repeat-region --include-syn-rvis --include-genome-asia --include-iranome --include-gme --include-top-med" 

QC="--exclude-artifacts --filter pass,likely,intermediate --exclude-evs-qc-failed --ccds-only --min-coverage 10 --include-qc-missing --qd 5 --qual 50 --mq 40 --gq 20 --snv-fs 60 --indel-fs 200 --snv-sor 3 --indel-sor 10 --rprs -3 --mqrs -10 --het-percent-alt-read 0.3-1"

CODING_SPLICE="HIGH:exon_loss_variant,HIGH:frameshift_variant,HIGH:rare_amino_acid_variant,HIGH:stop_gained,HIGH:start_lost,HIGH:stop_lost,HIGH:splice_acceptor_variant,HIGH:splice_donor_variant,HIGH:gene_fusion,HIGH:bidirectional_gene_fusion,MODERATE:3_prime_UTR_truncation+exon_loss_variant,MODERATE:5_prime_UTR_truncation+exon_loss_variant,MODERATE:coding_sequence_variant,MODERATE:disruptive_inframe_deletion,MODERATE:disruptive_inframe_insertion,MODERATE:conservative_inframe_deletion,MODERATE:conservative_inframe_insertion,MODERATE:missense_variant+splice_region_variant,MODERATE:missense_variant,MODERATE:splice_region_variant,LOW:5_prime_UTR_premature_start_codon_gain_variant,LOW:initiator_codon_variant,LOW:initiator_codon_variant+non_canonical_start_codon,LOW:splice_region_variant+synonymous_variant,LOW:splice_region_variant,LOW:start_retained,LOW:stop_retained_variant,LOW:synonymous_variant"

sh /nfs/goldstein/software/sh/atav.sh --list-var-geno \
--sample ${SAMPLE_FILE} \
$INCLUDES $QC --effect $CODING_SPLICE \
--gnomad-genome-pop global --gnomad-genome-af 0.01 \
--gnomad-exome-pop global --gnomad-exome-af 0.01 \
--exac-pop global --exac-af 0.01 \
--gene $PROJECT/FPF_candidate_genelist.txt \
--out $PROJECT/2022-05_All_WGS-dominantNone


sh /nfs/goldstein/software/sh/atav.sh --list-var-geno \
--sample ${SAMPLE_FILE} \
$INCLUDES $QC \
--effect HIGH:exon_loss_variant,HIGH:frameshift_variant,HIGH:rare_amino_acid_variant,HIGH:stop_gained,HIGH:start_lost,HIGH:stop_lost,HIGH:splice_acceptor_variant,HIGH:splice_donor_variant,HIGH:gene_fusion,HIGH:bidirectional_gene_fusion,MODERATE:3_prime_UTR_truncation+exon_loss_variant,MODERATE:5_prime_UTR_truncation+exon_loss_variant,MODERATE:coding_sequence_variant,MODERATE:disruptive_inframe_deletion,MODERATE:disruptive_inframe_insertion,MODERATE:conservative_inframe_deletion,MODERATE:conservative_inframe_insertion,MODERATE:missense_variant+splice_region_variant,MODERATE:missense_variant \
--polyphen probably --min-primate-ai 0.8 --min-revel-score 0.5 --ensemble-missense-2 \
--min-exac-vqslod-snv -2.632 --min-exac-vqslod-indel 1.262 \
--gnomad-genome-rf-tp-probability-snv 0.01 --gnomad-genome-rf-tp-probability-indel 0.02 \
--gnomad-exome-rf-tp-probability-snv 0.01 --gnomad-exome-rf-tp-probability-indel 0.02 \
--gnomad-genome-pop global --gnomad-genome-af 0.0005 \
--gnomad-exome-pop afr,amr,asj,eas,sas,fin,nfe --gnomad-exome-af 0.0005 \
--exac-pop afr,amr,nfe,fin,eas,sas --exac-af 0.0005 \
--gene $PROJECT/FPF_candidate_genelist.txt \
--out $PROJECT/2022-05_All_WGS-RareEnsemble2
