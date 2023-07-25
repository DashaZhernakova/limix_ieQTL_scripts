script_dir="scripts/"

zcat resources/ProbeAnnotation_STARv2.3.0e_Ensembl71.txt.gz | \
    awk 'BEGIN {FS=OFS="\t"}; {if (NR == 1) print "feature_id", "chromosome", "start", "end", "ensembl_gene_id", "gene_name"; else print  $2 ,$4,$5,$6,$2,$3}' \
    > limix_gene_annotation_Ensembl71.txt

c="LLD"
cohorts=("CODAM" "NTR_AFFY" "RS" "LLS_660Q")
script_dir="/groups/umcg-bios/tmp01/users/umcg-dzhernakova/sd-eQTL_scripts/"

for c in ${cohorts[@]}
do
    echo "$c"
    mkdir -p ${c}/data/
    awk -v cohort=${c} '{FS=OFS="\t"}; {if ($1 == cohort) print $0}' samples_per_cohort_geno_expr.txt > ${c}/data/${c}_samples.txt

    #
    # extract samples for the cohort from gene expression table
    #

    python3 ${script_dir}/subset_table_by_gte.py \
        /groups/umcg-bios/tmp01/users/umcg-dzhernakova/data/expression/gene_read_counts_BIOS_and_LLD_passQC.tsv.gz \
        ${c}/data/${c}_samples.txt \
        ${c}/data/${c}.gene_read_counts_BIOS_and_LLD_passQC.tsv \
        2

    # TMM normalize
    ml RPlus/4.0.3-foss-2018c-v21.12.10
    Rscript ${script_dir}/normalizeTMM_cutoffCPM.R \
        ${c}/data/${c}.gene_read_counts_BIOS_and_LLD_passQC.tsv \
        ${c}/data/${c}.gene_read_counts_BIOS_and_LLD_passQC.TMM.txt \
        0.01

    # TPM normalize
    #Rscript ${script_dir}/normalize_TPM_and_rename.R \
    #    ${c}/data/${c}.gene_read_counts_BIOS_and_LLD_passQC.tsv \
    #    resources/Homo_sapiens.GRCh37.75.gene_lengths.txt.gz \
    #    resources/ProbeAnnotation_STARv2.3.0e_Ensembl71.txt.gz


    gzip ${c}/data/${c}.gene_read_counts_BIOS_and_LLD_passQC.tsv
    gzip ${c}/data/${c}.gene_read_counts_BIOS_and_LLD_passQC.TMM.txt
    #gzip ${c}/data/${c}.gene_read_counts_BIOS_and_LLD_passQC.tsv.TPM.txt

    cut -f2,3 ${c}_samples.txt > ${c}_gte.txt


    # genotypes
    plink2 --bfile /groups/umcg-lld/tmp01/other-users/umcg-dzhernakova/genotypes/LLD/LL_combined --maf 0.05 --hwe 1e-6 --make-bed --out tmp_${c}_genotypes_raw
    plink2 --bfile tmp_${c}_genotypes_raw --indep-pairwise 250 50 0.2  --out tmp_${c}_genotypes_pruned
    plink2 --bfile tmp_${c}_genotypes_raw --extract tmp_${c}_genotypes_pruned.prune.in --make-king square --out ${c}

    awk 'BEGIN {FS=OFS="\t"}{ for ( i=1; i<=NF; i++ ) printf $i*2"\t" }{ print"" }' LLD.king > tmp.LLD.king
    cut -f2  ${c}.king.id | tr '\n' '\t'  > ${c}_kinship.txt
    echo "\n" >> ${c}_kinship.txt
    cut -f2  ${c}.king.id | tail -n+2 | paste - <(cat tmp.${c}.king) >> ${c}_kinship.txt

    rm tmp.${c}.king
    rm tmp_${c}_genotypes*
    
done