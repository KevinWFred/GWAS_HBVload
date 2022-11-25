# GWAS_HBVload
This project is to find SNPS (imputed from TOPMed) associated with hepatitis B viral load.  
1. Ordinal logistic regression was to estimate the association between SNPs and increasing levels of HBV DNA load.  
2. The HBV DNA load is coded as 1=<300 (Undetectable); 2=300-9999; 3=10000-99999; 4=100000-999999; 5= ≥1 million. This grouping is defined based on a previous publication using the same cohort (https://pubmed.ncbi.nlm.nih.gov/16391218/).
3. Sex, age group at baseline (29-39, 40-49, 50-59, and 60-79 years), and top 5 PCs were used as covariates.
