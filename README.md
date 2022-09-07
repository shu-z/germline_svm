# Germline SVM Classifier

Added preprocessing folder, which needs to be done in the following order:

1. bedpe_combination
2. filtering
3. fuzzy_matching
4. annotation

SVM Classifier to distinguish between somatic and germline SVs when matched normals are not available. Primarily for use with MANTA and SvABA structural variant callers.

**To Run SVM**

0. Filtered and Annotated Files
    - VCF (from SvABA/Snowman) to Bedpe conversion 
    - Filtered for high quality SVs
    - SV Breakpoints annotated for proximity to gNOMAD SVs, LINE/SINE elements
    - SV Breakpoints annotated for Exon and Whole Gene Impact 

1. filter_df_newmethod.R
2. add_features_testtrainsep.R
3. svm_main_newmethod.R
