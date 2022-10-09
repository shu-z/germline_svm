# Germline SVM Classifier

Added preprocessing folder, which needs to be done in the following order:

1. bedpe_combination
2. filtering
3. fuzzy_matching
4. annotation

SVM Classifier to distinguish between somatic and germline SVs when matched normals are not available. Primarily for use with MANTA and SvABA structural variant callers.

### To Run SVM ###

**0. Filtered and Annotated Files** <br /> 
Begin with outputs of preprocessing steps, which: <br />
    * Converted VCF files(from SvABA/Snowman) to BEDPEs
    * Filtered SV calls for high quality SVs
    * Annotated SV Breakpoints for proximity to gNOMAD SVs, LINE/SINE elements
    * Annotated SV Breakpoints for Exon and Whole Gene Impact 

**1. filter_df_newmethod.R** <br />
    This goes through each preprocessed SV file. It first creates a mapping file to match each sample ID to its corresponding tumor type. It then selects samples for training and testing. Finally, SVs greater than 1000bp from these samples are combined into train and test sets. 
    
**2. add_features_testtrainsep.R** <br />
    This converts columns in the SV file into a usable format for the SVM. Additional feature columns are also created. 
    
**3. svm_main_newmethod.R** <br />
    This creates the final somatic/germline training sets. Features are scaled here. The SVM is then run here. Hyperparameters are tuned. The train classifier is then used to predict germline/somatic for the test set. 


Additional analysis files for SVM performance, train/test splitting, misclassification, etc are found under the "analysis" folder. 

One thing to do in the future would be to integrate all these scripts so that the svm and analyses can be run together, and a folder of outputs is made. 
