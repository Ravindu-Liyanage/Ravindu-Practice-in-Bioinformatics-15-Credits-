# Ravindu: Practice-in-Bioinformatics-(15 Credits)
June 10-August 30
## Project Background
Macrophage-dependent immune remodeling and altered lipid metabolism are closely related to cancer progression and metastasis. Aim of the project: To refine and validate a prediction/risk-score model to be used for predicting risk of cancer progression and resistance to treatment. 

Transcriptomic and single-cell cohorts will be analyzed. Macrophage-related differentially expressed genes (MRDGs) will be extracted from GSE databases, combined with TCGA cancer genomic databases consisting of differentially expressed genes (DEGs) and lipid metabolism-related genes (LMRGs). Shared candidates will be screened through univariate Cox analysis and machine-learning procedures to obtain prognostic biomarker algorithms and build a risk-score model. Immune profiling, enrichment analysis, cell-cell communication inference, and pseudotime reconstruction will be conducted. 

The prediction models/algorithms will have great potential to be further developed as predictive biomarkers to predict patient clinical outcome and responsiveness to immunotherapy. This will support treatment decision-making and personalized medicine.
## Datasets Used
- TCGA Pan Cancer Atlas (TCGA-PRAD for prostate adenocarcnoma and TCGA-PAAD for pancreatic adenocarcinoma)
- MSKCC2010 for prostate adenocarcinoma
- Lab produced Novogen May 2026 Mouse RAW sequences

## Types of Analysis Done
### Data Normalization and Transformation
- **VST normalization:** Variance stabilizing transformation was used to reduce the dependence of variance on mean expression and make expression values more suitable for visualization and downstream analyses.
- **logCPM transformation:** Used to stabilize variance and improve comparability across samples before plotting and statistical testing.

### Differential Expression Analysis
- **DESeq2:** Used for differential gene expression analysis on count-based RNA-seq data, with normalization and statistical testing to identify significantly altered genes between groups.
- **edgeR / logCPM normalization:** Used to transform raw count data into log-counts per million for downstream comparison and visualization.
- **limma-voom:** Applied generalized linear modeling and precision-weighted differential expression analysis for robust comparison across sample groups.

### Statistical Testing
- **t-test:** Used to compare expression differences between two groups.
- **One-way ANOVA:** Used to test for expression differences across more than two groups.
- **Kruskal-Wallis test:** Used as a non-parametric alternative when data did not meet normality assumptions.
- **Spearman’s correlation:** Used to assess monotonic relationships between gene expression and other continuous variables.

### Survival Analysis
- **Overall Survival (OS), Disease-Free Survival (DFS), and Progression-Free Survival (PFS):** Used to evaluate the prognostic relevance of candidate genes.
- **Kaplan-Meier analysis:** Generated survival curves to compare patient outcomes between high- and low-expression groups.

### Visualization
- **Heatmaps:** Used to display expression patterns across samples and highlight clustering trends.
- **Boxplots:** Used to compare gene expression distributions between experimental groups.
- **Kaplan-Meier plots:** Used to visualize survival differences between patient subgroups.

### Functional Enrichment Analysis
- **Over-Representation Analysis (ORA):** Used to identify biological processes and pathways enriched in the gene lists.
- **GO enrichment analysis:** Performed to detect overrepresented Gene Ontology terms related to the candidate genes.
- **KEGG enrichment analysis:** Used to identify pathway-level associations and biological mechanisms.
- **Enrichment dot plots and enrichment maps:** Used to visualize enriched terms and the relationships between them.

### Immune Infiltration Analysis
- **MCP-counter:** Used to estimate the abundance of immune and stromal cell populations from transcriptomic data and assess tumor microenvironment composition.

### Gene Alignment and Mapping (Galaxy)
- **Falco-tool:** A high-speed emulation of the popular FastQC software for the Quality Control of sequencing data.
- **MultiQC tool:** Aggregates results from multiple bioinformatics analyses across many samples into a single report.
- **CutAdapter tool:** Designed to removed adapter sequences, primers, poly-A-tails and other unwanted sequences from high-throughput sequencing reads.
- **RNA-Star tool:** An aligner for RNA-seq data mapping using a strategy to account for spliced alignments.
- **FeatureCounts tool:** A lightweight read counting program to measure gene expression in RNA-seq experiments from SAM or BAM files.
- **ColumnJoin tool:** Join the columns on multiple databases.

## LOGS
### Week 1 (June 11–12)

#### June 11
Downloaded RNA-seq data files from the GDC data portal for TCGA-PRAD (prostate cancer) directly through R by using TCGAbiolinks package from Bioconductor. Then I used this data to conduct my expression analysis of 3 interested genes in prostate adenocarcinoma. Since TCGA cohorts main samples types were primary tumor and normal tissue conditions, I conducted my differential analysis between these two groups. Faced many issues after normalizing and variance stabilizing the data at first. The issue it seemed at the end was when creating the expression data for just the 3 interested genes. I was finally able to generate all boxplots needed and interpret the results. 

#### June 12
I ran an analysis to see if there was a significant difference between the expression of these 3 interested genes across their Gleason scores. I used the barcode of each sample for Gleason categorization, where 6 was low, 7 was medium and 9 or higher was considered high group. I then ran ANOVA statistical testing to check if these differences were statistically significant. I also carried out the survival study and created the KM plots. Since I did not contain DFS data in the meta data I downloaded, I focused on basing my survival study on Overall Survival (OS). I also plotted the heatmaps for each of the 3 interested genes, seeing their coexpression with other genes across the RNA-seq data. I was able to conclude my first analysis.

### Week 2 (June 15–18)

#### June 15
I created a report on my first analysis, interpreting the results and arriving at conclusions from the plots and heatmaps. Also compared my results with another colleague for any potential errors. I also had a pretty long and productive lab meeting with Jenny and the others in the lab, catching up on each individual's current work. At the end of the day, I also researched into AI platforms that could segment PSMA PET/CT scans online, as well as local models, which was an analysis planned for later on in the course and sent a list to Jenny.

#### June 16
Got access to my own profile and folder in the server and shifted all my workflow into the computer in the lab. Downloaded and installed all required software. I tried running the same scripts again in the new environment and faced issues with accessing Bioconductor packages for the whole day. So I shifted my attention further into the local AI models that can be accessed from GitHub. I specifically focused on a model named nnUNet. I went through the starter guide and set up directories in my local folder for pre-processing steps of the data. The next step is to receive PSMA PET/CT scans from Jenny to try out the code and train the model on it. 

#### June 17
Was finally able to get Bioconductor working again and installed all necessary packages. I redid the entire TCGA-PRAD 3-gene workflow, this time using an RDS file that contained the whole TCGA pan-cohort data from all types of cancers. Explored the dataset further for future analysis. 

#### June 18
I started to work on the same analysis as before, but now with the MSKCC2010 prostate adenocarcinoma dataset. I got access to the RDS file. I carried out all the steps and was able to replicate the results achieved by Martina. 

### Week 3 (June 22–26th)

#### June 22
Carried out expression analysis using generalized linear modeling and enrichment analysis on the DEGs of the TCGA-PRAD dataset. Plotted dot plots and enrichment maps in R. I attempted the same analysis with the MSKCC2010 dataset but the pathways enriched were not significant due to the low gene ratio in each pathway, so I didn't write a report on that. Also started an analysis to compare FOLH1 expression levels (responsible for PSMA) across different cancer types within the TCGA pan-cancer atlas dataset.

I also received a new task from Jenny, this time looking into TCGA-PAAD pancreatic adenocarcinoma, related to an ongoing manuscript in the lab undergoing revision. My task was to replicate and validate results obtained by Amjad, another bioinformatician. First carried out correlation testing for 31 genes of interest within the dataset against expression of PIP5K1A, and was able to replicate the results.  

#### June 23
Finished the FOLH1 expression comparison across all cancer types and obtained top expressed cancer types with statistical testing. Plotted barplots, dotplots, boxplots as well. Then I moved on to the pancreatic adenocarcinoma workflow where I conducted PFS survival analysis using the RDS on my computer locally. However, all PFS columns as well as the ones from which I could have derived PFS values were empty (NA). So I chose to use an online plotting tool, TIMER3, by compbio. I was able to generate KM plots for each of the immune infiltrate categories, and found two with significantly increased risk. All plots and statistics were saved and a report was made. I used the same TIMER3 tool to correlate PIP5K1A expression with each of the immune infiltrates. 

#### June 24
Received a second pancreatic adenocarcinoma request from Jenny, this time to look into PIP5K1A gene DFS association with classical and basal-like subtypes of cancer. In order to get DFS clinical data, I had to download a separate TSV file from cBioPortal. Gene expression was divided into High and Low groups based on the median values and KM curves were plotted followed by multiple methods of statistical testing. A report on this was created and sent to Jenny. 

#### June 25
Carried out a very similar analysis to the previous one, but this time looking into OS data. Had a meeting with Jenny to discuss my results and was advised to reduce DFS months when one group’s samples exceed the other group’s by a lot. Switching to the prostate adenocarcinoma, I did a full immune gene set analysis by using immune scores, enrichment and correlation of these immune genes with PIP5K1A expression. I then created a report.

#### June 26
Started going over all my TCGA-PRAD analyses to catch any mistakes and compare results with Martina. Next, did a new analysis, looking into immune categorization of TCGA-PRAD primary tumors. Had the weekly meeting with Jenny and others from the lab to quickly catch up, followed by another meeting, this time with a representative from a scRNA sequencing company, Parse, regarding their technology and how we could use it in our lab, which was very interesting. I also did Jenny’s request based on the feedback and cut off longer DFS months in the TCGA-PAAD basal-like and classical-like subtype survival analysis.

### Week 4 (June 29- July 3rd)

#### June 29
Further worked on the full prostate adenocarcinoma report, adding interpretations wherever missing. Carried out the same immune profiling as before, but this time looking into metastatic samples from MSKCC2020 dataset. Created a short report based on it. I also checked correlated gene expression with FOLH1 in primary tumors and normal tissues in the TCGA-PRAD dataset. Also carried out an analysis where I checked FOLH1 expression in primary tumors with known biomarkers of prostate adenocarcinoma. 

#### June 30
In TCGA-PRAD, carried out differential expression analysis between high and low FOLH1 expression groups and did enrichment analysis to identify top enriched pathways. I next carried out immune profiling, but this time, seeing which immune pathways are enriched with respect to high and low FOLH1 expression. Also, I did Microenvironment cell populations (MCP) counter analysis for TCGA-PRAD primary tumors with respect to FOLH1 expression levels. This was added as an extension to the previous script. Lastly, I compared previous scripts with Martina’s scripts and found some errors. Redid some of the works, like the Gleason categorization as well as the FOLH1 expression across various cancers analysis. 

#### July 1
Redid analysis where there were some errors, like the different cancer types analysis. Continued to double check everything, re-run all analyses and write more biological interpretations in the report, wherever it felt lacking. Later I received another request from Nattawan to restyle some of my pancreatic adenocarcinoma KM plots to match the style of a paper, as well as attempt to recreate results obtained from an online KM plotting tool. However, I ran into issues when doing that, since that tool used a combination of 16 databases, 12 from GEO and 4 other cohorts that I did not get access to. Therefore, I focused on re-styling just my previous plots, specifically the basal and classical-like subtypes DFS analysis.

#### July 2
Re-did the KM plots for Basal-like and classical-like subtypes with relation to high and low PIP5K1A for the TCGA-PAAD analysis as per Nattawan's request. There was a small mixup with the colours of the two curves so i rectified it. I completely re-structured my GitHub repository so that its fully up to dat and would reflect my actual work environment. Also expereinced some wet lab when Nattawan and Parsa, the two post-docs in the lab, carried out western blot on protein extracted from mice samples with knock-out and wild types for PIP5K1A. 

### Week 5 (July 6-July 10)

#### July 6
I redid the T-stage FOLH1 analaysis to not only include primary tumours but all types since the independant variable had to be the T stage, and then sent in my full report to Jenny before todays lab meeting. Afterwords, had the weekly lab meeting with everyone in the lab to discuss about what I've done so far, and scheduled a meeting for Wednesday. Afterwords, i double checked my pancreatic cancer workflow (TCGA-PAAD) from which two plots i generated will be used in an upcoming publication from the lab, to make sure thaty there are indeed only 150 samples of andenocarcinoma in the dataset, and was abl to confirm it by reading literature and checking data accessed in different methods. 

#### July 7
Updated my README file and pushed to the repository. Worked in the wet lab a bit by carrying out the gel transfer step of western blot under Nattawan's supervision. additonally did the primary antibody coating step and then left the membrane samples overnight at 4 degrees celcius. For the rest of the, i cleaned all of my R scripts to stick to a specific format and reuploaded everything to the git repository.

#### July 8
Continued the wet lab section. I coated the membrane sampled with the secondary antibody and then later looked at the western blot results under the chemiluminescent marker. Later i had the scheduled meeting with Jenny and Martina to discuss our results from the analysis carried out so far with regard to Pca. we were then joined by another senior professor, Felecia M, who provided advice and directions to focus on based on her expertise and experience in the field. We got numerous more biomarkers to focus our analysis on. After that, we discussed about the next analysis that i will move on to, where i will be using our own lab produced raw data from mice samples. I received access to the RAW FASTA files and began the alignment steps using the online platform Galaxy. The files were left to upload overnight. **

#### July 9
Continued uploading sequencing files to Galaxy for the allignment, trimming and mapping steps. Wrote a word document outlining my plans for the new datasets, as well as plans on how to improve my previous Prostate adenocarcinoma analyses with Felicia's new advice and bio-markers.Carried out the coreelation of FOLH1 with the new biomarkers and was able to find all bio-markers within all samples. Also, received information from Nattawan regarding the sequence data and structuring from the mice samples i would be allinging. I started making paried collections on Galaxy based on model of the mice, mostly based on treatement type. Started the quality control step, seperating based on each model, starting with "orginal tumor." Used the flatten tool as well as the falco tool.

#### July 10
Had weekly meeting with Jenny and discussed about my progress as well as my results from the new biomarkers analysis. I then wrote the report on that. Went further on with the Quality control steps with the raw sequence files on Galazy. Faced some issues on certain fasta files during the falco tool step. The error seems to either be during the flattening step, or soemthing inherently wrong with he galaxy system, and not corrupt files from my side. 

### Week 8 (July 24-July 31)

#### July 24
Was able to get all treatement groups to pass through the falco tooland check their QCs. Did hands on trimming on paried collections with Quality cutoff (R1) and min length (R1) both at 20. I faced a storage issues with running out of space on the galaxy server, so i decided to delete all my files and start over again and redo the whole analysis, this time on the temporary scratch storage that deletes files in 60 days. I also carried out the new biomarker analysis on the metastatic samples from the MSKCC2010 cohort, creates the report and sent it to Martina

#### July 27
Continued with the alignment steps on galaxy before it unexpectedly shut down for most of the day, before coming back up. Received Martina's analysis and started comparing our results for the primary tumor biomarkers. Afterwords, in Galaxy, i started facing a recurring issue with the GTF file during the mapping stage, therefore i emailed some other bioinformatian from the team regarding this issue.

#### July 28
Today I followed Josephine throughout the lab as she prepared for and carried out flow cytometry on blood samples she recieved from a healthy patient from the hospital. First, we separated the blood components and extracted only the PBMCs, which were later spiked with prostate cancer cells from our cell line. Later observed how the flow cytometer works and how it outputs results. Martuza came to the lab and i was able to get my issues with the mapping stage on the Galaxy platform resolved and so i continued mapping all of my sequences collections using the RNA STAR tool followed by the FeatureCounts tool in order to geenrate counts for each gene, which is where i ended up uplaoding the GTF file i ahd issues with previously. Martuza helped me find the appropriate GTF file needed for this. Was able to combine and create count tables for all of the treatement types, thereby finsihing the alignment stage of the analysis, so that i can move on to the downstream analyses as adviced by Jenny.   

#### July 29
uploaded my alignments (counts) to R locally and explored their structures. Had a one on one meeting with Jenny regarding my next steps. Received guidance and created a plan for the next couple of weeks. Was advised to work towards creating my own manuscript. 

#### July 30
Started doing the analysis of interested biomarkers on the Novogene sequences on our own mice models . I focused specifically on the original tumor sample collection first. looked at the expression differences between WT and KO. individual correlation plots with PSMA for each biomarker was also created. Generated results and plots but didnt analyse the results yet. Additionally i also carried out the polishing out of some of my previous KM plots to fit into the standard publishing read format.

#### July 31
Had the weekly meeting with Jenny and discussed my progress and current focuses in great detail. understood that i should have started my analysis with PCA clustering first to see if the samples show the expected trends and seperation. So i carried that out with the original tumor samples collection. The trends at first seemed to be wonky, especially given that a few WT samples in this collection had low assigned percentages after QC. Therefore i removed these samples from the PCA to see if there're patterns then. i then continued carrying out the same clustering analysis as well as the previously done biomarker analysis for each collection of treatment methods, also including the expression differences. was able to generate plots and results for each correlation. 

### Week 10 (August 3-August 7)

#### August 3
Carried out further analysis with the mice model data. Nattawan came back to the lab after vacation and we discussed about the samples with low assignmed percentages as well as the werid samples with 3 replicate count files. Nattawan sent an email to Novogene sequencing askign for clarification on why it is as so. I started doing basic analyses for each treatment category, looking into DEGs, clustering based on all significant DEGs, top 100 and top50 DEGs and creating heatmaps for each. I also did pathway enruchment analysis using GO BP, GSEA, HALLMARK as well as immunesigDB pathways. Later, i continued to make a powerpoint with my results to show Jenny in our next meeting.

#### August 4
Nattawan forwarded me the reply from the sequencing company regarding why one sample had 3 RAW files. they exlained the reasoning and suggested merging the 3 files for downstream analysis. So i did that locally on R by creating a new column in my raw counts matrix, then removing the original three, before normalising or VST the counts. Furthermore i realised my DEG should be more meaningfully done looking into both treatement control or not as well as type of KO, and compare within these groups. So i rewrote all the scripts to take that into account. For ISA treateemnt specifically, the study was designed in a way to also see if there was a difference btween genders (sexes). So i also considered sex wise comparisson of DEGS in just that collection.  
