# Pipeline for  pre-miRNAs discovery in the SARS-CoV-2 virus

This readme contains a description of the steps to obtain the main results of:

G.A. Merino, J. Raad, L.A. Bugnon, C. Yones, L. Kamenetzky, J. Claus, F. Ariel, D.H. Milone, G. Stegmayer, “**Novel SARS-CoV-2 encoded small RNAs in the passage to humans**,” *Bioinformatics*, 2020, [DOI](https://doi.org/10.1093/bioinformatics/btaa1002).

The following sections contain the steps and the libraries used in each of them. All datasets and libraries are open sourced and free to use. If you use any of the following in your pipelines, please cite them properly.

![Abstract](abstract.png)

## 1. Layout

Several intermediate results are provided in this repository, thus you don’t need to run the whole pipeline. 

```
sarscov2-mirna-discovery
│
├── genomes
│   └── NC_045512.2_Wuhan-Hu-1.fasta  -> The complete SARS-CoV2 isolate Wuhan-Hu-1 genome. 
│
├── sequences
│   ├── sars-cov2_hairpins.fasta      -> Hairpin sequences extracted from NC_045512.2_Wuhan-Hu-1.
│   ├── sars-cov2_hairpins.fold       -> sars-cov2_hairpins.fasta and its folding structure.
│   ├── pre-miRNAs_virus.fasta        -> Hairpin sequences extracted from known virus pre-miRNAs.
│   ├── pre-miRNAs_virus.fold         -> pre-miRNAs_virus.fasta and its folding structure. 
│   └── unlabeled_hairpins.fold       -> Hairpin sequences and its folding structure extracted from 
|                                           the human genome, sampled from sequences that do not code 
|   					    pre-miRNAs. Provided by external storage.
│
├── pre-miRNAs_features
│   ├── sars-cov2_hairpins.csv        -> Features extracted from sars-cov2_hairpins.fasta.
│   ├── pre-miRNAs_virus.csv          -> Features extracted from pre-miRNAs_virus.fasta.
│   └── unlabeled_hairpins.csv        -> Features extracted from unlabeled_hairpins.fasta. 
|					     Provided by external storage.
│
├── pre-miRNAs_models                 -> Trained models for pre-miRNA prediction.
│
├── pre-miRNAs_predictions            -> Predictions on SARS-CoV2 hairpin sequences by each model.
│
├── src
|   ├── train_pre-miRNA_models.ipynb  -> Source code to train the pre-miRNAs prediction models.
|   ├── predict_pre-miRNAs.ipynb      -> Source code to predict pre-miRNAs on SARS-CoV2 sequences.
|   ├── DifferentialExpression.Rmd    -> Notebook for differential expression analysis.
|   ├── deregulatedTargets.sh
|   ├── DETargetsExploration.Rmd
|   ├── DOWNDETargetsExploration.Rmd
|   ├── extractDianaIDs.sh
|   ├── joinTargets.sh
|   ├── mappingIDs.Rmd
|   ├── ORADEGenes.Rmd
|   ├── ORADETargets.Rmd
|   ├── targetsExploration.Rmd
|   └── link_Figs.R                   -> Notebook to generate manuscript figures
|
├── matures                           
|   └── miRNAs.csv                    -> Mature miRNAs sequences found. 
|
├── targets                           -> Target prediction files
|   ├── miRDB
|   ├── Diana
|   ├── TargetsDE
|   └── overlap_<mirna>.tab            -> Target overlap for each miRNA 
|
├── expression_data                    -> Gene expression data
|
└── functional_enrichment              -> Functional enrichment results
```

##  2. Data preparation

Note that all the results generated in this section are provided in  [`sequences/`](sequences)  and  [`pre-miRNAs_features/`](pre-miRNAs_features) directories.

### 2.1. Download complete genomes

SARS-CoV-2 genome was obtained from the [NCBI GenBank](https://www.ncbi.nlm.nih.gov/nuccore/1798174254) on March 2020.

These ".fasta" files can be also found in this repository, in the [`genomes/`](genomes) directory.

### 2.2. Extract hairpin-like sequences

The SARS-CoV-2 genome fasta file is cut in overlapping windows using the [Hextractor R package](https://cran.r-project.org/web/packages/HextractoR/index.html) [[2]](#ref2).

You need to have the following software installed in your system:
- [R](https://www.r-project.org/) 
- [RNAfold](https://www.tbi.univie.ac.at/RNA/)
- [BLAST](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)

In a R console, install and load the package with:
```R
> install.packages("snow")
> install.packages("HextractoR")
> library(HextractoR)
```

Then, run the main script with the following parameters:    

```R
> HextractoR("genomes/NC_045512.2_Wuhan-Hu-1.fasta", min_valid_nucleotides = 500, window_size = 600, window_step = 100, only_sloop = T, min_length = 60, min_bp = 16, margin_bp = 6, blast_evalue = 0.005, identity_threshold = 90, nthreads = 4, nworks = 4, filter_files = { })
```

This script generates a file with several hairpin-like sequences and its corresponding folding structure prediction. 

To model the positive class, well known pre-miRNAs from viruses hairpins were extracted from [mirbase](http://www.mirbase.org/). These sequences were 	folded with RNAfold using the following commands:

```bash
RNAfold --noPS --infile=sequences/pre-miRNAs_virus.fasta --outfile=sequences/pre-miRNAs_virus.fold
RNAfold --noPS --infile=sequences/sars-cov2_hairpins.fasta --outfile=sequences/sars-cov2_hairpins.fold
```

Additional sequences were used to train the mirDNN. A set of 1M hairpin-like sequences from the human genome, which are not pre-miRNAs, were used to model the negative set. The folding structure (unlabeled_hairpins.fold) and features (unlabeled_hairpins.csv) for these sequences can be downloaded from [this external repository](https://sourceforge.net/projects/sourcesinc/files/mirdata/sequences/unlabeled.tar.gz). You will only need to download the unlabeled_hairpins.fold if you want to train the mirDNN. This data and further details are available in [[3]](#ref3).

At this point, the sequence and folding prediction for each hairpin-like sequence in SARS-CoV-2, the known virus pre-miRNAs and some human non-pre-miRNA sequences are ready.

### 2.3. Hairpin features extraction

In order to extract meaningful features from the sequences, the [miRNAfe package](http://sourceforge.net/projects/sourcesinc/files/mirnafe/0.90/) was used [[4]](#ref4). Once installed, you can run the feature extraction with:

```matlab
miRNAfe('sars-cov2_hairpins.fasta', 'config/Example_prediction.yaml');
miRNAfe('pre-miRNAs_virus.fasta', 'config/Example_prediction.yaml');
```
The ".yaml" file is provided with the package.  

## 3. Training pre-miRNAs prediction models

The [training notebook](src/train_pre-miRNA_models.ipynb) is provided with instructions to train the OC-SVM [[5]](#ref5), deeSOM [[6]](#ref6) and mirDNN [[7]](#ref7) classifiers. Doing so may take several hours. 

If you want to directly predict the sequences, you can use the trained models that are provided in [`pre-miRNAs_models/`](pre-miRNAs_models) as explained in the following section.

## 4. Finding pre-miRNAs candidates in SARS-CoV-2

The [prediction notebook](src/predict_pre-miRNAs.ipynb) is provided with instructions to use the trained models to rank the SARS-CoV-2 sequences (a high rank means a higher chance of being a pre-miRNA). 

## 5. Predicting miRNAs targets

Once mature miRNAs were identified by combining MatureBayes predictions and the small RNA-seq reads profiles, their sequences (provided in  [`matures/miRNAs.csv`](matures/miRNAs.csv)) must be  submitted to [Diana MR MicroT](http://diana.imis.athena-innovation.gr/DianaTools/index.php?r=mrmicrot/index) and [miRDB (Custom prediction)](http://www.mirdb.org/custom.html) for predicting human gene targets. 

The prediction files from miRDB and Diana MR Micro T, one per SARS-CoV-2 miRNA, are provided in the [`targets/miRDB/`](targets/miRDB/) and [`targets/Diana/`](targets/Diana) directories, respectively. Once you have downloaded them, you can run the following bash scripts for extracting the Ensembl identifiers (ID) and the score for each transcript predicted by Diana as being targeted for the viral miRNAs, with a prediction score of 70 or higher. 

```bash
sh ../src/extractIDs.sh
```

After running this script, a new directory [`targets/Diana/70`](targets/Diana/70/) will be generated. Inside it, you should obtain the same files that are provided in this repository. Since Diana predicts transcripts and miRDB genes, Ensembl transcripts IDs were mapped to Gene names using the [mapping IDs](src/mappingIDs.Rmd) R notebook.

Once you have the Diana files ready (provided here in the [`targets/Diana/70/`](targets/Diana/70) folder, with the suffix `gene`), the following bash script will help you to combine them with miRDB predictions in order to obtain the set of targets that were predicted by both tools with scores predictions of 70 or higher.

```bash
sh ../src/jointargets.sh
```

After executing the code above, you will obtain eight files called `overlap_<mirna>.tab` with the set of human gene targets predicted for each  `mirna`  of the SARS-CoV-2 miRNA (provided here in the [`targets/`](targets) folder). You can explore the number of targets obtained for each miRNA and generate the Fig2**A** of [[1]](#ref1) by using the [targets exploration notebook](src/targetsExploration.Rmd).


## 6. Analyzing the down-regulation of the predicted targets
### 6.1 Differential expression and functional enrichment analyses

Expression data from RNA-seq experiments involving Calu3 cell-cultures infected with SARS-CoV-2 was downloaded from GEO-NCBI [GSE148729](ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE148nnn/GSE148729/suppl/GSE148729_Calu3_polyA_series1_readcounts.tsv.gz). The [differential expression notebook](src/DifferentialExpression.Rmd) is provided with all the instructions to identify the set of differentially expressed genes. 

After identifying the set of deregulated genes, you will be able to perform and analyze their functional enrichment by means of an overrepresentation analysis with [PANTHER](http://pantherdb.org/). For doing this, the [list of deregulated genes](expression_data/DEgenes_SARS-CoV-2_Full_FC1.txt) should be used as the analyzed list and the [list of expressed genes](expression_data/refGeneList.txt) as the reference list. The overrepresentation analysis is conducted by Fisher's exact test with the Bonferroni adjustment for p-value correction. The results for PANTHER GO-slim biological process (BP), PANTHER Pathways and Reactome Pathways annotation datasets are provided in the [‘functional_enrichment/DEgenes’](functional_enrichment/DEgenes) folder. The [overrepresentation analysis R notebook](src/ORADEGenes.Rmd) will guide you to obtain the graphical representation of the obtained results for the PANTHER GO-Slim BP set, shown in Fig2**C** of [[1]](#ref1).  

### 6.2 Identifying deregulated targets
Once DE analyses have been conducted, the following bash code is used  to find the list of deregulated targets and the set of those potentially being silenced by the predicted viral miRNAs. 

```bash
sh src/deregulatedTargets.sh
```

The previous bash script will create the ‘targets/targetsDE’ folder. Then, it will compare the miRNA targets and the RNA-seq differential expression (DE) results previously obtained (stored in the [‘expression_data/DE_SARS-CoV-2_Full_FC1.csv’](expression_data/DE_SARS-CoV-2_Full_FC1.csv) file). In the ‘targets/targetsDE’ folder, for each predicted SARS-CoV-2 miRNA, two files (‘overlapDE_overlap_SC2V-mir-<xxx>.tab’ and ‘overlapDOWN_overlap_SC2V-mir-<xxx>.tab’) will be generated, with the names of those genes that were found as deregulated and down-regulated for at least one of the differential expression analyses previously performed. The information contained in these files will be summarized, joining the results of all miRNA targets, in the files [‘targetsDE.txt’](targets/targetsDE/targetsDE.txt) and [‘targetsDEDOWN.txt’]](targets/targetsDE/targetsDEDOWN.txt). Finally, the script will also filter the [‘expression_data/DE_SARS-CoV-2_Full_FC1.csv’](expression_data/DE_SARS-CoV-2_Full_FC1.csv) file for keeping only those deregulated targets, generating the [‘targets/targetsDE/targetsDEgenes_SARS-CoV-2_Full_FC1.csv’](targets/targetsDE/TargetsDEgenes_SARS-CoV-2_Full_FC1.csv) file. 

The number of deregulated targets for each miRNA is explored to generate Fig3**A** of the manuscript by using the [differentially expressed targets exploration notebook](src/DETargetsExploration.Rmd).

## 6.3 Functional characterization of deregulated targets
As for analyzing functional enrichment of deregulated genes, the overrepresentation analysis of deregulated targets is also conducted by using [PANTHER](http://pantherdb.org/). Particularly, in this case, the full set of deregulated genes is used as reference list ([‘expression_data/DEgenes_SARS-CoV-2_Full_FC1.txt’](expression_data/DEgenes_SARS-CoV-2_Full_FC1.txt)) and the list of deregulated targets ([‘targets/targetsDE/targetsDE.txt’](targets/targetsDE/targetsDE.txt)) as analyzed list. The results for PANTHER GO-slim BP annotation dataset are provided in the ‘functionalEnrichment/DEtargets’ folder. The [targets overrepresentation analysis R notebook](src/ORADETargets.Rmd) will help you to obtain the graphical representation of the obtained results, shown in Fig3**D** of [[1]](#ref1).

### 6.4 Analysis of down-regulated targets

The list of the 28 targets that were identified as down-regulated potentially being silenced by the novel SARS-CoV-2 miRNAs can be functionally characterized by using the functional classification tool of [PANTHER](http://pantherdb.org/). The results for this step, provided here in the [‘functional_enrichment/DEDOWNTargets/PANTHERClassification.txt’](functional_enrichment/DEDOWNtargets/PANTHERClassification.txt), were manually processed and combined with the DE results to generate a summary file, [‘functional_enrichment/DEDOWNtargets/FunctionalCharacterization.csv’](functional_enrichment/DEDOWNTargets/FunctionalCharacterization.csv). 

The exploration of expression changes and functional characterization of down-regulated targets is carried out using the [down-regulated targets exploration notebook](src/DOWNDETargetsExploration.Rmd). Following it, you will be able to obtain Fig3*B* and Fig3*C* of [[1]](#ref1).  


## 7. Multiple sequence alignment

The candidate sequence found in the previous point were aligned using the [BLAST web platform](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) against the following coronaviruses:

- HKU1 (NCBI Reference Sequence DQ415897.1),
- OC43 (NCBI Reference Sequence KY983585.1),
- Bat-CoV RF1 (NCBI Reference Sequence DQ412042.1), 
- Bat-CoV RATG13 (NCBI Reference Sequence MN996532.1), 
- Malayan Pangolin (NCBI Reference Sequence MT072864.1), 
- SARS-CoV (NCBI Reference Sequence AY278741.1), 
- MERS-CoV (NCBI Reference Sequence MN120514.1), 
- Mice coronavirus (NCBI Reference Sequence GB KF294357.1),
- Ferret coronavirus (NCBI Reference Sequence 1264898). 

Select the mature miRNA extracted from the candidate sequence in the previous section. Then, perform the alignments with [MUSCLE](https://www.ebi.ac.uk/Tools/msa/muscle/) with ClustalW type output format and default settings.

## 8. Figures generation

The [Figures notebook](src/link_Figs.R) is provided with all the instructions to generate Figures 2, 3 and 4 of [[1]](#ref1). 

![Fig3](fig3ab.png)


## References

<a name=ref1></a>[1] G.A. Merino, J. Raad, L.A. Bugnon, C. Yones, L. Kamenetzky, J. Claus, F. Ariel, D.H. Milone, G. Stegmayer, "Novel SARS-CoV-2 encoded small RNAs in the passage to humans," Bioinformatics, 2020, [DOI](https://doi.org/10.1093/bioinformatics/btaa1002).

<a name=ref2></a>[2] C. Yones, N. Macchiaroli, L. Kamenetzky, G. Stegmayer, D.H. Milone, "HextractoR: an R package for automatic extraction of hairpins from genome-wide data", 2020, [DOI](https://doi.org/10.1101/2020.10.09.333898).

<a name=ref3></a>[3] L.A. Bugnon, C. Yones, J. Raad, D.H. Milone, G. Stegmayer,  "Genome-wide hairpins datasets of animals and plants for novel miRNA prediction," Data in Brief, 2019, [DOI](https://doi.org/10.1016/j.dib.2019.104209).

<a name=ref4></a>[4] C. Yones, G. Stegmayer, L. Kamenetzky, D.H. Milone, "miRNAfe: a comprehensive tool for feature extraction in microRNA prediction," BioSystems, 2015, [DOI](http://dx.doi.org/10.1016/j.biosystems.2015.10.003).

<a name=ref5></a>[5] L.A. Bugnon, C. Yones, D.H. Milone, G. Stegmayer, "Genome-wide discovery of pre-miRNAs: comparison of recent approaches based on machine learning," Briefings in Bioinformatics, 2020, [DOI](https://doi.org/10.1093/bib/bbaa184).

<a name=ref6></a>[6] L.A. Bugnon, C. Yones, D.H. Milone, G. Stegmayer, "Deep neural architectures for highly imbalanced data in bioinformatics," IEEE Transactions on Neural Networks and Learning Systems, 2020, [DOI](https://doi.org/10.1109/TNNLS.2019.2914471).

<a name=ref7></a>[7] C. Yones,  J. Raad,  L.A. Bugnon, D.H. Milone, G. Stegmayer, "High precision in microRNA prediction: a novel genome-wide approach based on convolutional deep residual networks," 2020, [DOI](https://doi.org/10.1101/2020.10.23.352179).

