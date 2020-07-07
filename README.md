# Automated consensus-based annotation of cell clusters (ConbAC)

This tool is based on the scRNAseq_Benchmark by tabdelaal 
-> https://github.com/tabdelaal/scRNAseq_Benchmark/tree/master/Snakemake
A Benchmarking classification tools for scRNA-seq data that can help you to select and evaluate the tools for this tool 
-> A combinates tool is planned 

### Sparse Matrix support
This snakemake has been optimzed for sparse matrices as loading and storing of 10x data is quite time and storage extensive. Feel free to update this tool in a way it supports sparse and dense matrices depending on the data! For csv support see the original tool by tabdelaal :) 

## What is it for? 
The Annotation of single cell Transcriptomics defines all further downstream analysis as it is the first step that gives a sense and a meaning to the data. 
Still single cell transcritomics data is very difficult to annotate: Manual marker based annotation requieres good data, expertise and time; automated tools can provide a faster annotation but overfit the data in most of the time and strongly vary between tools.
Especially for malignanbt cells like cancer no tool and neither manual annotation can provide a satisfying annotation. 
CAC should overcome this problem by efficiently combine multiple automated annotation tools and compare the results with prior calculated clusters. 

This tool takes these as an input and returns a 3D confusion table for Tool, cluster and label. From there on I recommend to manualy check the proportions of each cluster for different cell types and annotation tools and choose the annotation suited the best for your data.

I hope this output can help you to provide a better and more suited annotation for your malignant and benign datasets, I recommend to select a deep and well-annotated reference dataset for best results. A subsetting of the anotated reference as well as test sample should improve calculation time without significantly compromising the results.   

## Input
Reference and Sample data in 10xformat as well as a Reference label file (CSV) and a Cluster file for your sample. 
Currently the Reference label file is limited to one column, different reference labels will soon be supported but will require complete new calculations (beside count aggreagtion) for each label. 
Multiple cluster annotation to test different clusters is not supported yet either, but will be soon. This will enable to rerun the ClusterSummary with multiple cluster tables at the same time!

## How to use
[snakemake](https://snakemake.readthedocs.io/en/stable/index.html) and
[singularity](https://www.sylabs.io/docs/) need to be available on your 
system. You will need to run this on a linux system, as singularity
only supports linux.

Singularity is not provided for all functions yet, only those from the original tool. 

From the root of this repository:
```
snakemake \
  --configfile <configfile> \
  --use-singularity
  --cores <N>
```

If your data or output directory is not located under the root of this
repository, be sure to tell snakemake to mount the appropriate directories
in singularity:
```
snakemake \
  --configfile <configfile> \
  --use-singularity \
  --singularity-args '--bind <location of inputs>:<location of inputs> --bind <output directory>:<output directory>'
  --cores <N>
```

#### The config file
```YML
output_dir: <path to outputs directory>
refdatafile: <path to refdatadir in 10x format (matrix.mtx, genes.tsv, barcodes.tsv)>
reflabelfile: <csv with refernece labels per cell>
testdatafile: <path to testdatadir in 10x format (matrix.mtx, genes.tsv, barcodes.tsv)>E
testclusterfile: <csv with cluster per cell for test sample>

number_of_features: <number of features to be used as input for the classification methods, 0 means all, defaults to 0>
tools_to_run:
  - <tool 1>
  - <tool 2>
  - <...>
  
Not supported in this version: 
genes: <path to gene name list, only needed for garnett_CV and Garnett_Pretrained>
column: <The index of the column in the labels file which ought to be used, defaults to 1>


##### Tool specific inputs
Some tools require specific inputs. These tools have not been implemented yet as they lack sparse matrix support!
        one of these tools:
        - Garnett_CV
          ```YML
          Garnett_CV:
            markers: <path to Gernett marker gene file>
          ```
        - Garnett_Pretrained
          ```YML
          Garnett_Pretrained:
            classifier: <path to Gernett classifier>
          ```
        
        <!-- TODO explain these input files -->

## currentyl included  and tested tools/methods
  - scVI
  - scmapcell 
  - Seurat_CCA
  - Seurat_PCA
  - SVM
  - SVM_rejection
  - LDA
  - LDA_rejection
  - RF
  - singleCellNet
  - CHETAH
  - scmapcluster
  - SingleR
  - scID


## Adding new tools
In order to add a tool to this benchmarking workflow, a rule for this tool
needs to be added to the `Snakefile`. This rule should produce as output:
- a table of predicted label (`<output directory/<tool>/<tool>_pred.csv`).
- a table of true labels (`<output directory/<tool>/<tool>_true.csv`).
- a tables of testing, prediction and/or total time:
  - `<output directory>/<tool>/<tool>_test_time.csv`
  - `<output directory>/<tool>/<tool>_training_time.csv`
  - `<output directory>/<tool>/<tool>_total_time.csv`

The input to this rule should be:
- a count table (specified as the `datafile` in the config).
- a true labels file (specified as the `labfile` in the config).

You will want to write a wrapper script for the tool you want to
add to facilitate this. The `"{output_dir}/CV_folds.RData"` input may be
used to provide your wrapper script with folds for cross_validation.
It is recommended to make a docker image containing all dependencies for both
the tool and any wrappers for the tool.  
This wrapper script should also make a selection of the features to be used.
This selection should be based on ranking which can be accessed by providing
`feature ranking` as input to the wrapper script. The number of features to be
used should be configurable and settable through the 'number_of_features' field
in the config.

The following can be used as a template for new rules. Replace everything
surrounded by (and including the) `<>` with appropriate values.

    
rule SingleR:
  input:
    datafile = rules.generate_datafile.output,
    labfile = "{output_dir}/labfile.csv",
    folds = "{output_dir}/CV_folds.RData",
    ranking = feature_ranking
  output:
    pred = "{output_dir}/SingleR/SingleR_pred.csv",
    true = "{output_dir}/SingleR/SingleR_true.csv",
    total_time = "{output_dir}/SingleR/SingleR_total_time.csv"
  log: "{output_dir}/SingleR/SingleR.log"
  params:
    n_features = config.get("number_of_features", 0)
  singularity: "docker://scrnaseqbenchmark/singler:{}".format(dockerTag)
  shell:
    "Rscript Scripts/run_SingleR.R "
    "{input.datafile} "
    "{input.labfile} "
    "{input.folds} "
    "{wildcards.output_dir}/SingleR "
    "{input.ranking} "
    "{params.n_features} "
    "&> {log}"

```
