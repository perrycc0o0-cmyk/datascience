# Rat infralimbic cortex snRNA-seq analysis

Analysis of single-nucleus RNA sequencing data from the rat infralimbic cortex after withdrawal alone (WD), extinction (Ext) or retrieval-extinction (Ret-Ext). Each condition has three pooled samples, with tissue from five rats per pool.

[Infralimbic_snRNA_R](Infralimbic_snRNA_R/) contains the numbered analysis notebooks, from library QC and cell annotation to differential expression, trajectory, enrichment, subtype proportion and neuronal communication analyses.

## Usage

1. See the [analysis guide](Infralimbic_snRNA_R/README.md) for inputs, package versions and notebook order.
2. Open the required `.Rmd` in RStudio. Install the packages listed in its first section, then edit the paths and check the sample metadata and settings near the top. Where present, set `script_path` to the full path of that notebook.
3. In a fresh R session, run the R chunks in order. For the gene-category network, edit the paths in the Python script and run it with Python 3 before opening the plotting notebook.
4. Use the saved RDS, CSV and figure files in the specified output directories. Keep sequencing matrices, external references and generated results outside the code directory.

The original code and documentation are available under the [MIT License](LICENSE).
