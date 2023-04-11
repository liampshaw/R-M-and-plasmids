# R-M-and-plasmids

Supplementary code for analysis of R-M systems and plasmids. 

For more information, see associated preprint:

Restriction-modification systems have shaped the evolution and distribution of plasmids across bacteria (2022)   
Liam P. Shaw, Eduardo P. C. Rocha, R. Craig MacLean. **biorxiv**  
doi: [10.1101/2022.12.15.520556v1](https://www.biorxiv.org/content/10.1101/2022.12.15.520556v1) 

Data files are large and stored externally at figshare, and should be downloaded into a `data` directory in this directory to run notebooks: [10.6084/m9.figshare.21923121.v1](https://doi.org/10.6084/m9.figshare.21923121.v1)

By default, all output files (tables and figures) will be date-stamped with the date they are run on e.g. `output-figures/2023-03-29-{file_name}.pdf`.

This can be run with `./run_all.sh` which should take about 15 minutes on a standard laptop.

## Data generation

The data generated that is available via figshare was generated on a computing cluster. The scripts in `scripts` are not designed to be used by others. However, it should be possible to use parts of them, and to understand the basic approach that we use. 

First, we download all of the genomes for a species from NCBI with `ncbi-genus-download`. Then, we run our analyses. 

We use `makePipelineScript.py` to make the version of the pipeline script for a given species (changing directory names). An example of this is provided for `example_pipeline_Escherichia_coli.py` which can be inspected for the key steps:
* PanACoTA filtering, annotation, and construction of persistent/core genome, as well as fastas for each genome containing only the core genes 
* Running R'MES on the fasta files to get exceptionality scores for every k-mer (k=4,5,6) at multiple values of subsampling the fastas to control for any possible effect of size. (N.B. the word `score` is used for 'exceptionality score'). See `runRMES.py` for a wrapper script to run R'MES.
* Merging those R'MES outputs
* Running rmsFinder on the proteomes for every genome to detect putative R-M systems
* Merging the rmsFinder output

This generates both exceptionality scores for every k-mer and putative R-M systems for a species. This analysis is conducted at multiple subsampling lengths to control for size.  

We then need to combine these together to understand how avoidance of k-mers correlates with R-M system presence across bacterial taxonomy. `calculateAvoidanceKmerGroupsAcrossDatasetv2.py` was used to combine the datasets together and create a 'dictionary' of k-mer targets for each value of k. This dictionary is then used to calculate the median avoidance of different subsets of k-mers for different taxonomic levels and subsamplings. 

# Correspondence with manuscript

The correspondence between the manuscript and the notebooks is as follows:
* `01_avoidance_of_palindromes.Rmd` - Avoidance of 6-bp palindromes is stronger in plasmid genes than in core genes. 

| Item 	| File names 	| Notes 	|
|---	|---	|---	|
| Figure 1 	| `figure-1-k6-subsampling-50000.pdf`  and  `variation-pangenome-components.pdf` 	|  	|
| Figure 2 	|  	|  	|
| Figure 3  | `escherichia-coli-plasmids.pdf` (a, b) and   `plasmid-boxplots.pdf` (c, d, e). But now `2023-04-04-plasmid-boxplots-figure-3.pdf`         |  `04_plasmid_density.Rmd`     | 
| Figure 4 	| `ptus-palindromic.pdf` 	| In `05_plasmid_host_range.R` 	|
| Figure 5 	|  	|  	|
| Figure 6 	| | |  
