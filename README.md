# R-M-and-plasmids

Supplementary code for analysis of R-M systems and plasmids. 

For more information, see associated preprint:

Restriction-modification systems have shaped the evolution and distribution of plasmids across bacteria (2022)   
Liam P. Shaw, Eduardo P. C. Rocha, R. Craig MacLean. **biorxiv**  
doi: [10.1101/2022.12.15.520556v1](https://www.biorxiv.org/content/10.1101/2022.12.15.520556v1) 

The `data-generation-scripts` contains python scripts that were used in a cluster environment to generate the intermediate data needed for the main analysis e.g. exceptionality scores for all k-mers across all genomes at multiple subsampling sizes. The resulting data files are large and stored externally at figshare, and should be downloaded into a `data` directory in this directory to run notebooks: [10.6084/m9.figshare.21923121.v1](https://doi.org/10.6084/m9.figshare.21923121.v1)

The `analysis-scripts` folder contains R scripts to generate output figures and tables in the manuscript. These can be run with `./run_all.sh` which should take about 15 minutes on a standard laptop.

By default, all output files (tables and figures) will be date-stamped with the date they are run on e.g. `output-figures/{file_name}-{date}.pdf`.

## Data generation

The data generated that is available via figshare was generated on a computing cluster. The scripts in `data-generation-scripts` were not designed to be used by others, and the data they generate is provided on figshare. However, it should be possible to inspect them to understand the basic approach that we use, and to use and adapt sections of them for your own analysis. 

First, we download all of the genomes for a species from NCBI with `ncbi-genus-download`. Then, we run our analyses. 

We use `makePipelineScript.py` to make the version of the pipeline script for a given species (changing directory names). An example of this is provided for `example_pipeline_Escherichia_coli.py` which can be inspected for the key steps:
* PanACoTA filtering, annotation, and construction of persistent/core genome, as well as fastas for each genome containing only the core genes 
* Running R'MES on the fasta files to get exceptionality scores for every k-mer (k=4,5,6) at multiple values of subsampling the fastas to control for any possible effect of size. (N.B. the word `score` is used for 'exceptionality score'). See `runRMES.py` for a wrapper script to run R'MES.
* Merging those R'MES outputs
* Running rmsFinder on the proteomes for every genome to detect putative R-M systems
* Merging the rmsFinder output

This generates both exceptionality scores for every k-mer and putative R-M systems for a species. This analysis is conducted at multiple subsampling lengths to control for size.  

We then need to combine these together to understand how avoidance of k-mers correlates with R-M system presence across bacterial taxonomy. `calculateAvoidanceKmerGroupsAcrossDatasetv2.py` was used to combine the datasets together and create a 'dictionary' of k-mer targets for each value of k. This dictionary is then used to calculate the median avoidance of different subsets of k-mers for different taxonomic levels and subsamplings. 

## Analysis scripts

We conduct all downstream analysis in R. Scripts in `analysis-scripts` will produce all figures (including supplementary) and tables in the manuscript. To find the code related to a particular figure, consult the list below. Note that some figures were processed in Inkscape for formatting before inclusion in the final paper. 

```
01_avoidance_of_palindromes.R
Figure1
FigureS2
FigureS3

02_heatmaps_targets.R
FigureS1
FigureS4
FigureS5
FigureS6

03_RM_targets.R
Figure2 
FigureS7
FigureS8

04_plasmid_density.R
Figure3
FigureS9

05_plasmid_size_RM.R
FigureS10

06_plasmid_host_range.R
Figure4
FigureS11

07_PTU_modelling.R
Figure5
Figure6
FigureS12
FigureS13
FigureS14
```
