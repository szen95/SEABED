# SEABED
The repository for the SEABED project. SEABED stands for SEgmentation And Biomarker Enrichment for defining Differential drug response.

![logo](https://github.com/szen95/SEABED/blob/master/website_s1/img/SEABED_logo.png)

### Workflow

This file covers the steps to get from the public data files (Table S2) to the segmented subpopulations and the relevant biomarkers generated (Table S3) by SEABED.

The scripts are named as they are used.

1. Get the data (from GDSC, CCLE, CTRP) via Table S2

2. Segmentation code

3. Feature enrichment to nominate biomarkers

4. Classification of pair-wise drug responses
* `drug_pair_classification.R`

5. Visualization of pair-wise drug responses
* `heatmap.py` takes the matrix in `MAPK_AKT_drug-drug_matrix.csv` (this is found in `example_inputs`) and creates a heatmap that shows the distinct drug response types.
*  The relevant files to generate the supplemental website are `index.html` and in the folder `website_s1`.
* The SEABED supplemental website can be found [here](https://szen95.github.io/SEABED/).  

6. 2-D visualization of drug response profiles
* `scatter_plot.py` takes the output data after segmentation and creates a scatter plot showing the differential drug response of the subpopulations.
* It also annotates subpopulations with significantly enriched biomarker(s).

7. Tree diagram visualization of subpopulations
* To generate tree diagrams, call the `create_tree` method in `tree_diagram.py`. This creates tree diagrams showing the hierarchical sequence in how the subpopulations are formed during the segmentation process.
* `create_tree_image.py` and `modify_tree_image.py` are used to customize the appearance of the tree diagram, including the parent and terminal nodes.
