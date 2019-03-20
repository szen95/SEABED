# SEABED
The repository for the SEABED project. SEABED stands for SEgmentation And Biomarker Enrichment for defining Differential drug response.

### Workflow

This file covers the steps to get from the public data files (Table S2) to the segmented subpopulations and the relevant biomarkers generated (Table S3) by SEABED.

The scripts are named as they are used.

1. Get the data (from GDSC, CCLE, CTRP) via Table S2

2. Segmentation code

3. Feature enrichment to nominate biomarkers

4. Classification of pair-wise drug responses

5. Visualization of pair-wise drug responses
* 'heatmap.py'

6. 2-D visualization of drug response profiles
* 'scatter_plot.py' takes the output data after segmentation and creates a scatter plot showing the differential drug response of the subpopulations.
* It also annotates subpopulations with enriched biomarker(s).

7. Tree diagram visualization of subpopulations
* To generate tree diagrams, call the 'create_tree' method in 'tree_diagram.py'. This creates tree diagrams showing the hierarchical sequence in how the subpopulations are formed during the segmentation process.
* 'create_tree_image.py' and 'modify_tree_image.py' are used to customize the appearance of the tree diagram, including the parent and terminal nodes.
