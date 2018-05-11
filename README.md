# MeTeOR
## Contents

 - [Project Description](#project-description)
 - [Running from Scratch](#running-from-scratch)
 - [Raw Data and Results](#raw-data-and-results)

## Project Description
The scientific literature is vast, and valuable information connecting findings from disparate works is easily missed. Teams of collaborators address this problem up to a point but could still benefit from systematic “big data” approaches that mine the entire literature to generate testable hypotheses on a large scale.

MeTeOR, or the MeSH Term Objective Reasoning network, mines the PubMed literature, revealing knowledge previously hidden in a sea of information. Given one biological entity (a gene, drug, or disease), it can give a ranked list of associations with other biological entities, and it can highlight papers pertaining to any two biological entities.

This MeTeOR network was assembled with python 3 and it was assessed and predicted upon using MATLAB.

A website serving the resulting network can be found [here](http://meteor.lichtargelab.org/).

## Running from Scratch
There is a [shell script file](src/pipeline.sh) that can be run to assemble MeTeOR and to assess the resulting network. This may be relevant if you wish to have the latest PubMed articles or if you wish to modify some aspect of the creation process. For example, you could create a custom weighting process or create a subnetwork based only on a certain part of the literature. 
### Notes before running the pipeline
#### Requirements
 - [python3](https://www.python.org/downloads/release/python-360/)
 - [MATLAB](https://www.mathworks.com/products/matlab.html) (for prediction and network assessment)
 - [Graphviz](https://www.graphviz.org/) (for the network visualization)

[pyupset](https://github.com/ImSoErgodic/py-upset) is required, but a modified version is already provided in this repository.

#### Time, space, and memory

 - PubMed Data: The data takes a very long time to obtain from PubMed using the provided script, allow 2-3 days depending on download speedsHowever, this script downloads in a query specific manner, allowing the user to customize the query. For bulk download, this project can also be modified to run on [the NLM bulk downloads](https://www.nlm.nih.gov/databases/download/pubmed_medline.html). We provide a basic [python script](MEDLINE/downloadMEDLINE.py) to download this data. The raw XML can take upwards of 200 GB of space.
 - All code was run on an Intel® Core™ i7-4820K CPU @ 3.70GHz × 8 with 64 GB RAM. From start to finish, everything should complete within a week.
 - The Non-negative Matrix Factorization (NMF) conducted in the analysis part of the pipeline and run in MATLAB can be very time and memory intensive. If you chose to, you can download [pre-computed results](#raw-data-and-results) to greatly increase the speed of analysis.

#### Obtaining the data
In order to run the pipeline or assemble the network, certain files are needed. These are available as discussed [below](#raw-data-and-results). 
### Running the whole pipeline
```bash
chmod +x ./run.sh
./run.sh
```
### Generating and characterizing only the network
```bash
python main.py
```
Because this downloads a lot of data, you can specify the storage directory for the PubMed XML.
```bash
python main.py --storagedir /path/to/dir
```
### Run only the assessment part of the pipeline
```bash
cd matlabpipeline
matlab -r pipeline
cd ../EGFR
python runEGFR.py
```

## Raw Data and Results
To host the data and results of the project, I used [dat](https://datproject.org/) and [bulk downloads](http://meteor.lichtargelab.org/download).
### Installing Node and dat
To download data using dat, ensure node (version  >= 4) is installed
```bash
node -v
```
If it needs to be installed. go to their [website](https://nodejs.org/en/download/).

To install dat:
```bash
npm install -g dat
```
### Downloading the data with dat
Navigate to the directory you wish to download, either data or results, and use:
```bash
dat clone ./
```
from the directory of interest.
<!--stackedit_data:
eyJoaXN0b3J5IjpbMjA2MTcwMjQwMF19
-->