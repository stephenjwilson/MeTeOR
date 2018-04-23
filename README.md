# MeTeOR
## Contents

 - [Project Description](#project-description)
 - [Running from Scratch](#running-from-scratch)
 - [Raw Data and Results](#raw-data-and-results)

## Project Description
The scientific literature is vast, and valuable information connecting findings from disparate works is easily missed. Teams of collaborators address this problem up to a point but could still benefit from systematic “big data” approaches that mine the entire literature to generate testable hypotheses on a large scale.

MeTeOR, or the MeSH Term Objective Reasoning network, mines the PubMed literature, revealing knowledge previously hidden in a sea of information. Given one biological entity (a gene, drug, or disease), it can give a ranked list of associations with other biological entities, and it can highlight papers pertaining to any two biological entities.

This project is built off of python 3 for the creation of the MeTeOR network and MATLAB for the assessment of the network as well as the unsupervised link prediction.

A website serving the resulting network can be found [here](http://meteor.lichtargelab.org/).

## Running from Scratch
There is a [shell script file](src/pipeline.sh) that can be run to assemble MeTeOR and to assess the resulting network. This may be relevant if you wish to have the latest PubMed articles or if you wish to modify some aspect of the creation process. For example, you could create a custom weighting process or create a subnetwork based only on a certain part of the literature. 
### Notes before running the pipeline
#### Requirements
#### Time and memory
 - Time
	 - PubMed Data: This data takes a very long time (2-3 days) depending on download speeds to obtain from PubMed using the script provided. However, this script downloads in a query specific manner, allowing the user to customize the query. For bulk download, that this project can be modified to run on see [the NLM download [page](https://www.nlm.nih.gov/databases/download/pubmed_medline.html].
	 - 

### Running the whole pipeline
From the src folder:
```bash
./pipeline.sh
```
If necessary:
```bash
chmod +x ./pipeline.sh
```
### Run the python part of the pipeline, generating and characterizing the network
```bash
python main.py
```
### Run the assessment part of the pipeline, testing the network against other networks. Includes MATLAB code.
```bash
cd matlabpipeline
matlab -r pipeline
cd ../EGFR
python runEGFR.py
```

## Raw Data and Results
To host the data and results of the project, I used [dat](https://datproject.org/) and [bulk downloads](http://meteor.lichtargelab.org/download).
### Installing Node and dat
To download data using dat, ensure node is installed and the version is higher than 4
```bash
node -v
```
If it needs to be installed got to their [website](https://nodejs.org/en/download/).
To install dat:
```bash
npm install -g dat
```
### Downloading the data with dat
Navigate to the directory you wish to download, either data or results, and use:
```bash
dat clone ./
```
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTI0MDA4Mjk1MSwxMTIwNTQyNDg0LDU3ND
M4NDkyMywtMTAwNDk0ODI1NywxMDA2OTUxODYwLDExNTEyMDY3
MTIsLTE3NDM2NTg0MjIsMTIzMTg4Mzg2NywtMTkzNjQwMDIwMS
wtMTEzNTYwNDkzXX0=
-->