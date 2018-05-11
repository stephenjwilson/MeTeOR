# MeTeOR
## Contents

 - [Project Description](#project-description)
 - [Webserver](#webserver)
 - [Running from Scratch](#running-from-scratch)
 - [Download Data and Results](#download)

## Project Description
The scientific literature is vast, and valuable information connecting findings from disparate works is easily missed. Teams of collaborators address this problem up to a point but could still benefit from systematic “big data” approaches that mine the entire literature to generate testable hypotheses on a large scale.

MeTeOR, or the MeSH Term Objective Reasoning network, mines the PubMed literature, revealing knowledge previously hidden in a sea of information. Given one biological entity (a gene, drug, or disease), it can give a ranked list of associations with other biological entities, and it can highlight papers pertaining to any two biological entities.

This MeTeOR network was assembled with python 3 and it was assessed and predicted upon using MATLAB.

## Webserver
A website serving the resulting network can be found [here](http://meteor.lichtargelab.org/).

## Running from Scratch
There is a [shell script file](run.sh) that can be run to assemble MeTeOR and to assess the resulting network. This may be relevant if you wish to have the latest PubMed articles or if you wish to modify some aspect of the creation process. For example, you could create a custom weighting process or create a subnetwork based only on a certain part of the literature. 

Alternatively, you can download the results and use those for your project. This can be done as described [below](#raw-data-and-results):

### Notes before running the pipeline
#### Requirements
 - [python3](https://www.python.org/downloads/release/python-360/)
 - [MATLAB](https://www.mathworks.com/products/matlab.html) (for prediction and network assessment)
 - [Graphviz](https://www.graphviz.org/) (for the network visualization)
 - python3-dev([Ex.](https://packages.ubuntu.com/search?keywords=python3-dev))
 - [npm](#installing-dat)
 - [dat](#installing-dat)

[pyupset](https://github.com/ImSoErgodic/py-upset) is required, but a modified version is already provided in this repository.

#### Time, space, and memory

 - PubMed Data:  This project runs on [the NLM bulk downloads](https://www.nlm.nih.gov/databases/download/pubmed_medline.html). The raw XML can take upwards of 280 GB of space. The pipeline can also be run off specific queries, as was done for the publication version of the code; however, this method takes a very long time to obtain data from PubMed, so please allow 2-3 days depending on download speeds. 
 - All code was run on an Intel® Core™ i7-4820K CPU @ 3.70GHz × 8 with 64 GB RAM. From start to finish, everything should complete a day for bulk downloads or within a week otherwise.
 - The Non-negative Matrix Factorization (NMF) conducted in the analysis part of the pipeline and run in MATLAB can be very time and memory intensive. If you chose to, you can download [pre-computed results](#raw-data-and-results) to greatly increase the speed of analysis.

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

## Download
The MeTeOR network can be download via the [results bulk download](http://static.lichtargelab.org/MeTeOR/results.tar.gz) or as [flat files](http://static.lichtargelab.org/MeTeOR/MeTeORFlatFiles.tar.gz). All downloads are available at [here](http://meteor.lichtargelab.org/download).

For smooth integration into the git project, I also used [dat](https://datproject.org/). The run script automatically downloads the dat data for ease of use. See below for use of dat.
### Installing dat
To download data using dat, ensure node (version  >= 4) is installed:
```bash
node -v
```
If it needs to be installed. go to their [website](https://nodejs.org/en/download/). or:
```bash
curl -sL https://deb.nodesource.com/setup_10.x | sudo -E bash -
sudo apt-get install -y nodejs
```

To install dat:
```bash
sudo npm install -g dat
```
### Downloading with dat
Navigate to the directory you wish to download, either data or results, and use:
```bash
dat clone ./
```
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTEyMzE3MjEzNjIsLTE4NjcwMzYyNjEsLT
IwOTE1NDM3NzIsLTE5MjU5MDA3NjcsLTE0NjY2MDQ5NzBdfQ==

-->