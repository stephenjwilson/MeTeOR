# MeTeOR
## Contents

 - [Project Description](#project-description)
 - [Running from Scratch](#running-from-scratch)
 - [Raw Data and Results](#raw-data-and-results)

## Project Description
MeTeOR, or the MeSH Term Objective Reasoning network, mines the PubMed literature, revealing knowledge previously hidden in a sea of information.

The scientific literature is vast, and valuable information connecting findings from disparate works is easily missed. Teams of collaborators address this problem up to a point but could still benefit from systematic “big data” approaches that mine the entire literature to generate testable hypotheses on a large scale.

MeTeOR mines the PubMed literature, revealing knowledge previously hidden in a sea of information. Given one biological entity (a gene, drug, or disease), it can give a ranked list of associations with other biological entities, and it can highlight papers pertaining to any two biological entities.

This project is built off of python 3 for the creation of the MeTeOR network and MATLAB for the assessment of the network as well as the unsupervised link prediction.

A website serving the resulting network can be found [here](http://meteor.lichtargelab.org/).

## Running from Scratch

## Raw Data and Results
To host the data and results of the project, I used [dat](https://datproject.org/) and [bulk downloads](http://meteor.lichtargelab.org/download).
### Installing Node and dat
To download data using dat, ensure node is installed and the version is higher than 4
```
node -v
```
If it needs to be installed got to their [website](https://nodejs.org/en/download/).
To install dat:
```bash
npm install -g dat
```
### Downloading the data
Navigate to the directory you wish to download, either data or results, and use:
```bash
dat clone ./
```
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTY0NDE5NzUzMiwtMTc0MzY1ODQyMiwxMj
MxODgzODY3LC0xOTM2NDAwMjAxLC0xMTM1NjA0OTNdfQ==
-->