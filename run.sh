
# Check for the data
file='./data/c2018.bin'
if [ ! -e "$file" ]; then
    echo "Please Download the raw data need to run the code"
    echo "Attempting to download now..."
    cd data
    dat clone
fi 

cd src
# Run the python part of the pipeline, generating and characterizing the network
python main.py --storagedir /mammothshared/ProjectStructures/MeTeOR/MEDLINE

file='../results/pmid_ui.npz'
if [ ! -e "$file" ]; then
    echo "Python pipeline failed"
    exit 1
fi


# Run the MATLAB part of the pipeline, testing the network against other networks
cd matlabpipeline
matlab -r pipeline

file='../results/validation/mmr300_MeTeOR1956.txt'
if [ ! -e "$file" ]; then
    echo "MATLAB pipeline failed"
    exit 1
fi

cd ../EGFR
python runEGFR.py
