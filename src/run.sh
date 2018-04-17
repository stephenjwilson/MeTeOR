

# Run the python part of the pipeline, generating and characterizing the network
python main.py

# Run the MATLAB part of the pipeline, testing the network against other networks
cd matlabpipeline
matlab -r pipeline
cd ../EGFR
python runEGFR.py
