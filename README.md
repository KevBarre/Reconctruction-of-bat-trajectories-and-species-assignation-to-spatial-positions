# Reconctruction of bat trajectories and species assignation to spatial positions
This repository proposes a reconstruction of bat trajectories from spatial positions previously collected applying trigonometry equations to time difference of echolocation call arrivals between differents microphones (01_trajectoryReconstruction.R). After that, a second step assigns automatically identified species (from raw sound files used to compute spatial coordinates) to each position (from any automated identification software) (02_autoIdentificationMatch.R). Here we do not present a method to compute spatial coordinates, but only how to group positions which come from a same trajectory, and link a species assignation from a file containing automatically identified bat passes.

# Scripts used
- 01_trajectoryReconstruction.R
- 02_autoIdentificationMatch.R

# Corresponding data
- 01_dataSample.csv # example to run the script "01_trajectoryReconstruction.R" 
- TriTrajAGarderFull.csv # example of an output of the script "01_trajectoryReconstruction.R" used in the script                          "02_autoIdentificationMatch.R" to assign automatically identified species from a software
- AutoID.csv # automtically identified bat passes (recording were splited in file of 5 seconds duration max to improve the accuracy of the link between a given spatial position and the sound file in which it is contained.
- dataFinal_ProbForest.csv # a dataset obtained from these steps thanks to which we investigated the relationship between the probability of bats being inside the forest (versus open area) in relation with the distance to the light for different spectra (white, red and unlit).
