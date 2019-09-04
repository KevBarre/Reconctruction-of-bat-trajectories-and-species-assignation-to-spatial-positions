# Reconctruction of bat trajectories and species assignation to spatial positions
This repository proposes a reconstruction of bat trajectories from spatial positions previously collected applying trigonometry equations to time difference of echolocation call arrivals (TDOA) between differents microphones (01_trajectoryReconstruction.R). After that, a second step assigns automatically identified species (from raw sound files used to compute spatial coordinates) to each position (from any automated identification software, here Tadarida software) (02_autoIdentificationMatch.R). Here we do not present a method to compute spatial coordinates (see Ing et al. 2016 for that), but only how to group positions which come from a same trajectory, and link a species assignation from a file containing automatically identified bat passes.

Reference: Ing RK, Colombo R, Gembu G-C, Bas Y, Julien J-F, Gager Y, Hassanin A. 2016 Echolocation Calls and Flight Behaviour of the Elusive Pied Butterfly Bat (Glauconycteris superba ), and New Data on Its Morphology and Ecology. Acta Chiropterologica 18, 477–488. (doi:10.3161/15081109ACC2016.18.2.014)

# Scripts used
- 01_trajectoryReconstruction.R
- 02_autoIdentificationMatch.R

# Corresponding data
- 01_dataSample.csv # example to run the script "01_trajectoryReconstruction.R" 
- TriTrajAGarderFull.csv # example of an output of the script "01_trajectoryReconstruction.R" used in the script                          "02_autoIdentificationMatch.R" to assign automatically identified species from a software
- AutoID.csv # automtically identified bat passes (recording were splited in file of 5 seconds duration max to improve the accuracy of the link between a given spatial position and the sound file in which it is contained.


