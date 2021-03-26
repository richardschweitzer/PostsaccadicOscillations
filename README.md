# Postsaccadic Oscillations

This repository contains R code developed in the context of the book chapter "Definition, modeling and detection of saccades in the face of post-saccadic oscillations", by Richard Schweitzer and Martin Rolfs. You can find the pre-print here: https://doi.org/10.1101/2021.03.24.436800

At this point it contains:
- Implementations of models to fit saccades with post-saccadic oscillations, as proposed by Bouzat et al. (2018) and Del Punta et al. (2019). 
- The "saccade simulator".
- Direction- and velocity-based PSO detection add-on for the Engbert-Kliegl algorithm for microsaccade detection.
- Saccade and PSO detection using linear classifiers trained on simulated saccade data.

To see how the functions work, see the Rmarkdown file [ModelingDetectionOfPSOs.md](../main/ModelingDetectionOfPSOs.md). 

Please cite as:
R. Schweitzer and M. Rolfs. Definition, modeling and detection of saccades in the face of post-saccadic oscillations. *bioRxiv*, 2021.03.24.436800, 2021. doi: 10.1101/2021.03.24.436800
