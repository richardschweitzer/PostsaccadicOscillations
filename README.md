# Postsaccadic Oscillations

This repository contains R code developed in the context of the book chapter (pre-print https://doi.org/10.1101/2021.03.24.436800):

Schweitzer, R., Rolfs, M. (2022). Definition, Modeling, and Detection of Saccades in the Face of Post-saccadic Oscillations. In: Stuart, S. (Ed) Eye Tracking. Neuromethods, vol 183. Humana, New York, NY. https://doi.org/10.1007/978-1-0716-2391-6_5. 

At this point the repository contains:
- Implementations of models to fit saccades with post-saccadic oscillations, as proposed by Bouzat et al. (2018) and Del Punta et al. (2019). To play around with the model parameters, you can have a look at this shiny visualization: https://richardschweitzer.shinyapps.io/pso_fitting_example/
- The "saccade simulator".
- Direction- and velocity-based PSO detection add-on for the Engbert-Kliegl algorithm for microsaccade detection.
- Saccade and PSO detection using linear classifiers trained on simulated saccade data.

To see how the functions work, see the Rmarkdown file [ModelingDetectionOfPSOs.md](../main/ModelingDetectionOfPSOs.md). 

If you use or adapt any of the methods, please cite the above book chapter - thank you!
