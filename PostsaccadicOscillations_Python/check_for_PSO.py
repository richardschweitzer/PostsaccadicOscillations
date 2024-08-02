#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Minimal PSO onset detection according to the algorithm described in Figure 5 of
Schweitzer & Rolfs "Definition, modeling and detection of saccades in the face 
of post-saccadic oscillations" (https://doi.org/10.1101/2021.03.24.436800),  

Created on Thu Mar 30 15:35:06 2023

@author: richard schweitzer
"""

import numpy as np


class check_for_PSO_class:
    
    def __init__(self, 
                 SAMPLING_RATE=1000, FIVE_POINT_SMOOTH=True, 
                 DIRECTION_CUTOFF=90):
        # set these global parameters, we don't need to set them every time
        print("Welcome to PSO onset detection based on direction inversion!")
        print("Please acknowledge the work by citing the original paper:\n", 
              "Schweitzer, R., & Rolfs, M. (2022). Definition, Modeling, and Detection of Saccades in the Face of Post-saccadic Oscillations. In Eye Tracking: Background, Methods, and Applications (pp. 69-95). New York, NY: Springer US.\n")
        self.set_params(SAMPLING_RATE, FIVE_POINT_SMOOTH, DIRECTION_CUTOFF)
        

    def set_params(self, 
                   SAMPLING_RATE=1000, FIVE_POINT_SMOOTH=True, DIRECTION_CUTOFF=90):
        """
        Sets the global parameters for PSO detection
        """
        self.SAMPLING_RATE = SAMPLING_RATE
        self.FIVE_POINT_SMOOTH = FIVE_POINT_SMOOTH
        self.DIRECTION_CUTOFF = DIRECTION_CUTOFF
        print("Set parameters for PSO detection:\n SAMPLING_RATE =", str(self.SAMPLING_RATE), "\n", 
              "DIRECTION_CUTOFF =", str(self.DIRECTION_CUTOFF), "\n", 
              "FIVE_POINT_SMOOTH =", str(self.FIVE_POINT_SMOOTH), "\n")
    
    
    def vecvel(self, x, y):
        """
        Computes a 2D velocity vector. Code is translated from the original R function 
        by Ralf Engbert, Petra Sinn, Konstantin Mergenthaler, and Hans Trukenbrod
        """
        # create a Nx2 matrix to match the existing code
        xy = np.matrix(np.vstack((np.array(x), np.array(y)))) 
        xy = np.matrix.transpose(xy) 
        # determine size and preallocate the velocity vector
        d = np.shape(xy)
        N = d[0] # Python, unlike R, starts at 0 as first vector index - this is fixed below
        v = np.matrix(np.zeros(d))
        # compute velocity
        if self.FIVE_POINT_SMOOTH==True:
            v[2:(N-3),] = self.SAMPLING_RATE/6 * (xy[4:(N-1),] + xy[3:(N-2),] - xy[1:(N-4),] - xy[0:(N-5),])
            v[1,] = self.SAMPLING_RATE/2 * (xy[2,] - xy[0,])
            v[(N-2),] = self.SAMPLING_RATE/2 * (xy[(N-1),] - xy[(N-3),])   
            # this has to be added for python compatibility, as indexing does not include the last element
            v[(N-3),] = self.SAMPLING_RATE/6 * (xy[(N-1),] + xy[(N-2),] - xy[(N-4),] - xy[(N-5),])
        else:
            v[1:(N-2),] = self.SAMPLING_RATE/2 * (xy[2:(N-1),] - xy[0:(N-3),])
            # this has to be added for python compatibility:
            v[(N-2),] = self.SAMPLING_RATE/2 * (xy[(N-1),] - xy[(N-3),])
        return(v)   
    
    
    def vecdir(self, x, y):
        """
        Computes a vector of sample-to-sample directions
        """
        v = self.vecvel(x=x, y=y)
        N = np.shape(v)[0] - 1
        # compute directions
        deg = np.array(np.arctan2(v[0:N, 1], v[0:N, 0]) * 180/np.pi)
        deg[deg < 0] = deg[deg < 0] + 360
        return(deg)
    
    
    def get_PSO_onset(self, x, y, onset_index, offset_index):
        """
        This gets us the onset of the PSO according to the algorithm described 
        in Schweitzer & Rolfs "Definition, modeling and detection of saccades 
        in the face of post-saccadic oscillations", mostly in Figure 5. 
        onset_index and offset_index are the indeces demarcating the saccade.
        """
        eye_x = np.array(x)
        eye_y = np.array(y)
        assert(len(eye_x)==len(eye_y))
        # what is the overall direction of the saccade?
        overall_sac_direction = np.arctan2(eye_y[offset_index]-eye_y[onset_index], 
                                           eye_x[offset_index]-eye_x[onset_index]) * 180/np.pi
        if overall_sac_direction < 0: 
            overall_sac_direction = overall_sac_direction + 360
        # what is the direction criterion?
        opposite_direction_limits = np.array([overall_sac_direction+self.DIRECTION_CUTOFF, # max
                                              overall_sac_direction-self.DIRECTION_CUTOFF]) # min
        opposite_direction_limits[opposite_direction_limits>360] = \
            opposite_direction_limits[opposite_direction_limits>360] - 360
        opposite_direction_limits[opposite_direction_limits<0] = \
            opposite_direction_limits[opposite_direction_limits<0] + 360
        opposite_direction_limits[opposite_direction_limits==0] = \
            opposite_direction_limits[opposite_direction_limits==0] + np.nextafter(0, 1)
        # where should we start looking for the inversion point?
        start_looking = round(np.floor(onset_index+(offset_index-onset_index)/4.0))
        # compute direction vector
        sac_directions = self.vecdir(x = x, y = y)
        sac_directions = np.squeeze(sac_directions)
        # compute velocity vector
        sac_velocities = self.vecvel(x = x, y = y)
        sac_velocities_x = np.squeeze(np.array(sac_velocities[0:,0]))
        sac_velocities_y = np.squeeze(np.array(sac_velocities[0:,1]))
        sac_velocities_abs = np.sqrt(np.square(sac_velocities_x) + 
                                     np.square(sac_velocities_y))
        # determine the minimum velocity within the range that we'll search
        # we'll not consider the last 3 samples of the saccade, as they'll have
        # low velocity anyway, but are unlikely to signal a PSO onset
        absolute_min_vel = np.min(sac_velocities_abs[start_looking:(offset_index-2)])
        # preallocate:
        direction_inversion = np.empty(len(eye_x))
        direction_inversion[:] = np.nan
        minimum_velocity = np.empty(len(eye_x))
        minimum_velocity[:] = np.nan
        # now extract the point, where direction deviates from the overall 
        # saccade direction, i.e., likely the PSO
        for it in range(start_looking, offset_index):
            # check directions
            if opposite_direction_limits[0] > opposite_direction_limits[1]: # max > min, basically any saccade
                if (sac_directions[it] < opposite_direction_limits[0] and
                    sac_directions[it] > opposite_direction_limits[1]): 
                    direction_inversion[it] = 0
                else:
                    direction_inversion[it] = 1
            else: # min > max, the case for a rightward saccade
                if ((sac_directions[it] >= 0.0 and sac_directions[it] < opposite_direction_limits[0]) or
                     (sac_directions[it] <= 360.0 and sac_directions[it] > opposite_direction_limits[1])):
                    direction_inversion[it] = 0
                else:
                    direction_inversion[it] = 1
            # check velocities, is the current velocity the minimum velocity, but not the last above-threshold sample?
            minimum_velocity[it] = sac_velocities_abs[it]==absolute_min_vel
        # DIRECTION: what's the first direction mismatch? (if there's any)
        found_inversion = (np.isnan(direction_inversion)==False) & (direction_inversion==1)
        if np.any(found_inversion):
            first_dir_mismatch = np.min(np.where(found_inversion))
        else:
            first_dir_mismatch = np.nan
        # VELOCITY: where's our minimal velocity (if not the last sample)
        found_min_vel = (np.isnan(minimum_velocity)==False) & (minimum_velocity==1)
        if np.any(found_min_vel):
            first_min_vel = np.min(np.where(found_min_vel))
        else:
            first_min_vel = np.nan
        # prepare output
        output_now = {'PSO_dir_index': first_dir_mismatch,  # PSO onset index according to direction criterion
                      'overall_sac_dir': overall_sac_direction, # saccade direction (offset re onset positions)
                      'PSO_minvel_index': first_min_vel,  # PSO onset index according to minimum velocity criterion (not recommended)
                      'min_sac_vel': absolute_min_vel,  # minimum velocity found for PSO min vel criterion
                      'sac_onset_index': onset_index,  # saccade onset index
                      'sac_offset_index': offset_index, # saccade offset index
                      'sac_directions': sac_directions, # vector of sample-to-sample directions
                      'sac_velocities_xy': sac_velocities, # 2D vector of directional saccade velocities
                      'sac_velocities_abs': sac_velocities_abs # absolute saccade velocities
                      }
        # output
        return(output_now)
    


### DEMO 
if __name__ == '__main__': 
    
    # packages we need for the demo
    import matplotlib.pyplot as plt
    import pandas as pd
    
    # load hand-labeled saccade data, just one saccade.
    # it's the saccade from Fig. 5 in Schweitzer & Rolfs (2022):
    #    Definition, modeling and detection of saccades in the face of post-saccadic oscillations
    # (handlabeled data from Nystroem et al, see: https://osf.io/n36fx/)
    sac_data = pd.read_csv("some_saccade.csv")
    x = sac_data.x
    y = sac_data.y
    samp_rate = 1000.0/np.median(sac_data.samp_freq)
    print("sampling rate:", samp_rate, "Hz")
    # plot the x-y saccade trajectory
    plt.plot(x, y)
    plt.show()
    
    # determine saccade onset and offset based on labels
    sac_onset_index = np.min(np.where(sac_data.what=="SAC"))
    sac_offset_index = np.max(np.where(sac_data.what=="PSO"))
    sac_PSO_index = np.min(np.where(sac_data.what=="PSO"))
    
    # import the detection module and set sampling rate and smoothing
    import check_for_PSO
    check_PSO = check_for_PSO.check_for_PSO_class(SAMPLING_RATE=samp_rate, FIVE_POINT_SMOOTH=True)
    
    # compute 2D velocity to see whether the output is correct
    v = check_PSO.vecvel(x, y)
    # ... and plot the x-y velocity
    plt.plot(v[0:, 0], v[0:, 1])
    plt.show()
        
    # now find out where the PSO onset is ...
    PSO_output = check_PSO.get_PSO_onset(x, y, sac_onset_index, sac_offset_index)
    PSO_output
    print("Found PSO at", PSO_output['PSO_dir_index'], "where handlabeled onset was at", sac_PSO_index, " - see the plots...")
    
    # plot the direction over time
    plt.plot(PSO_output['sac_directions'])
    plt.vlines([sac_onset_index, sac_PSO_index, sac_offset_index], 0, 360, colors="black") # hand-labeled
    plt.vlines([PSO_output['PSO_dir_index'], PSO_output['PSO_minvel_index']], 0, 360, colors="red") # detected
    plt.show()
    
    # plot the absolute velocity over time
    plt.plot(PSO_output['sac_velocities_abs'])
    plt.vlines([sac_onset_index, sac_PSO_index, sac_offset_index], 
               0, max(PSO_output['sac_velocities_abs']), colors="black") # hand-labeled
    plt.vlines([PSO_output['PSO_dir_index'], PSO_output['PSO_minvel_index']], 
               0, max(PSO_output['sac_velocities_abs']), colors="red") # detected
    plt.show()
    
    # plot distance traveled over time
    dist_traveled_over_time = np.sqrt(np.square(x-x[sac_onset_index]) + np.square(y-y[sac_onset_index]))
    plt.plot(dist_traveled_over_time)
    plt.vlines([sac_onset_index, sac_PSO_index, sac_offset_index], 
               0, max(dist_traveled_over_time), colors="black") # hand-labeled
    plt.vlines([PSO_output['PSO_dir_index'], PSO_output['PSO_minvel_index']], 
               0, max(dist_traveled_over_time ), colors="red") # detected
    plt.show()
    print("Done.")
    
    