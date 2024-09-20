# MLSE-Basketball-FT-Data-Challenge
A basketball free throw physics based trajectory model, that uses optimization to calculate trajectory and will eventually be able to predict basektabll freethrow (FT) success using ML.

# SPL Open Data
The SPL Open Data repository acts as a collection of biomechanics datasets collected by Maple Leaf Sports & Entertainment's (MLSE) Sport Performance Lab (SPL) in Toronto, Ontario, Canada. Through the open-sourcing of this data, SPL's goal is to provide raw markerless motion capture data typically used by sports biomechanists to the general public in an effort to improve data equity and analytical biomechanical skills in the community.

The original dataset contains 125 shots and can be found at (https://github.com/mlsedigital/SPL-Open-Data). 

# Motivation and Description
The data provided by the MLSE SPL has a limited capture volume, which emphasizes the participant's biomechanics during repeated FT attempts. This means that the ball is in the capture volume for a very limited amount of time, with only its landing position in the plane of the hoop being known once the ball leaves the camera's focus. As such, there is no information about the ball's path of flight.

This basketball trajectory model is designed to simulate and analyze the flight path of a basketball shot based on initial conditions such as position and velocity, where these values are extracted or calculated from the last known position of the ball in the capture volume, respectively. The model incorporates key physical factors including gravity, air drag, and the Magnus effect, which are the forces exerted on the ball during flight. The primary goal of the model is to accurately predict the trajectory of the basketball to ensure it meets specified final conditions, such as the final target position and velocity at the end of the shot. To address challenges such as unrealistic apex heights for low initial velocities and ensure the trajectory aligns with desired final conditions, the model includes constraints and penalties to maintain physical realism. Optimization techniques are employed to adjust initial velocities and other parameters, ensuring that the basketball's flight path is realistic and meets the end conditions. Additionally, a machine learning component is currently being developed to predict shot success based on initial conditions, using data from previous trajectories to train the model. This combined approach helps in fine-tuning the model and improving its predictive accuracy for real-world applications. The basketBallTrajectoryCalculation.m script was used to find the trajectories, which uses the MLSE provided data (linked earlier) as input. The trajectories are then saved into a desired folder, along with other shot information, and a shot chart is produced by running the shotChart.m script.

# Data Set Filtering

About 40% of the data provided was discarded for this trajectory model since they produced unrealistic flight paths. When inaccurate trajectories occured, they were often a result of errors in the initial conditions, particularly in the inital velocity velocity. For the faulty trials, the initial calculations for velocity, which were found by dividing the known discplacements by the time steps at the final two known ball positions, produced insufficient velocities, which were much smaller than those reported in literature. These inaccuracies can cause the optimizer to compensate by predicting excessively large velocities, leading to unrealistically high apexes or trajectories that intersect the ground. This highlights the need for improved ball tracking data, to enhance the accuracy of initial velocity calculations. Utilizing high-speed cameras and more sophisticated physics models can provide better estimates of both linear and angular velocities, enabling more precise calculation of spin and drag forces within the trajectory model. 

All trajectories that were deemed unrealistic (with flight apexes greater than 20 ft, not reaching the known final position, or hitting the ground) were discarded and only 71 trajectories were deemed valid. An exmaple of these trajectories is shown below, which were first extracted from the basketBallTrajectoryCalculation.m script which calcualtes the trajectories which are then used as inputs for the shotChart.m script.

# Shot Chart (shotChart.m output) 
 We can see that increasing the shot arc really does lead to greater shot success (even the optimizer agrees). Conclusions like these can be taken directly to players and coaches to help them increase their performance. 
![shotChart_Subject01](https://github.com/user-attachments/assets/0c291b37-64e7-4592-a18d-04d84a994793)

# Individual Shot Analysis (basketBallTrajectoryCalculation.m output for 1 shot)

The newest feature allows us to visualize the various conditions that the ball expereinces in a useful 3D tool! We can see the shot launch and entry angles, in addition to the exit velocity (including the axial components), shot result and total shot time. This tool can be used by teams, atheletes, coaches and trainers to evaluate their freethrow shooting performance and allow the player to re-live their shots after the games. It can ultimately be used to help players gain a better understanding of why their shots miss (or even go in) and how they can imprave. For example, if a player is consistently short on tehir freethrows, that would be reflected by a smaller magnitude of inital exit velocity in the x-direction (towards the basket), indicating a need for them to add more power to their shot, without compromising the arc of the shot (since we can see from the shot chart that increasing the shot arc really does lead to greater shot success).

https://github.com/user-attachments/assets/4adf36ad-10d0-443d-a01e-dc596bcdc6fa



