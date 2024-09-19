# MLSE-Basketball-FT-Data-Challenge
A basketball free throw physics based trajectory model, which will eventually be able to predict basektabll freethrow (FT) success using ML.

# SPL Open Data
The SPL Open Data repository acts as a collection of biomechanics datasets collected by Maple Leaf Sports & Entertainment's (MLSE) Sport Performance Lab (SPL) in Toronto, Ontario, Canada. Through the open-sourcing of this data, SPL's goal is to provide raw markerless motion capture data typically used by sports biomechanists to the general public in an effort to improve data equity and analytical biomechanical skills in the community.

The original dataset contains 125 shots and can be found at (https://github.com/mlsedigital/SPL-Open-Data). 

# Motivation and Description
The data provided by the MLSE SPL has a limited capture volume, which emphasizes the participant's biomechanics during repeated FT attempts. This means that the ball is in the capture volume for a very limited amount of time, with only its landing position in the plane of the hoop being known once the ball leaves the camera's focus. As such, there is no information about the ball's path of flight.

This basketball trajectory model is designed to simulate and analyze the flight path of a basketball shot based on initial conditions such as position and velocity, where these values are extracted or calculated from the last known position of the ball in the capture volume, respectively. The model incorporates key physical factors including gravity, air drag, and the Magnus effect, which are the forces exerted on the ball during flight. The primary goal of the model is to accurately predict the trajectory of the basketball to ensure it meets specified final conditions, such as the final target position and velocity at the end of the shot. To address challenges such as unrealistic apex heights for low initial velocities and ensure the trajectory aligns with desired final conditions, the model includes constraints and penalties to maintain physical realism. Optimization techniques are employed to adjust initial velocities and other parameters, ensuring that the basketball's flight path is realistic and meets the end conditions. Additionally, a machine learning component is currently being developed to predict shot success based on initial conditions, using data from previous trajectories to train the model. This combined approach helps in fine-tuning the model and improving its predictive accuracy for real-world applications.

