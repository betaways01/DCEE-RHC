# DCEE + RHC Path Planning

## A particle filter implementation for an Unmanned Aerial Vehicle (UAV), attempting to locate a radioactive source in a 3D environment. 

## High-level overview:

The UAV is defined with a specific starting position, velocity, and movement distance.

The environment (i.e., the search area) is defined as a 3D grid. The radioactive source is defined with specific characteristics, including its coordinates, the rate of radioactive release, wind speed, and specific constants related to the dispersion model.

The UAV begins a loop (for 100 iterations, or until the battery limit is reached, or until the uncertainty in the source location is below a threshold).

In each loop, it first generates a 'measurement' of the radioactive level at its current location by using a radioactive dispersion model and then adds some random noise to this 'measurement'. This model represents the actual physical process of how the radioactive material spreads out from the source.

Using the measurement, the UAV then updates its beliefs about the location and characteristics of the source using a particle filter approach. This involves having multiple 'particles' which each represent a hypothesis of the source's state (e.g., location, rate of radioactive release). These particles are updated based on how well they predicted the measurement, with particles that did a better job being given more weight in future predictions.

The UAV then evaluates potential new locations to move to, which are adjacent to its current location. It evaluates these based on the expected reduction in uncertainty about the source's location if a measurement were to be taken at these locations. This is essentially trying to determine the next step that would give the most information.

The UAV chooses the location that results in the greatest reduction in uncertainty and moves there.

This process repeats until the loop ends.
