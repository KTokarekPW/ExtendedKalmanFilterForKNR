# What is the code even about?

My code is a simple simulation of all the hardware for idealized flight of a drone along the ZY plane, along with my KF trying to determine the position, velocity as well as the bias of the IMU as 3 variables in the state space.

# How does it work?

I feed it p_true and all the derivatives, meaning i feed the programme the path of the drone that it's gonna take, then i simulate the sensors readings on the basis of what the manufacters stated in their documentation.

Then i simulate my filter 10000 times to get really good mean and median errors of the filter in order to test it's effectiveness.

# Why did you use a pitot tube for speeds of above 2m/s?

I'm not using a pitot tube, im doing something a bit out of ordinary i plan on using a very sensitive pressure sensor instead of the pitot tube. This allows me for very good readings in speeds lower than 10m/s.

# Does that resemble the real world even slightly?

No, but we have to start somewhere. 

# How will i improve my model to better work outside of simulation?

I will be adding state spaces for temperature gradient that i will simply model as a change of up to 0.1 degree kelvin per meter travelled randomly. To simulate non even temperature of air outside. I tiwll affect air density and the pitot tube.

Additionally i will be adding wind, which will work very similarly but will be a sine wave function over time, that will change the reading of my pressure sensor by the pressure change from the wind.

