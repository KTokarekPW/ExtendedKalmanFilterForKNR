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


## 31.01.2026 COMMIT Wind Added

Added wind as the 4th variable in the state space, simulated wind using simplified model of a sine wave. Overhauled the simultion massively to aid in tuning all the parameters, managed to achieve MAE of 0.3m, which considering the subpar gps module and IMU is in my opinion and impresssive feat.

# Whats the new code about? It is pretty well documented, hence i will not waste my time explaining the details, ill explain my choices here:

# Why the pseude doppler and not discrete integration when calculating v_gps? 

Discrete integration, has it's drawbacks, we integrate a very noised out signal. Meaning that it leads to enourmous errors, hence i used this algorithm. 
If you are curios how much of a difference does it make, the MAE of position goes up to 0.8 from 0.3 when using the discreet integration so it's basically worthless.

# Why the sine wave as a wind approximation, when there are better approximations?

Yes, but they don't change fundamentally what my filter does and it was quick and easy to code. Essentially my filter doesn't know it's a sine wave, it only knows it changes.

# Why don't you add temperature into the equation?

Temperature only changes the readings (substantially) of one sensor the pitot tube, but mine has a buiilt in system to combat the temperature noise(from what the manufacturer wrote atleast). 

We will see about that in real world testing for my engineering thesis in a year. I highly doubt it's as effective as they proclaim.

# Why the sudden use of functions? 

My code is unreadable without them.

# Why the sudden overhaul and addition of scripts to make plots?

So i can troubleshoot better

# What are you planning to add/do next?

I don't know at this point, i will be probably compiling my code into C so it's faster to run. And then try running it on the Rasperry Pi 5 we have onboard our drone, it will probably lag it to death though.

And i will be doing research for the next few weeks and doing my regular tasks for the software team until then i will pause the work on EKF. Once i know more and find more solutions to make the filter more accurate im gonna tweak it occasionally and that would be it. 

# What code does what?

NaPrzEKF is the old code with 3 states, keeping it mostly for writing my thesis sake to show the evolution.

NaPrzEKFWindTest is the new overhauled 4 state version, which will be tweaked in the coming weeks to improve it.



