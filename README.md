# ImagingAnalysis

Requirements:
* numpy
* matplotlib
* cv2
* scipy
* parmap
* image_registration



The scripts function as modules that one can import and utilize as they see fit.

The code relies heavily on the scientific python libraries scipy and numpy to do the averaging, filtering, 
matrices computations and seed pixel correlations. We use cv2 to import mp4 frames as a numpy matrices.
Then parmaps is used to speed up the calculations by allowing us to use multiple cpu cores to do the computations.
Finally image_registration is a great library written by an astronomer (keflavich). It works really well! and it was
used to align the movies between trials.
