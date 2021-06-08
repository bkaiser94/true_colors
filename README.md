# true_colors
Ever wanted to know what color that star (or galaxy) would be if you could look through that telescope that's so big it no longer has an eyepiece? Well now using a flux-calibrated spectrum covering (3750-7200 angstroms) and true_colors you can see a representation of the (approximate) color it would be on your computer screen.

Created by Benjamin C. Kaiser (UNC-Chapel Hill) 2021-06-03

Please cite me if you use this for images with something like: "made using true_colors from 
Benjamin C. Kaiser"

**** This is not a "science" code. I make no guarantees that this actually is the real color you 
would see if looking at a given star or blackbody. This code is a good faith attempt to do so,
and I am reasonably confident that the end result is roughly the color a person with normal
color vision would see. Again, it could be wrong; I did it in like 5 days.

***** I was inspired to make this code after seeing this PheT https://phet.colorado.edu/sims/html/blackbody-spectrum/latest/blackbody-spectrum_en.html

source code: https://github.com/phetsims/blackbody-spectrum/blob/4d8b70c78ee85210d7b79661c65c5ba8732c8273/js/blackbody-spectrum/model/BlackbodyBodyModel.js

and the helpful support staff directed me to the source code, but I did *not* use the 
same method they used. I did use their constants and formula for calculating f_lambda for a 
given black body, but that's it (as far as I can remember). It looked like they went directly from the blackbody spectrum to an approximate RGB (which ignores large portions of the visible range and doesn't really calibrate for human eye sensitivity).

***********

My code does the following:

Take an optical spectrum and generate an RGB circle that shows the "true" color of the
spectrum. This could probably be generalized to spectra that do not cover the full visible 
range such that the code infers the slope at those other parts, but initially it's just going to to 
work for my Goodman data of J1644, and it will probably actually use the stitched spectrum 
that I used to make the figure 1 of the Science paper.

The spectrum will be rescaled to the CIE 1931 XYZ system
(https://en.wikipedia.org/wiki/CIE_1931_color_space ;wikipedia link because it's plain-ish 
language, but I didn't use the math here), using linear interpolation between the grid points of the same CIE 1931 XYZ curves as in Wyman et al. 2013
(http://jcgt.org/published/0002/02/01/). CIE 1931 XYZ curves were retrieved from https://www.rit.edu/cos/colorscience/rc_useful_data.php


I'm then using skimage.color.xyz2rgb() to convert these XYZ colors to the sRGB color space, 
which I'm pretty sure is the endzone for the color transformation.

Then I just have to make a circle with that color. I'm 80% sure this last step will inexplicably be 
the most difficult for me based on my extreme difficulties with simple tasks *shrug*.

There are parts of this code that I left in deliberately to show how to edit it for contrived 
spectral inputs to make interesting colors. Basically, I would just make flat spectra at 0 with 
sections with flux=1 to see what colors I could make. 


***** These are not hard and fast requirements for the code to run, but they are the versions I used:

Python 3.6.5
astropy 3.2.3
skimage 0.14.2
numpy 1.19.3