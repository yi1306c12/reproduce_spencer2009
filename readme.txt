This is the readme for the model associated with the paper:

The Functional Consequences of Cortical Circuit Abnormalities on Gamma
Oscillations in Schizophrenia: Insights from Computational Modeling
Kevin M. Spencer
Front Hum Neurosci. 2009; 3: 33.
doi:  10.3389/neuro.09.033.2009
PMCID: PMC2769552

These are the files needed to run the model and analyze the
output. Below is a brief description of the files which run under IDL.

----------------------------------------------------------------------
Kevin M. Spencer, Ph.D.  Director, Neural Dynamics Laboratory
(http://ndl.hms.harvard.edu) Research Health Scientist, VA Boston
Healthcare System Associate Professor of Psychiatry, Harvard Medical
School
----------------------------------------------------------------------

setup.txt:
The user reads in this batch file to set some basic parameters before
doing anything else.
 
ctx:
Master program. Initializes random number generation, sets up various
variables and constants, then calls "init_net" to initialize the
network structure and "noise_gen" to generate the random noise
input. Certain simulations of SZ changes are implemented here (such as
NMDA input strength). Then runs the model. Integration is Euler, .001
ms time step, which is downsampled to 1 ms. The argument "str" is used
to give the run a unique name. 3 files are created: str_rand.dat (seed
of random number generator), str_spk.dat (spike data), and str_Vm.dat
(membrane voltage data).
 
init_net:
Sets up the weight matrix. Simulations of SZ changes to various
aspects of the weight matrix (patterns, weights) are done here. 2
files are created: str_param.dat (various parameters of the run), and
str_weights.dat (the weight matrix).
 
noise_gen:
Generates random noise input to network, saves in str_noise.dat.
 
netanal:
This takes the output of a run and processes the data, displaying
things like power spectra and firing rates.
 
results_new:
This is used to create the condition X frequency plots in the
paper. It basically does the analyses in "netanal" for each condition
(say 10-90% reduction of FSI output). It puts out files that are
displayed by "werstf".
 
werstf:
Displays the contour plots of condition X frequency. This was adapted
from the program I use to display EEG time/frequency plots.
 
wers_thresh:
This isn't used, but werstf will complain if it's not there.
 
c_sims.txt:
This is an example batch file for running several of the simulations
in the paper, so you can see the syntax.
