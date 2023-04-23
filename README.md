# CaBMI_PrairieView
Code for running CaBMI with Prairie_View and Matlab.
The first initial version of this code was written by Vivek, Ines and Nuria. Nuria has since adapted and modified a bit here and there.
 
## Prairie view setup
Prairie view is set to run in resonant scanner. Matlab should be able to connect with the API and send commands. For that, there is a small script that needs to be set up in the Prairie software as below:
```
srd 50
lbs 5
```
check prairie link documentation for more information.

On top of that, Prairie needs to record voltage rec samples to synchronize the BMI frames with the prairie frames.

There are some environment files that need to be defined a priori. See some examples on utils.
At least 2 environments need to exist, one for baseline (of about 15 min) and one for BMI (30-60min).
 
## Protocol
The code is run through the "main_protocol" script. They are a step-by-step instructions on how to run the BMI.
Matlab needs to have a clean slate, or at the very least, make sure that an arduino instance doesn't exist.
At the very beginning, there are some variables (animal, day, date, experiment type) to create a path that can be change from here. 
So the "parameter's files" of the BMI can be kept constant and no need to be changed. 
The only thing that will change from day to day will be those variables in the main protocol script.
The expt_str variable is the type of experiment to run (go to expt2bmi_flags.m for definitions).
Parameters can also be changed on the main matlab workspace. For example:
```
tset.flag_delay = 1;
fbset.fb_bool = 0;
```
The first line will introduce a delay of tset.delay_time seconds.
The second line turns off the audio to remove the feedback to the mouse.
In the first run, define the main folders for environment or to store the results of the BMI.
Then the code will connect to the Prairie API to obtain the size of the image, the first image, etc.
Some "pauses" in the code are there just to let Prairie work on its own time. If you remove them, the connection with matlab may be compromised.

For the first image, it is recommendable that you obtain a good quality image, the code will use it to identify rois.
So use an averaged image (max 128) with enough power laser. Then run the first 3 sections until you reach the template matching.
Remember to go back to "low and fast" quality before running the baseline.
If the template matching doesn't do a good job you should change tset.roi parameters. The template matching is an algorithm from Ohki et al. 2005.
The code allows you to remove or add neurons manually.
Then you run the baseline.
The script will then plot some of the more variable neurons for you to choose which ones belong to the BMI.
Chose some (2 and 2 works), and put them on:
```
E1_base = sort([XX YY], 'ascend'); 
E2_base = sort([ZZ WW], 'ascend'); % 
```
Where XX, YY, ZZ, WW are the neurons you have selected as direct neurons for the BMI.
The code will then run the calibration of the BMI based on the data collected during baseline. See calibration files for more info.
After that, it only remains to run the BMI. GOOD LUCK!

 
## Hardware attached (Besides the 2p of course)
 - The code uses an arduino to generate audio as feedback for the mice.
 - Jetball (phenosys) that delivers water reward and records the movement of the mice.
 - Nidaq:
     - to send a TTL pulse to the DORIC dual light source (blue and UV). Blue gets active with D5 and UV with D3.
     - to send a TTL pulse every time the BMI has obtained and works with a frame. To synchronize with prairie images a posteriori.
 - Dual light source to stim D1Rs.

## Folders and what they contain
- calibration -> All the code you need to calibrate the BMI
- parameters -> parameters needed to run the BMI, if you need to change a parameter, please change it here.
- plots -> self-explanatory.
- rois -> all the code to identify and select rois.
- simulation -> to simulate the BMI a posteriori, recalculate T1, obtain T2, etc., following the same code that would have been used in a real experiment.

## Main files for the BMI
- Script to run step by step the BMI  - > Main_Protocol.m
- Baseline acquitisition (see Main_Protocol.m)  - > Baseline_Acqnvs_Prairie.m
- BMI acquisition (see Main_Protocol.m) - > BMI_Acqnvs_Prairie.m

Please let us know if there is anything missing or if documentation needed!
