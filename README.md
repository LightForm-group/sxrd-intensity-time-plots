sxrd-intensity-time-plots
-----------

A Python notebook for creating intensity-time plots (intensity versus time and azimuthal angle) from a dataset of synchrotron X-ray diffraction (SXRD) pattern images, for multiple individual lattice plane rings, using the pyFAI and FabIO packages. 

Intensity-time plots show how the variation in intensity around each lattice plane ring changes over time, throughout an experiment, which can be used to discern in-situ qualitative evolution of material properties, such as changes in crystallographic texture and activation of continuous/discontinuous recrystallization. 

The package works with diffraction pattern image data in the form of .cbf or .tiff images, and uses a .poni calibration file to calibrate the data. See the [pyFAI-integration-caking](https://github.com/LightForm-group/pyFAI-integration-caking) package for creating a .poni calibration file from a calibrant diffraction pattern image.

The intensity-time plots are saved as .png images, along with .txt files containing the intensity values for replotting.

Development
--------------

This package was developed by Christopher S. Daniel at The University of Manchester, UK, and was funded by the Engineering and Physical Sciences Research Council (EPSRC) via the LightForm programme grant (EP/R001715/1). LightForm is a 5 year multidisciplinary project, led by The Manchester University with partners at University of Cambridge and Imperial College, London (https://lightform.org.uk/).

Contents
-----------

**It is recommended the user works through an example in the notebook:**
    
`intensity-time-plots.ipynb` Contains an example using inputs from a .yaml text file demonstrating the production of intensity-time plots for plotting changes in intensity variation around multiple individual diffraction pattern rings over time.

*Note, the `example-data/`, `example-calibration/` and `example-analysis/` folders contain data that can be used as an example analysis, but a clear external file structure should be setup to support the analysis of large synchrotron datasets.*

Installation and Virtual Environment Setup
-----------

Follow along by copying / pasting the commands below into your terminal (for a guide on using a python virtual environments follow steps 4-7).

**1. First, you'll need to download the repository to your PC. Open a unix command line on your PC and navigate to your Desktop (or GitHub repository):**
```unix
cd ~/Desktop
```
**2. In your teminal, use the following command to clone this repository to your Desktop:**
```unix
git clone https://github.com/LightForm-group/sxrd-intensity-time-plots.git
```
**3. Navigate inside `Desktop/sxrd-intensity-time-plots/` and list the contents using `ls`(macOS) or `dir`(windows):**
```unix
cd ~/Desktop/sxrd-intensity-time-plots/
```
**4. Next, create a python virtual environment (venv) which contains all of the python libraries required to use sxrd-intensity-time-plots.
Firstly, use the following command to create the venv directory which will contain the necessary libraries:**
```unix
python -m venv ~/Desktop/sxrd-intensity-time-plots/venv
```
**5. Your `sxrd-intensity-time-plots/` directory should now contain `venv/`. Install the relevant libraries to this venv by first activating the venv:**
```unix
source ~/Desktop/sxrd-intensity-time-plots/venv/bin/activate
```
*Note, the beginning of your command line should change from `(base)` to include `(venv)`.*

**6. Install the python libraries to this virtual environment using pip and the `requirements.txt` file included within the repository:**
```unix
pip install -r ~/Desktop/sxrd-intensity-time-plots/requirements.txt
```
**7. To ensure these installed correctly, use the command `pip list` and ensure the following packages are installed:**
```unix
pip list
# Check to ensure that all of the following are listed:
#notebook
#pyFAI
#fabio
#numpy
#tqdm
#pyyaml
```
**8. If all in step 7 are present, you can now run the notebook.
Ensure the venv is active and use the following command to boot jupyter notebook (using all libraries installed in the venv)
Warning - using just `jupyter notebook` without `python -m` can result in using your default python environment (the libraries may not be recognised):**
```unix
python -m jupyter notebook
```
**9. Work through the notebook and setup yaml text files for reproducible intensity-time plots from multiple individual diffraction pattern rings across large synchrotron datasets.**

**10. When you're finished using the virtual environment, deactivate it!
This will avoid confusion when using different python libraries that are not installed within the virtual environment:**
```unix
deactivate
```

Required Libraries
--------------------

The required libraries are listed in requirements.txt