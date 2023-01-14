import matplotlib.pyplot as plt
import pyFAI.gui.jupyter
import pyFAI
import fabio
from pyFAI.test.utilstest import UtilsTest
from pyFAI.calibrant import CALIBRANT_FACTORY
from pyFAI.goniometer import SingleGeometry
import math
import numpy as np
import pathlib
from tqdm import tqdm
import warnings
import yaml
import os
import re

figure_size = 10

def get_config(path: str) -> dict:
    """Open a yaml file and return the contents."""
    with open(path) as input_file:
        return yaml.safe_load(input_file)
    
def find_nearest(array, value):
    array = np.asarray(array)
    index = (np.abs(array - value)).argmin()
    return array[index], index
    
def plot_single_slice(ai, input_filepath: str, two_theta_min: int, two_theta_max: int, 
                      number_of_points: int = 1000, number_of_cakes: int = 36000, 
                      image_number: int = 0, v_max: int = 10):
    
    """ Plot an azimuthal angle versus two-theta slice from a diffraction pattern image.
    
    :param input_file: input file name of the diffraction pattern image.
    :param two_theta_min: minimum two-theta value for the slicing of the data.
    :param two_theta_max: maximum two-theta value for the slicing of the data.
    :param number_of_points: number of radial points in two-theta (default is 1000).
    :param number_of_cakes: number of azimuthal points in chi, default gives 0.01 degree resolution (default is 36000).
    :param image_number: diffraction pattern image number to plot.
    :param v_max: colormap data range for intensity plotting.
    
    :return: numpy array containing the azimuthal angles, and the intensity values for each time increment.
    """ 
    
    image_list = sorted(pathlib.Path(input_filepath).glob("*"))
    
    input_file = image_list[image_number]
    
    image = fabio.open(input_file).data
    
    result = ai.integrate2d(image, number_of_points, number_of_cakes, unit="2th_deg")
    
    result_intensity = result.intensity 
    result_radial = result.radial
    result_azimuthal = result.azimuthal
    
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize = (15,10))
    
    image1 = ax1.imshow(result_intensity,
               origin='lower',
               extent=[result_radial.min(), 
                       result_radial.max(), 
                       result_azimuthal.min(), 
                       result_azimuthal.max()],
               cmap='viridis',
               interpolation='none',
               aspect='auto',
               vmax = v_max)
#     plt.colorbar(image1)
    plt.xlabel(r"Scattering angle ${\theta}$ ($^{o}$)", fontsize = 25)   
    plt.ylabel(r"Azimuthal angle ${\chi}$ ($^{o}$)", fontsize = 25);

    two_theta_min_value, two_theta_min_index = find_nearest(result_radial, two_theta_min)
    print("The closest two-theta value = ", two_theta_min_value, ", with index = ", two_theta_min_index)

    two_theta_max_value, two_theta_max_index = find_nearest(result_radial, two_theta_max)
    print("The closest two-theta value = ", two_theta_max_value, ", with index = ", two_theta_max_index)

    result_intensity_cropped = result_intensity[:,two_theta_min_index:two_theta_max_index]
    result_radial_cropped = result_radial[two_theta_min_index:two_theta_max_index]
    
    image2 = ax2.imshow(result_intensity_cropped,
               origin='lower',
               extent=[result_radial_cropped.min(), 
                       result_radial_cropped.max(), 
                       result_azimuthal.min(), 
                       result_azimuthal.max()],
               cmap='viridis',
               interpolation='none',
               aspect='auto',
               vmax = v_max)
#     plt.colorbar(image2)
    plt.xlabel(r"Scattering angle ${\theta}$ ($^{o}$)", fontsize = 25)   
    plt.ylabel(r"Azimuthal angle ${\chi}$ ($^{o}$)", fontsize = 25);

    fig.tight_layout()

def get_intensity_time(ai, input_filepath: str, output_filepath:str, lattice_planes: list, 
                       two_theta_min: list, two_theta_max: list, 
                       number_of_points: list, number_of_cakes: list, image_step: int): 
    
    """ Get the intensity data around a ring for a given two-theta slice from a set of diffraction pattern 
    images, and save as a text file.
    
    :param ai: pyFAI detector calibration.
    :param input_filepath: input filepath containing the diffraction pattern images.
    :param output_filepath: output filepath for saving the analysis text file.
    :param lattice_planes: list of hkil or hkl indices of the lattice planes.
    :param two_theta_min: list of minimum two-theta value for the slicing of the data.
    :param two_theta_max: list of maximum two-theta value for the slicing of the data
    :param number_of_points: list of number of radial points in two-theta (default is 1000).
    :param number_of_cakes: list of number of azimuthal points in chi, default gives 0.01 degree resolution (default is 36000).
    
    :return: dictionary of numpy arrays containing the azimuthal angles, and the intensity values for each time increment.
    """ 
    # get a list of the files
    image_list = sorted(pathlib.Path(input_filepath).glob("*"))

    # define dictionary for storing data
    result_dict = dict()
    
    # cake data over multiple images for multiple lattice planes
    for count, lattice_plane in enumerate(lattice_planes):
        start = True
        
        for image_path in tqdm(image_list[0::image_step]):

            # create an image array and cake the data
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                image = fabio.open(image_path)
            pattern_image_array = image.data
            result2d = ai.integrate2d(pattern_image_array,
                                      number_of_points,
                                      number_of_cakes,
                                      unit="2th_deg",
                                      polarization_factor=0.99,
                                      method='full_csr')
            # get the caked results
            result_intensity = result2d.intensity 
            result_radial = result2d.radial
            result_azimuthal = result2d.azimuthal

            # slice the bounds of a ring using the two-theta min and max values
            two_theta_min_value, two_theta_min_index = find_nearest(result_radial, two_theta_min[count])
            two_theta_max_value, two_theta_max_index = find_nearest(result_radial, two_theta_max[count])
            result_intensity_sum = result_intensity[:,two_theta_min_index:two_theta_max_index].sum(axis=1)

            if start:
                # create an array to store the results
                result_array = np.array([result_azimuthal,result_intensity_sum])
                start = False

            else:
                # then stack the results in the array
                result_array = np.vstack((result_array,result_intensity_sum))

        # check output folder exists
        output_folder = f"{output_filepath}/intensity-time-analysis/"
        CHECK_FOLDER = os.path.isdir(output_folder)

        if not CHECK_FOLDER:
            os.makedirs(output_folder)
            print("Created folder : ", output_folder)
        
        # write out the data to a text file
        result_array = np.transpose(result_array)
        np.savetxt(f"{output_folder}intensity_time_{lattice_plane}.txt", result_array)
        
        # store the result array in the result dictionary
        result_dict[lattice_plane] = result_array
        
        # empty the result array for the next iteration
        result_array = np.zeros(np.shape(result_array))
        
    return (result_dict)

def plot_intensity_time(result_dict: dict, lattice_plane_number: int, 
                       lattice_planes: list, v_max, 
                       acquisition_frequency: int, image_step: int, 
                       output_filepath: list, test = None):

    # define the lattice plane
    lattice_plane = lattice_planes[lattice_plane_number]

    # result dictionary contains azimuthal values in first column, and intensity values with time in subsequent columns
    intensity_array = result_dict[lattice_plane][:,1:]
    azimuthal_array = result_dict[lattice_plane][:,0]
    
    time_step = image_step / acquisition_frequency
    time = np.arange(0,np.size(intensity_array,1)) * time_step
    
    if isinstance(v_max, int):
        print("v_max is a single integer value of: ", v_max)
        vmax = v_max
    elif isinstance(v_max, list):
        print("v_max is a list object.")
        vmax = v_max[lattice_plane_number]
    
    plt.rc('xtick', labelsize = 30)
    plt.rc('ytick', labelsize = 30)
    plt.rc('legend', fontsize = 30)
    plt.rc('axes', linewidth = 2)
    plt.rc('xtick.major', width = 2, size = 15)
    plt.rc('xtick.minor', width = 2, size = 10)
    plt.rc('ytick.major', width = 2, size = 15)
    plt.rc('ytick.minor', width = 2, size = 10)

    plt.figure(figsize = (15,10))

    plt.imshow(intensity_array,
               origin='lower',
               extent=[time.min(), 
                       time.max(), 
                       azimuthal_array.min(), 
                       azimuthal_array.max()],
               cmap='viridis',
               interpolation='none',
               aspect='auto',
               vmin = 0,
               vmax = vmax)
    plt.colorbar()
    plt.title(f"{lattice_plane} Intensity-Time", fontsize = 35)
    plt.xlabel(r"Time (s)", fontsize = 35)   
    plt.ylabel(r"Azimuthal angle ${\chi}$ ($^{o}$)", fontsize = 35);

    if not test:
        # check output folder exists
        output_folder = f"{output_filepath}/intensity-time-analysis/"
        CHECK_FOLDER = os.path.isdir(output_folder)

        if not CHECK_FOLDER:
            os.makedirs(output_folder)
            print("Created folder : ", output_folder)      
              
        plt.tight_layout()
        plt.savefig(f"{output_folder}intensity_time_{lattice_plane}.png")
        print(f"Figure saved to: {output_folder}intensity_time_{lattice_plane}.png")

    else:
        print("Image has not been saved because 'test' has been specified.")
        
def save_intensity_time(result_dict, lattice_planes, v_max, 
                        acquisition_frequency, image_step, output_filepath):

    for count, lattice_plane in enumerate(lattice_planes):
        lattice_plane_number = count
        plot_intensity_time(result_dict, lattice_plane_number, lattice_planes, v_max, 
                            acquisition_frequency, image_step, output_filepath)
        
def load_intensity_time_data(output_filepath: str):
    """ Load intensity versus time profiles for different lattice planes from a series of text files.
    
    :param output_filepath: Output filepath to saved text files.
    
    :return: dictionary of numpy arrays containing the azimuthal angles, 
    and the intensity values for each time increment,
    as well as list of lattice planes, which can be used as dictionary keys.
    """
    
    input_folder = f"{output_filepath}/intensity-time-analysis/"
    CHECK_FOLDER = os.path.isdir(input_folder)

    if not CHECK_FOLDER:
        print("Input filepath does not exist: ", input_folder)
        return

    # get a list of the files
    file_list = sorted(pathlib.Path(input_folder).glob("intensity_time_*.txt"))

    # define dictionary and list for storing data
    result_dict = dict()
    lattice_planes = []

    for file_name in file_list:
        # extract lattice plane from file name
        lattice_plane = re.findall(r"intensity_time_(.*)", file_name.stem)
        # remove square bracket and quote from findall and convert to string
        lattice_plane = str(lattice_plane)[2:-2]
        lattice_planes.append(lattice_plane)

        # load array and store in dictionary
        result_array = np.loadtxt(file_name)
        result_dict[str(lattice_plane)] = result_array
   
    return result_dict, lattice_planes