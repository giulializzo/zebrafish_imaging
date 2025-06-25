# -*- coding: utf-8 -*-

"""
This script processes TIFF images to create montages, detect fish objects, and extract specific locations.
It writes the extracted information into an INI file for further analysis and documentation.

Find the requirements in environment.yml

Author:
Guillaume Jacot (guillaume.jacot@rd.nestle.com), 2023
"""
import os
import glob
import yaml
from collections import defaultdict
import math
import sys
import warnings
warnings.filterwarnings("ignore")

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import numpy as np
import pandas as pd

from skimage.util import montage
from skimage import filters
from skimage import morphology, measure
from skimage.exposure import equalize_adapthist as clahe
from skan import csr
from scipy import spatial, ndimage
import imageio as iio
import time
from timeit import default_timer as timer
try:
    from configparser import ConfigParser
except ImportError:
    from ConfigParser import ConfigParser  # ver. < 3.0
import re
import warnings
import datetime

# =============================================================================
# Functions
# =============================================================================

def get_date():
    """
    Get the current date in YYYY-MM-DD format.

    Returns:
        str: Current date in YYYY-MM-DD format.
    """
    now = datetime.datetime.now()
    return now.strftime("%Y-%m-%d")

def collect_image_metadata(image, path):
    """
    Collect metadata from the image and its file path.

    Args:
        image: Image object from which metadata is to be extracted.
        path (str): Path of the image file.

    Returns:
        dict: Dictionary containing metadata.
    """
    md = dict()# Initialize an empty dictionary to store metadata

    # Regular expression to extract PlateID and Barcode from the file path
    searchObj = re.search( r'.+4x[\\](.*)[\\]\d{4}-\d{2}-\d{2}[\\](\d+)[\\].*$', path, re.M|re.I)

    if searchObj:
        # Extract PlateID and Barcode if the regular expression matches
        md['PlateID'] = searchObj.group(2)
        bc = searchObj.group(1)
    else:
        # Alternative regular expression if the first one fails
        searchObj_2 = re.search( r".+4x[\\]\d{4}-\d{2}-\d{2}[\\](\d+)[\\].*$", path, re.M|re.I)
        if searchObj_2:
            # Extract PlateID if the alternative regular expression matches
            md['PlateID'] = searchObj_2.group(1)
            bc = ""
        else:
            print("Nothing found!!")
            md['PlateID'] = 'None'
            bc = ""

    # Clean and split the description metadata
    dsc = (image.meta['description']
           .replace('&#13','')
           .replace('&#10','')
           .replace('\n<prop',';')
          )

    lst=dsc.split(';')  # Split the description into a list of properties

    for prop in lst:
        if prop.find('Barcode') > -1:
            md['barcode'] = prop[8:]

        newSearchObj = re.search( r'(.*\d+_.*)$',prop)
        if newSearchObj:
            md['barcode'] = newSearchObj.group(0)

        if prop.find('stage-label') > -1:
            temp = prop[prop.find('value')+7:len(prop)-3]
            md['stage_label'] = temp.rsplit(' : ')[0]

        if prop.find('spatial-calibration-x') > -1:
            md['cal_x'] = float(prop[prop.find('value')+7:len(prop)-3])

        if prop.find('spatial-calibration-y') > -1:
            md['cal_y'] = float(prop[prop.find('value')+7:len(prop)-3])

        if prop.find('stage-position-x') > -1:
            md['stage_x'] = float(prop[prop.find('value')+7:len(prop)-3])

        if prop.find('stage-position-y') > -1:
            md['stage_y'] = float(prop[prop.find('value')+7:len(prop)-3])

        if prop.find('pixel-size-x') > -1:
            md['size_x'] = int(prop[prop.find('value')+7:len(prop)-3])

        if prop.find('pixel-size-y') > -1:
            md['size_y'] = int(prop[prop.find('value')+7:len(prop)-3])

        #print('------')

    # Calculate origin coordinates
    md['origin_x'] = md['stage_x']-(md['size_x']/2)*md['cal_x']
    md['origin_y'] = md['stage_y']-(md['size_y']/2)*md['cal_y']

    # Set barcode to 'None' if not found
    if not md.get('barcode'):
        if len(bc) > 0:
            md['barcode']= bc
        else:
            md['barcode']= 'None'

    return md

def convert_to_absolute_position(row, col, metadata):
    """
    Convert image coordinates to absolute stage coordinates.

    Args:
        row (int): Row coordinate in the image.
        col (int): Column coordinate in the image.
        metadata (dict): Metadata dictionary containing calibration and origin information.

    Returns:
        tuple: Absolute x and y coordinates.
    """
    x = np.round(metadata['origin_x'] + col * metadata['cal_x'], 1)
    y = np.round(metadata['origin_y'] + row * metadata['cal_y'], 1)
    return x, y

def distanceToZero(df, center):
    """
    Calculate Euclidean distance from a point to the center.

    Args:
        df (DataFrame): DataFrame containing coordinates.
        center (tuple): Center coordinates.

    Returns:
        float: Euclidean distance.
    """
    u = (df['image-coord-src-0'], df['image-coord-src-1'])
    return spatial.distance.euclidean(u, center)

def distanceToOne(df, center):
    """
    Calculate Euclidean distance from a point to the center.

    Args:
        df (DataFrame): DataFrame containing coordinates.
        center (tuple): Center coordinates.

    Returns:
        float: Euclidean distance.
    """
    u = (df['image-coord-dst-0'], df['image-coord-dst-1'])
    return spatial.distance.euclidean(u, center)

def findStartingPoint(bd, center):
    """
    Find the starting point for skeleton traversal.

    Args:
        bd (DataFrame): DataFrame containing skeleton data.
        center (tuple): Center coordinates.

    Returns:
        tuple: Starting node and its coordinates.
    """
    center = center[1][0]  # Extract the center coordinate
    df = pd.DataFrame() # Initialize a DataFrame to store distances

    # Calculate Euclidean distance to the center for each node
    df['EdToZero'] = bd.apply(distanceToZero, args=(center,), axis=1)
    df['EdToOne'] = bd.apply(distanceToOne, args=(center,), axis=1)

    # Find the maximum distance and its index
    df['maxEd'] = df.values.max(1)
    maximum = df['maxEd'].max()
    idx = df['maxEd'].idxmax()

    # Determine the starting node based on the maximum distance
    if maximum in df['EdToZero'].unique():
        return (bd['node-id-0'].iloc[idx],
                bd['image-coord-src-0'].iloc[idx],
                bd['image-coord-src-1'].iloc[idx]
                )
    else:
        return (bd['node-id-1'].iloc[idx],
                bd['image-coord-dst-0'].iloc[idx],
                bd['image-coord-dst-1'].iloc[idx]
                )

def GetAngle (l1, l2):
    ''' compute angle (in degrees) between two lines
    Inputs:
        l1, l2 - lines in the form [[x0,y0],[x1,y1]]
    '''
    v1 = np.array([l1[1][0]-l1[0][0],l1[1][1]-l1[0][1]])
    v2 = np.array([l2[1][0]-l2[0][0],l2[1][1]-l2[0][1]])

    v1_norm = np.linalg.norm(v1)
    v2_norm = np.linalg.norm(v2)
    scalar = np.dot(v1,v2)
    angle = 180*np.arccos(scalar/(v1_norm*v2_norm))/math.pi
    return angle
				
def getNextNode(next_indices, current_node, previous_nodes, coords):
    """
    Determine the next node to move to in a skeletonized image based on the angle formed with the previous nodes.

    Args:
        next_indices (list): List of indices of the next possible nodes.
        current_node (int): Index of the current node.
        previous_nodes (list): List of indices of the previously visited nodes.
        coords (ndarray): Array of coordinates of the nodes.

    Returns:
        int: Index of the next node to move to.
    """
    # Determine the reference line l0 based on the previous nodes
    if len(previous_nodes) >= 22:
        l0 = [list(coords[previous_nodes[-21]].astype(int)),
              list(coords[current_node].astype(int))]
    else:
        l0 = [list(coords[previous_nodes[0]].astype(int)),
              list(coords[current_node].astype(int))]

    angles = pd.Series()

    # Calculate the angle for each possible next node
    for index in next_indices:
        if index not in [current_node, previous_nodes[-2]]:
            row = coords[index][0].astype(int)
            col = coords[index][1].astype(int)
            lx = [list(coords[index].astype(int)),
                  list(coords[current_node].astype(int))]
            angle = GetAngle(l0, lx)
            temp = pd.Series(angle, index=(index,))
            angles = angles.append(temp)

    # Return the index of the node that forms the largest angle
    return angles.idxmax()

def walk_the_line(region_label, majorAxisLength, label, metadata):
    """
    Traverse the skeleton of a segmented fish image to find the end point.

    Args:
        region_label (int): Label of the region to be processed.
        majorAxisLength (float): Major axis length of the region.
        label (ndarray): Labeled image array.
        metadata (dict): Metadata dictionary containing calibration and origin information.

    Returns:
        ndarray: Coordinates of the end point of the skeleton.
    """
    # Isolate the single fish region
    single_fish = (label == region_label)

    # Compute the centroid of the region
    center = csr.compute_centroids(single_fish)

    # Skeletonize the region
    skeleton = morphology.skeletonize(single_fish)

    # Summarize the skeleton
    bd = csr.summarise(skeleton)

    # Convert the skeleton to a graph
    graph, coords, degrees = csr.skeleton_to_csgraph(skeleton,
                                                     spacing=(metadata['cal_x'],
                                                              metadata['cal_y'])
                                                     )

    # Find the starting point for traversal
    start_node, start_row, start_col = findStartingPoint(bd, center)

    # Initialize traversal variables
    points_list = [start_node]
    neigh, count, distance = 2, 0, 0
    next_node = start_node
    stop = False

    # Traverse the skeleton
    while not ((neigh == 1) or
               (count > 5000) or
               (distance > majorAxisLength * metadata['cal_x'] / 3) or
               stop):
        indices = graph.getrow(next_node).indices
        data = graph.getrow(next_node).data
        i = 0

        if indices.shape[0] <= 2:
            for index in indices:
                if index not in points_list:
                    next_node = index
                    points_list.append(index)
                    row = coords[next_node][0].astype(int)
                    col = coords[next_node][1].astype(int)
                    neigh = degrees[row, col]
                    distance += data[i]
                    i += 1
        else:
            next_node = getNextNode(indices, next_node, points_list, coords)
            points_list.append(next_node)
            row = coords[next_node][0].astype(int)
            col = coords[next_node][1].astype(int)
            neigh = degrees[row, col]
        count += 1

    return coords[next_node]

def write_ini_file(roi_row, roi_col, configfile, metadata):
    """
    Write the region of interest (ROI) information to an INI file.

    Args:
        roi_row (int): Row coordinate of the ROI.
        roi_col (int): Column coordinate of the ROI.
        configfile (ConfigParser): ConfigParser object to write the INI file.
        metadata (dict): Metadata dictionary containing calibration and origin information.
    """
    try:
        # Add 'Experiment' section if it doesn't exist
        configfile.add_section('Experiment')
        configfile.set('Experiment', '4x-ID', metadata['PlateID'])
    except:
        pass

    sl = metadata['stage_label']

    # Add stage label section if it doesn't exist
    if not configfile.has_section(sl):
        configfile.add_section(sl)
        reg_num = 0
    else:
        reg_num = configfile.getint(sl, 'region_number')

    # Update the region number and origin coordinates
    configfile.set(sl, 'region_number', str(reg_num + 1))
    configfile.set(sl, 'OriginX', "{0:.2f}".format(metadata['origin_x']))
    configfile.set(sl, 'OriginY', "{0:.2f}".format(metadata['origin_y']))

    # Set the ROI coordinates
    temp = "{}:{}".format(roi_row, roi_col)
    configfile.set(sl, "region_{}_ROIMicro".format(reg_num + 1), temp)

def fish_detection(root_dir, imageList, name, config):
    """
    Detect fish in a set of images, create a montage, process the montage to find fish,
    and write the region of interest (ROI) information to an INI file.

    Args:
        root_dir (str): Root directory containing the images.
        imageList (list): List of image file names.
        name (str): Name for the output files.
        config (ConfigParser): ConfigParser object to write the INI file.

    Returns:
        tuple: Plate ID and barcode extracted from the metadata.
    """
    # =============================================================================
    # Load Images
    # =============================================================================
    well_image = []
	
	#Folder to save montage images for quality control, to be changed accordingly
    out_dir = 'C:' + os.sep + 'Images'

    # Read and append images to the well_image list
    in_path = os.path.join(root_dir, imageList[0])
    raw = iio.imread(in_path)
    metadata = collect_image_metadata(raw, in_path)
    well_image.append(raw)

    in_path = os.path.join(root_dir, imageList[1])
    raw = iio.imread(in_path)
    well_image.append(raw)

    in_path = os.path.join(root_dir, imageList[2])
    raw = iio.imread(in_path)
    well_image.append(raw)

    in_path = os.path.join(root_dir, imageList[3])
    raw = iio.imread(in_path)
    well_image.append(raw)

    out_dir = os.path.join(out_dir, str(metadata['PlateID']))

    try:
        os.mkdir(out_dir)
    except FileExistsError:
        pass

    # =============================================================================
    # Create Montage from the 4 images
    # =============================================================================
    total = montage(well_image)

    # =============================================================================
    # Process montage image to find the fish
    # =============================================================================
    gaussed = filters.gaussian(total, sigma=14, mode='mirror', preserve_range=True)
    raw_mask = gaussed > filters.threshold_triangle(gaussed)
    raw_mask = ndimage.morphology.binary_fill_holes(raw_mask.astype(int))
    raw_mask = morphology.binary_erosion(raw_mask, morphology.disk(7))

    # =============================================================================
    # Measure the ROI found on the processed mask
    # =============================================================================
    label = measure.label(raw_mask)
    regions = measure.regionprops(label, total)

    # =============================================================================
    # Loop through raw objects, select them based on area and other measurements.
    # Once selected, send the image of single fish to be skeletonized and ROI to be found
    # =============================================================================
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.imshow(clahe(total), cmap='gray',
              vmin=np.percentile(clahe(total), 50),
              vmax=np.percentile(clahe(total), 99.5)
              )

    for prop in regions:
        if (prop.area > 35000) and (prop.area < 150000):
            roi_row, roi_col = walk_the_line(prop.label,
                                             prop.major_axis_length,
                                             label,
                                             metadata)
            print("Region {} ({},{})".format(prop.label, roi_row, roi_col))
            circle = mpatches.Circle((roi_col, roi_row), radius=20,
                                     fill=False, edgecolor='red', linewidth=2)
            ax.add_patch(circle)
            r, c = convert_to_absolute_position(roi_row, roi_col, metadata)
            write_ini_file(r, c, config, metadata)

    ax.set_axis_off()
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, name + '.png'))
    plt.close()

    return metadata['PlateID'], metadata['barcode']


	def check_mappedDrive():
    """
    Check if the required network drives are mapped and accessible.

    Returns:
        str: Root directory path if the drives are mapped correctly.

    Raises:
        SystemExit: If the required drives are not mapped.
    """
    # Check if Drive1 (image storage location) is mapped, Drive1 is usally mapped to X:
    if os.path.exists('Drive1:' + os.sep + "ImageStorage" + os.sep + 'Subfolder' + os.sep + 'Data'):
        root_dir = ('Drive1:' + os.sep + "ImageStorage" + os.sep + 'Subfolder' +
                    os.sep + 'Data')
    else:
        s = ("\n****\n"
             "I cannot find this path:\n{}\n"
             "Please MAP the network drive:\n{}\n"
             "to the letter Drive1 (image storage location)\n"
             "****\n").format(
                     'Drive1:' + os.sep + "ImageStorage" + os.sep + 'Subfolder' + os.sep + 'Data',
                     '\\\\network_path\\image_storage'
                     )
        print(s)
        sys.exit()

    # Check if Drive2 (INI file storage location) is mapped, Drive2 is usually mapped to Y:
    if not os.path.exists('Drive2:'):
        s = ("\n****\n"
             "I cannot find this path:\n{}\n"
             "Please MAP the network drive:\n{}\n"
             "to the letter Drive2 (INI file storage location)\n"
             "****\n").format(
                     'Drive2:',
                     '\\\\another_network_path\\ini_file_storage'
                     )
        print(s)
        sys.exit()

    return root_dir
	
def main():
    """
    Main function to process images, detect fish, and write ROI information to an INI file.
    """
    config = ConfigParser()
    config.optionxform = str

    # Get the current date and script directory
    today = get_date()
    script_dir = os.getcwd()
    print(script_dir)

    # Define the processed file path
    processed_file = script_dir + os.sep + today + '_processed.yaml'
    print(processed_file)

    # Check if the required network drives are mapped
    root_dir = check_mappedDrive()
    print(root_dir)

    # Change the current working directory to the root directory
    os.chdir(root_dir)

    processed = set([])
    to_process = defaultdict(set)

    # Load the processed files if the processed file exists
    if os.path.exists(processed_file):
        with open(processed_file, 'r') as fh:
            data = yaml.safe_load(fh)
            processed = set(data['files'])

    os.chdir(root_dir)
    print("Let's go ...")

    # =============================================================================
    # *** PLEASE UPDATE ACCORDING TO YOUR PLATE PLAN ***
    #    Be careful:    odd columns finish at row H
    #                   even columns finish at row A
    # =============================================================================
    last_well = 'A12'
    # =============================================================================
    #
    # =============================================================================

    finish_loop = False

    while not finish_loop:
        new_file = False

        for root, dirs, files in os.walk(root_dir + os.sep + today):
            for f in files:
                f = os.path.join(root[root.find("4x") + 3:], f)
                if (f.endswith(".tif") and
                    f.find("thumb") < 0 and
                    f not in processed):
                    plate_id = f.rsplit(os.sep)[-3]
                    well_id = os.path.basename(f).split("_")[1]
                    group_id = plate_id + "_" + well_id
                    to_process[group_id].add(f)

                    if len(to_process[group_id]) == 4:
                        to_process_list = list(to_process[group_id])
                        to_process_list.sort()
                        print("processing: {}".format(group_id))
                        plateID, barcode = fish_detection(root_dir,
                                                          to_process_list,
                                                          group_id,
                                                          config)
                        for image in to_process[group_id]:
                            processed.add(image)
                            new_file = True
                        to_process.pop(group_id)

                        if well_id == last_well:
                            finish_loop = True

            if new_file:
                with open(os.path.join(script_dir, processed_file), 'w') as fh:
                    yaml.dump({'files': list(processed)}, fh)
            try:
                if barcode.find('None') < 0:
                    toWrite = barcode
                else:
                    toWrite = plateID
            except NameError:
                pass

    with open('W:' + os.sep + '{}.ini'.format(str(toWrite)), 'w') as configfile:
        config.write(configfile)
    print('INI file <{}> is created'.format(toWrite))

if __name__=="__main__":
    main()

