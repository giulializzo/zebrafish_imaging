FishDetection Code Walkthrough
================

The FishDetection.py script is designed to be initiated before the live acquisition of fish at 4x magnification using the ImageXpress confocal microscope. It processes a set of TIFF images to detect fish objects and extract specific locations, which are then documented in an INI file for subsequent 20x acquisition with the same microscope. Below is a step-by-step explanation of the process:

Step 1: Drive Mapping
----------------------
- The function `check_mappedDrive` ensures that the necessary network drives are mapped correctly.
- It checks for the MolecularDevices ImageXpress Confocal image storage location as well as location to store the resulting INI file.
- The paths in this function are explicit and need to be modified according to the user's environment, specifying where to look for images and where to store the resulting INI file.

Step 2: Image Collection
------------------------
- The function `fish_detection` gathers four TIFF images from a specified directory.
- These images are combined to create a single montage for processing.

Step 3: Montage Creation
------------------------
- Using the `montage` function from the `skimage.util` module, the four images are merged into a single montage image.
- This montage serves as the basis for object detection.

Step 4: Object Detection and Filtering
--------------------------------------
- The montage image is processed to detect objects using various image processing techniques, including Gaussian filtering and thresholding.
- Detected objects are filtered based on area criteria to identify potential fish.

Step 5: Skeletonization
-----------------------
- The function `walk_the_line` is employed to skeletonize the identified fish objects.
- This involves tracing the skeleton of each fish to determine its structure.

Step 6: Position Extraction
---------------------------
- The `walk_the_line` function further traverses the skeleton to find specific positions on each fish.
- These positions are crucial for subsequent analysis and documentation.

Step 7: INI File Creation
-------------------------
- The function `write_ini_file` writes the extracted position information into an INI file.
- This file serves as a record of the detected fish and their respective positions, which will be read during the subsequent 20x acquisition.

Additional Information
----------------------
For detailed information on the Python and package requirements, please refer to the `FishDetection_environment.yml` file included in this repository.