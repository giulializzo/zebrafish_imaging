Complete Image Acquisition, Fish Detection, and Analysis Pipeline
=================================================================

This document provides a comprehensive overview of the process for acquiring, detecting, and analyzing images of zebrafish larvae using the MolecularDevices ImageXpress Confocal microscope.

Software Requirements
---------------------
- MetaXpress Software Version: 6.7.2
- python Version 3.6.4

Author
------
- Guillaume Jacot (guillaume.jacot@rd.nestle.com), 2023

Overview
--------

1. **Preparation of Fish Larvae:**
   - Fish larvae at 3 days post-fertilization (dpf) are placed 3 by 3 in each well of a 96-well plate, totaling 288 fish per plate.
   - The larvae are anesthetized to limit their movement.

2. **Initial Setup:**
   - The plate is positioned inside the MolecularDevices ImageXpress Confocal microscope.

3. **4x Magnification Acquisition:**
   - A 4x magnification acquisition of the plate is conducted, capturing 4 sites per well.
   - Concurrently, the Python script "FishDetection.py" executes on another (or the same) computer to detect fish in a 2x2 montage of the well images.

4. **Fish Detection and ROI Identification:**
   - During the 4x acquisition, the detection code identifies fish and determines the region of interest (ROI).
   - These regions (3 per well, 1 per fish) are documented in an INI file generated upon completion of the code execution.

5. **Completion of 4x Acquisition:**
   - The 4x acquisition concludes when the entire plate is imaged, while the Python code processes the final images.

6. **Transition to 20x Acquisition:**
   - A monitoring loop running on the microscope checks for the presence of the INI file.
   - Once the INI file is accessible, the 20x acquisition is initiated.

7. **20x Magnification Acquisition:**
   - The acquisition positions at 20x are determined by the ROI specified in the INI file.
   - A Z-stack is acquired through the muscle of the larvae to image autophagosomes.

8. **Image Analysis:**
   - Upon completion of the 20x acquisition, the analysis script "FishAnalysis.jzp" can be executed to detect autophagosomes.

Files and Scripts
-----------------
- **Acquisition Setup and Procedures:** "FishAcquisition.jzp"
  - Contains multiple files and README that must be reviewed to understand their usage.
- **Python Detection Code:** "FishDetection.py"
- **Analysis Script:** "FishAnalysis.jzp"
  - Contains multiple files and README that must be reviewed to understand their usage.

Additional Information
----------------------
For detailed information on the Python and package requirements, please refer to the `FishDetection_environment.yml` file included in this repository.