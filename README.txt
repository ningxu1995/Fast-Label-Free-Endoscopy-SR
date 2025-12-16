GENERAL INFORMATION

1. Title of Dataset: 
   Source Code for Iterative Phase Retrieval and Pattern Generation in Super-Resolution Structured Illumination Microscopy

2. Author Information
   A. Principal Investigator / Corresponding Author
      Name: Ning Xu and Sarah E. Bohndiek
      Institution: Department of Physics, Cavendish Laboratory, University of Cambridge, JJ Thomson Avenue, Cambridge, CB3 0HE, UK; State Key Laboratory of Precision Measurement Technology and Instruments, Department of Precision Instrument, Tsinghua University, Beijing 100084, China;  Cancer Research UK Cambridge Institute, University of Cambridge, Robinson Way, Cambridge, CB2 0RE, UK
      Email: nxu@mail.tsinghua.edu.cn; seb53@cam.ac.uk

3. Date of Collection: 
   2021-2022

4. Geographic Location of Data Collection: 
   Cambridge, UK & Beijing, China 

-----------------------------------------------------------------------------------------

SHARING/ACCESS INFORMATION

1. Licenses/restrictions placed on the data: 
   Open access for academic research purposes. Attribution required.

2. Recommended citation for this dataset: 
   Xu, N., et al. (2025). Code for Fast label-free point-scanning super-resolution imaging for endoscopy. http://arxiv.org/abs/2512.13432.

-----------------------------------------------------------------------------------------

DATA & FILE OVERVIEW

1. File List: 
   A. SRspotN_N_Optimized.m
      - Implementation of the Weighted Gerchberg-Saxton (GS) algorithm.
      - Used to generate the phase mask for high-density multi-spot arrays (10x10) approaching the diffraction limit.
   
   B. Get_Phase_Estimation.m
      - Algorithm for estimating raw phase shifts in point-by-point intensities images using cross-correlation in the Fourier domain.
      - Includes sparsity checks for cell samples.

   C. Get_System_Parameters.m
      - Configuration file defining the optical system constraints (NA, wavelength, pixel size) and generating the Optical Transfer Function (OTF).

   D. Sequence_Phase_Generator.m
      - Generates the final grating patterns for the Spatial Light Modulator (SLM) by combining the GS-calculated phase with linear phase ramps for scanning.

2. Methodological Notes:
   - The phase retrieval logic (SRspotN_N) utilizes a modified IFTA (Iterative Fourier Transform Algorithm) with an amplitude constraint at the Fourier plane to ensure spot uniformity.

3. Software-specific information:
   - Platform: MATLAB R2021b or later.
   - Dependencies: Image Processing Toolbox, Signal Processing Toolbox.