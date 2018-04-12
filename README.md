# Angle_SelfCalibration
Self calibration algorithms for illumination angle recovery.

Illumination angle self-calibration in Fourier ptychography
Matlab code implementation of circle-finding algorithm to calibrate the angles of illumination for brightfield images 
and extend the calibration to darkfield images. Circular edge detection image processing is done on the Fourier transform 
of the data to calibrate the brightfield illumination directly. Python code implementation of spectral correlation 
calibration method inside the Fourier ptychography algorithm. The angles are iteratively refined within the FPM algorithm 
by correlating overlapping spectra. Sample data is provided.

Please cite as:
R. Eckert, Z.F. Phillips, L. Waller. Efficient illumination angle self-calibration in Fourier ptychography. Manuscript submitted for publication. (2018) 

Available in preprint at: https://arxiv.org/abs/1804.03299.
