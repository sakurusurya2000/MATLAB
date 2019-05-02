Note: EMD source code can be download from 
http://perso.ens-lyon.fr/patrick.flandrin/emd.html. 

If you have any problems in running this code, please contact Heng Li via wudaomana@gmail.com.

Main Files:

-	tvfemd.m
	Computes the TVF-EMD.

-	splinefit.m
	Computes the B-spline least squares approximation (B-spline order = 26). 

Optional Files:

-	INST_FREQ_local.m:
	This function computes the Hilbert-Huang spectrum using the Hilbert transform, the instantaneous amplitude and the instantaneous 
	
-	spectrogram_emd.m
	transforms the 2D instantaneous amplitude and frequency into 3D spectrum matrix 
	 to plot the instantaneous amplitude contours on the time-frequency plane. 

-	disp_hhs.m 
	plots the 3D HHS spectrum of the spectrum matrix obtained as an output from spectrogram_emd.m.

Example Files:


-	examples\signal_decomposition.m : muti-component signal decomposition using TVF-EMD				    