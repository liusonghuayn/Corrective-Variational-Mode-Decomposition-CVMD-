MATLAB code for **Corrective Variational Mode Decomposition (CVMD)**, an adaptive extension of VMD that automatically selects the optimal number of modes for robust time–frequency signal decomposition.

This repository contains:
- `CVMD.m` – main implementation of CVMD (returns optimal mode number and IMFs).
- `VMD.m` – standard Variational Mode Decomposition implementation used inside CVMD.
- `maincode.m` – example script demonstrating how to apply CVMD to a test signal.
- `platCVMD_IMFs.m` – helper function to plot the decomposed IMFs.

If you use this code in your research, please cite our paper:

> Liu, S., Lang, X., Wu, J., Zhang, Y., Lei, C., & Su, H.  
> **“Corrective variational mode decomposition to detect multiple oscillations in process control systems”**,  
> *Control Engineering Practice*, 2025, 154: 106123.
