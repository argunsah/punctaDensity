# punctaDensity
Puncta density estimation code for Gesuita et al., 2022, Cell Reports.

%%%%%

For synaptic boutons counting, the mSyp::GFP channel was filtered using the Wiener filtering followed by bilateral filtering for denoising (the Wiener filter performs adaptive noise removal through low-pass filtering of the image that has been corrupted by a stationary noise (Lim, 1990); the bilateral filter smoothens the signal while preserving the edges (Sokullu et al., 2020; Tomasi and Manduchi, 1998)). Filtered images were segmented using Otsu thresholding (Otsu, 1979). Connected puncta were separated using distance transformation based watershed algorithm (Cuisenaire, 1999). The tdTomato channel was used for tracing processes using accurate fast marching (Hassouna and Farag, 2007; Kroon Dirk-Jan, 2021). Tracing was performed both from source to sink and from sink to source to improve the centerline tracing. We calculated the puncta density by finding puncta profiles falling onto 3D centerline profiles (Figure S1G). The distribution of puncta sizes is plotted in Figure S1H. Puncta smaller than 0.5µm (yellow column) were removed from the total count as they appear as background signal. Any observation classified as a suspected outlier according to the interquartile range criterion (red columns) corresponds to a cluster of puncta too close to be distinguishable by the program. In order to be able to use these data points and estimate the number of puncta composing the cluster, we divided each outlier by the mean size of true detected puncta (white columns). Puncta analysis was performed in MATLAB (MathWorks); 

%%%%%
How to Run:

Download Bioformats Matlab Library and Put in this same folder with this code
(https://downloads.openmicroscopy.org/bio-formats/5.3.4/artifacts/bfmatlab.zip)

Run punctaDensity.m

Part 1 and 2 reads oif files and converts them in matlab format after extracting some features.
Part 3 is the manual tracing part. 

- When contrast enhancement windows pops up, adjust accordingly and click X to close (Don't click Adjust).
- Click end of one of the processes. Use 'z' to zoom, 'x' to unzoom. Left click multiple points along the process. Right click determines the last point.
- If you are happy with tracing, click 'Yes' to accept, 'No' to reject, 'Done' to complete the tracing step,

Part 4 is to read all files to convert R format.

Please leave a comment in the Issues section if you want to use the code but couldn't figure out how to use it exactly.
