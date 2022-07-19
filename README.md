# punctaDensity
Puncta density estimation code for Gesuita et al., 2022, Cell Reports.

Download Bioformats Matlab Library and Put in this same folder with this code
(https://downloads.openmicroscopy.org/bio-formats/5.3.4/artifacts/bfmatlab.zip)

Run punctaDensity.m

Part 1 and 2 reads oif files and converts them in matlab format after extracting some features.
Part 3 is the manual tracing part. 

- When contrast enhancement windows pops up, adjust accordingly and click X to close (Don't click Adjust).
- Click end of one of the processes. Use 'z' to zoom, 'x' to unzoom. Left click multiple points along the process. Right click determins the last point.
- If you are happy with tracing, click 'Yes' to accept, 'No' to reject, 'Done' to complete the tracing step,

Part 4 is to read all files to convert R format.
