# VBET 2.0
Valley Bottom Extraction Tool (upated: 02/22/2021)

VBET uses a stream network shapefile and digital elevation model to derive a valley bottom polygon based on two lines of evidence: slope and inundation depth. The values for these parameters are adjusted based on the drainage area of a given segment of the network based on the generalizations that at lower drainage areas (small streams), slopes within valley bottoms are generally steeper than at high drainage areas (larger rivers), where slopes tend to be very flat in valley bottoms; and also that high in the network, at low drainage areas, flood depths tend to be shallower than in larger portions of the network.

As such, the user selects to drainage area threshold values that are used to split the network into 'large', 'medium' and 'small' portions. The user then selects threshold values for slope and inundation depth for each of the three portions of the network.

### Known issues
When deriving a stream network from a DEM, segments passing through lakes or resersvoirs sometimes are completely straight with no slope. These segments often give the tool trouble. For now, the best solution is to combine segments so that a segment encompasses more than the flat area, and therefore has some slope, and to manually add vertices and create some small amount of curvature so that it is not perfectly straight. 

If additional issues are encountered, please report them in the GitHub Issues page for VBET 2.0.


## Data preparation
Because the tool relies on drainage area, a drainage area raster is required to run the tool. To be accurate, this raster must encompass the entire watershed upstream of the area for which a valley bottom is being delineated. In cases where a small LiDAR dataset is being used to delineate valley bottoms for a portion of a drainage network, a coarser DEM of the entire basin (e.g. the 10m DEMs of the NED from USGS) should be used to generate this raster. This can be accomplished using the common GIS workflow of filling the pits in the DEM, generating flow directions, and then a flow accumulation raster. The flow accumulation raster can be converted to drainage area in a raster calculator using the equation (flow_accumulation * (flow_accumulation_raster_resolution^2))/1000000.

For the best results, the drainage network should align well with the DEM used. If the network was generated using the DEM, this isn't a problem, but when using external datasets such as NHD, there can be a large spatial disparity between the location of the drainage network in the DEM and in the actual dataset. If this is the case, to the extent possible, the user should edit the network, adjusting it so that it aligns with the network location in the DEM. Additionally, the segment length of the network should be suitable given the DEM resolution. That is, with higher resolution DEMs the network should be segmented at smaller lengths than for coarser DEMs (e.g. 10m). A good length is generally 100-500 times the DEM resolution (e.g. for a 10m DEM, we use 1000m segments, and for a 1m DEM we use 200-500m segments, depending on river size). Stream networks *must contain ONLY SINGLEPART FEATURES. If it contains multipart features*, use a "multipart to singlepart" GIS tool to fix this issue. Finally, to work with the detrending algorithm, network segments *must have 5+ vertices*. Very short, or perfectly straight segments usually don't have enough and the tool will throw an exception, notifying the user of which segments need more vertices.

## Running the model
To run the model, go to the releases tab and download the release appropriate for your OS (Windows 64 bit and Linux 64 bit available). Unzip the download folder in a chosen location. For windows, open the folder and double click on the VBET.exe to run the program. For Linux, open a command prompt inside the folder and type the command: ./VBET to run the program. 

In addition to a shapefile of the valley bottom, the tool produces a text file with the same name as the output with "metadata" appended on the end that contains the files and parameter values used for that particular run of the tool.

If you are using a different OS, or prefer to run the code manually in an IDE, you can download the repository and run it in an IDE. Configure a Python environment with the packages listed below, fill out the parameter values in the 'run_VBET.py' script and then run the script.

### Required Python packages
#### Python 3
- geopandas 
- rasterio 
- shapely
- rasterstats
- numpy 
- scipy
- scikit-image 

## Parameters 

- Stream Network: navigate to the stream network shapefile. This can be e.g. the National Hydrologic Dataset (NHD), extracted from a DEM, or manually delineated by the user. The stream network should be segmented so that each reach is approximately the desired segment length. The network should also be projected into a projected coordinate reference system.
- Digital Elevation Model: navigate to the DEM from which you wish to derive the valley bottom polygon (.tif). The DEM should be in the same projected CRS as the stream network.
- Scratch Work Folder: Navigate to a path to a folder where temporary files will be stored. 
- Drainage Area Raster: navigate to the drainage area raster for the basin (.tif). This should be in the same CRS as the stream network and DEM.
- Output: enter a path and file name to store the valley bottom shapefile produced by the tool.
- Large Drainage Area Threshold: enter a value (int) for a drainage area threshold above which the network will be considered 'large'.
- Medium Drainage Area Threshold: enter a value (int) for a drainage area threshold above which the network will be considered 'medium', and below which the network will be considered 'small'.
- Large Slope Threshold: enter a value for the slope threshold in the large portion of the network. 
- Medium Slope Threshold: enter a value for the slope threshold in the medium portion of the network.
- Small Slope Threshold: enter a value for the slope threshold in the small portion of the network.
- Large Buffer: enter a value representing the maximum distance from the large portion of the network within which valley bottoms occur.
- Medium Buffer: enter a value representing the maximum distance from the medium portion of the network within which valley bottoms occur.
- Small Buffer: enter a value representing the maximum distance from the small portion of the network within which valley bottoms occur.
- Minimum Buffer: enter a value representing a minimum valley bottom width/channel width for the entire network.
- Large Depth: enter a value for the depth threshold of the large portion of the network.
- Medium Depth: enter a value for the depth threshold of the medium portion of the network.
- Small Depth: enter a value for the depth threshold of the small portion of the network.

Note: the depth parameters are not required; VBET can run based only on slope (the first version did), but results are generally significantly better by including depth (hence the update).

![VBET Output image](/pics/vbet_output.png)

![VBET Reach image](/pics/vbet_bitterroot.png)

![VBET Basin image](/pics/vbet_basin.png)

