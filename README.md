# VBET 2.0
Valley Bottom Extraction Tool (upated: 11/25/2020)

VBET uses a stream network shapefile and digital elevation model to derive a valley bottom polygon based on two lines of evidence: slope and inundation depth. The values for these parameters are adjusted based on the drainage area of a given segment of the network based on the generalizations that at lower drainage areas (small streams), slopes within valley bottoms are generally steeper than at high drainage areas (larger rivers), where slopes tend to be very flat in valley bottoms; and also that high in the network, at low drainage areas, flood depths tend to be shallower than in larger portions of the network.

As such, the user selects to drainage area threshold values that are used to split the network into 'large', 'medium' and 'small' portions. The user then selects threshold values for slope and inundation depth for each of the three portions of the network.

## Required Python packages
### Python 3
- geopandas (v. 0.8.0)
- rasterio (v. 1.1.1)
- shapely
- rasterstats
- numpy 
- scipy
- scikit-image (v. 0.16.2)

## Data preparation
Because the tool relies on drainage area, a drainage area raster is required to run the tool. To be accurate, this raster must encompass the entire watershed upstream of the area for which a valley bottom is being delineated. In cases where a small LiDAR dataset is being used to delineate valley bottoms for a portion of a drainage network, a coarser DEM of the entire basin (e.g. the 10m DEMs of the NED from USGS) should be used to generate this raster. This can be accomplished using the common GIS workflow of filling the pits in the DEM, generating flow directions, and then a flow accumulation raster. The flow accumulation raster can be converted to drainage area in a raster calculator using the equation (flow_accumulation * (flow_accumulation_raster_resolution^2))/1000000.

For the best results, the drainage network should align well with the DEM used. If the network was generated using the DEM, this isn't a problem, but when using external datasets such as NHD, there can be a large spatial disparity between the location of the drainage network in the DEM and in the actual dataset. If this is the case, to the extent possible, the user should edit the network, adjusting it so that it aligns with the network location in the DEM. Additionally, the segment length of the network should be suitable given the DEM resolution. That is, with higher resolution DEMs the network should be segmented at smaller lengths than for coarser DEMs (e.g. 10m). A good length is generally 100-500 times the DEM resolution (e.g. for a 10m DEM, we use 1000m segments, and for a 1m DEM we use 200-500m segments, depending on river size). Stream networks must contain ONLY SINGLEPART FEATURES. If it contains multipart features, use a "multipart to singlepart" GIS tool to fix this issue.

## Running the model
To run, download the repository, and open the 'run_VBET.py' script in a Python IDE.
Fill out the parameters as described below, and run the script.

## Parameters 

- network: enter the path to the stream network shapefile. This can be e.g. the National Hydrologic Dataset (NHD), extracted from a DEM, or manually delineated by the user. The stream network should be segmented so that each reach is approximately the desired segment length.
- dem: enter the path to the DEM from which you wish to derive the valley bottom polygon (.tif). 
- out: enter a path and file name to store the valley bottom shapefile produced by the tool.
- scratch: enter a path to a scratch workspace folder where temporary files will be stored. If the entered folder does not exist it will be created.
- lg_da: enter a value (int) for a drainage area threshold above which the network will be considered 'large'.
- med_da: enter a value (int) for a drainage area threshold above which the network will be considered 'medium', and below which the network will be considered 'small'.
- lg_slope: enter a value for the slope threshold in the large portion of the network. 
- med_slope: enter a value for the slope threshold in the medium portion of the network.
- sm_slope: enter a value for the slope threshold in the small portion of the network.
- lg_buf: enter a value representing the maximum distance from the large portion of the network within which valley bottoms occur.
- med_buf: enter a value representing the maximum distance from the medium portion of the network within which valley bottoms occur.
- sm_buf: enter a value representing the maximum distance from the small portion of the network within which valley bottoms occur.
- dr_area: enter the path to the drainage area raster for the basin (.tif).
- lg_depth: enter a value for the depth threshold of the large portion of the network.
- med_depth: enter a value for the depth threshold of the medium portion of the network.
- sm_depth: enter a value for the depth threshold of the small portion of the network.

Note: the depth parameters are not required; VBET can run based only on slope (the first version did), but results are generally significantly better by including depth (hence the update).

![VBET Output image](/pics/vbet_output.png)

![VBET Reach image](/pics/vbet_bitterroot.png)

![VBET Basin image](/pics/vbet_basin.png)

