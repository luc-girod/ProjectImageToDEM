# ProjectImageToDEM
Luc Girod, Simon Filhol, January 2021

## 0. TODO
- [x] Implement distortion correction in projection to DEM
- [ ] test implementing weight on GCPs based on distance from camera
- [ ] we now have a two step least square. The ideal case woulddo all at once, but we need to handle the spread in order of magnitude of the different unknowns (in particular the distortion coefficients vs. the rest)
- [ ] Clean up script RunScript.py. 
- [ ] add option in projection to either correct distortion or assume the image is calibrated (current version)
- [ ] go through class and function names for better clarity
- [ ] update README examples after changes

## 1. Description

**projimdem** is Python package to 

1. compute camera resection from a list of Ground Control Points (GCP) expressed in image (x,y) with their respective world coordinates (x,y,z). 
2. Project the image on a DEM.

## 

## 2. Installation
```sh
# pull the repo locally from github
git pull git@github.com:luc-girod/ProjectImageToDEM.git

# install in developement mode
pip install -e .
```
**Required package:** gdal (v3.1), numpy, pandas, maptlotlib, scipy, json, osgeo
## 3. Usage
### 3.1. Resection from GCPs

#### 3.1.1. Inputs

Input needed for the resection are: 
1. create an initialization JSON file for the camera 

```json
{
"eop":{
    "omega":0,
    "phi":0,
    "kappa":0,
    "X_ini":209.89527614679403023,
    "Y_ini":91.205,
    "Z_ini":107.031846453497209},
"iop":{
    "x0":0,
    "y0":0,
    "Foc":2011.887
    }   
}
```
`eop` for external orientation parameters. `phi, omega, kappa` are the orientatino angle in radians. `iop` are the internal orientation parameter with `Foc` the focal length in pixels, `x0` and `y0` 

2. save the GCPs image and world coordinates in a csv file as follow:
```csv
name x_img y_img x_world y_world z_world lstsq_IO
p1 -557.5 251 188.439 108.181 9.922 1
p2 56.5 286 218.21 114.282 4.2 0
p3 -353.5 347 196.444 115.668 3.113 1
```
#### 3.1.2. Compute Resection

```python
from projimdem import resection as rs
from projimdem import ProjectImageOnDEM as pi
import numpy as np

cam_file = './example/FinseFromPhoto4D/CamFinseInit.json'
GCP_file = './example/FinseFromPhoto4D/GCPs.csv'
dem_file = './example/FinseFromPhoto4D/DEM_2m.tif'
viewshed_file = './example/FinseFromPhoto4D/viewshed_test.tif'
image_file = './example/FinseFromPhoto4D/DSC03111.JPG'
output_file = './example/FinseFromPhoto4D/finse_proj_2m.tif'

#================================
# 1. Compute transformation
finse = rs.resection(cam_file, GCP_file, image_file, delimiter_GCP=' ',
                    free_param=['omega', 'phi', 'kappa','X_ini','Y_ini','Z_ini', 'Foc'],
                    param_bounds=([-3.15, -3.15, -3.15,-np.inf,-np.inf,-np.inf,4000], [3.15,3.15,3.15,np.inf,np.inf,np.inf,6000]))
#finse.GCPs
finse.estimate_cam()
finse.print_residuals()

#================================
# 2. Visualize transform
finse.proj_GCPs2img()
```

Resection uses `scipy.optimize.least_square` to estimate the new camera orientation. Some optoins can be passed to `estimate_cam()` (see source code).

### 3.2. Project image on DEM

```python
#================================
# 3. Projection
finse_proj = pi.ProjIm2dem(dem_file = dem_file,
                          viewshed_file = viewshed_file,
                          image_file = image_file,
                          cam_param = finse.new_cam.proj_param,
                           output_file = output_file)
finse_proj.ProjectImage2DEM(return_raster = True, epsg=32632)
```
