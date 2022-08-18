# ProjectImageToDEM
Luc Girod, Simon Filhol, January 2021

## 0. TODO
**Short Term:**
- [x] Implement distortion correction in projection to DEM
- [ ] test implementing weight on GCPs based on approximative distance from camera
- [ ] we now have a two step least square. The ideal case would do all at once, but we need to handle the spread in order of magnitude of the different unknowns (in particular the distortion coefficients vs. the rest)
- [ ] add option in projection to either correct distortion or assume the image is calibrated (current version)
- [x] go through class and function names for better clarity
- [x] Clean up script RunScript.py. 
- [x] update README examples after changes

**Long term:**
- [ ] add quality assessment based on spliting GCPs in two
- [x] add plotting functionalities for various steps of the processing
- [ ] add method with OpenCV to extract camera parameter from Chessboard, See if can use calibration from Micmac
- [x] make a "perfect" test with a camera of known calibration distortion parameters and location. 
- [x] add method to project DEM in the image (all variable avail., simply needs proper plotting function)
- [ ] implement function to classify image snow/no snow based on thresholding
- [ ] add correction in elevation due to Earth curvature [Corripio (2004)](https://www.tandfonline.com/doi/pdf/10.1080/01431160410001709002?needAccess=true) page 12
- [ ] implement snow/no-snow classification functions in `classification.py`
- [x] add interpolation of the undistorted image for pixel with no value

## 1. Description

**projimdem** is Python package to 

1. compute camera resection from a list of Ground Control Points (GCP) expressed in image coordinates (x,y) with their respective world coordinates (x,y,z). 
2. Project the image on a DEM.

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
    "omega":1.49,
    "phi":-0.39,
    "kappa":-0.39,
    "X_ini":419169.2,
    "Y_ini":6718421.3,
    "Z_ini":1212},
"iop":{
    "Foc":1484,
    "DCx":0,
    "DCy":0,
    "R1":0,
    "R3":0,
    "R5":0,
    "K1":0,
    "K2":0,
    "K3":0,
    "K4":0,
    "K5":0,
    "P1":0,
    "P2":0,
    "P3":0,
    "P4":0,
    "P5":0,
    "P6":0,
    "P7":0
    }   
}
```
`eop` for external orientation parameters. `phi, omega, kappa` are the orientatino angle in radians:
- omega is the angle from looking straight up (pi/2=1.57 is the camera being level, the value should be in [0:pi])
- phi the angle from the camera being level (pi/2=1.57 is the camera turned to portrait mode with the right side of the image pointing up)
- kappa is the angle from looking straight north (pi/2=1.57 is pointing west, pi=3.14 is pointing south, -pi/2=-1.57 is pointing East)

`iop` are the internal orientation parameter with `Foc` the focal length in pixels, `DCx` and `DCy` the pricipal point (distortion centers), and R1-5, K1-5 and P1-7 the supported distortion parameters. 

2. save the GCPs image and world coordinates in a csv file as follow:
```csv. The image coordinate system has its origin in the center of the image. 


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
