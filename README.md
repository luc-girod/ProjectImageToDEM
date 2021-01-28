# ProjectImageToDEM
Luc Girod, Simon Filhol, January 2021

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
name x_img y_img x_world y_world z_world
p1 -557.5 251 188.439 108.181 9.922
p2 56.5 286 218.21 114.282 4.2
p3 -353.5 347 196.444 115.668 3.113
```
#### 3.1.2. Compute Resection

```python
from projimdem import resection as rs
cam_file = 'CuczaDemoData/CamCucza.json'
GCP_file = 'CuczaDemoData/GCPs_Centered.csv'
test = rs.resection(cam_file, GCP_file)
test.estimate_cam()
```

Resection uses `scipy.optimize.least_square` to estimate the new camera orientation.

### 3.2. Project image on DEM