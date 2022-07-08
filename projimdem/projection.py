# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 16:27:46 2021

@author: lucg
"""
import numpy as np
from osgeo import gdal
from osgeo import osr
from matplotlib import pyplot
import time
import os, pdb
import pandas as pd

# ignore pandas warnings
import warnings
warnings.filterwarnings("ignore")


class Projection():
    def __init__(self, dem_file, viewshed_file, image_file, cam_param, output_file, dem_nan_value=-9999):
        '''
        dem_file: geotif of the DEM
        viewshed_file: name for the viewshade file (will be computed by gdal based on camera parameters). If None, it will project ofr all points of the DEM within the theoriticial FoV
        image_file: image file to project on to DEM. JPG format, RGB
        cam_param: camera paramters expressed as in the output of Resection.
        output_file:
        dem_nan_value:
        
        TODO: 
        - review the viewshed file. add option to not have to compute it each time
        '''

        self.output_file = output_file
        self.cam_param = cam_param
        # Load in DEM
        self.dem_raster = gdal.Open(dem_file)
        self.geot = self.dem_raster.GetGeoTransform()
        self.Xsize = self.dem_raster.RasterXSize
        self.Ysize = self.dem_raster.RasterYSize
        self.dem_data = self.dem_raster.GetRasterBand(1).ReadAsArray(0, 0, self.Xsize, self.Ysize)
        self.dem_data[self.dem_data == dem_nan_value] = np.nan
        # define extent and resolution from geotiff metadata
        self.extent = [self.geot[0], self.geot[0] + np.round(self.geot[1], 3) * self.Xsize, self.geot[3],
                       self.geot[3] + np.round(self.geot[5], 3) * self.Ysize]
        self.Xs = np.linspace(self.extent[0] + np.round(self.geot[1], 3), self.extent[1], self.Xsize)
        self.Ys = np.linspace(self.extent[2] + np.round(self.geot[5], 3), self.extent[3], self.Ysize)
        self.viewshed_file = viewshed_file
        # Compute DEM Z at Camera XY
        col = int((self.cam_param[0][0] - self.geot[0]) / self.geot[1])
        row = int((self.cam_param[0][1] - self.geot[3]) / self.geot[5])
        ZDEMatCamera = self.dem_data[row][col]
        print("Z_DEM_at_Camera", ZDEMatCamera)
        print("Z_Camera", self.cam_param[0][2])
        ZCameraOverDEM = self.cam_param[0][2] - ZDEMatCamera + 20 # +20 so the viewshed is fuller rather than a bit underestimated

        # Compute viewshade file using GDAL and the new camear prosition
        if self.viewshed_file is not None:
            command = 'gdal_viewshed -ox ' + str(self.cam_param[0][0]) + ' -oy ' + str(
                self.cam_param[0][1]) + ' -oz ' + str(ZCameraOverDEM) + ' ' + str(dem_file) + ' ' + str(viewshed_file)
            print(command)
            os.system(command)
            self.viewshed_raster = gdal.Open(viewshed_file)
            self.viewshed_data = self.viewshed_raster.GetRasterBand(1).ReadAsArray(0, 0, self.Xsize, self.Ysize)
            self.viewshed_raster = None

        self.image = pyplot.imread(image_file)

        # Create output object
        self.ortho = np.zeros((self.Ysize, self.Xsize, 3), dtype=np.uint8)

        # Compute remove radial and tangential distortion based on camera parameters
        self.X_undistort, self.Y_undistort, self.image_undistort = self.img_correct_distortion()
        self.pt_proj = pd.DataFrame()

    def XYZ_to_img(self, pts_world):
        imsize = self.image.shape
        pts_cam = pts_world.copy()
        Foc = self.cam_param[2]
        [DCx, DCy] = self.cam_param[3]
		
        # World to camera coordinate
        XYZ_world = pts_cam[['X_world', 'Y_world', 'Z_world']].values
        pts_cam[['X_cam', 'Y_cam', 'Z_cam']] = pd.DataFrame(
            (self.cam_param[1].dot(np.subtract(XYZ_world, self.cam_param[0]).T)).T)

        XYZ_world = None

        # Removing points behind the camera
        pts_cam.Z_cam.loc[pts_cam.Z_cam > 0] = np.nan
        pts_cam = pts_cam.dropna(axis=0)
		
        # Projecting points to sensor (through a non distorting lens)
        pts_cam['X_proj'] = pts_cam.X_cam.values / pts_cam.Z_cam.values
        pts_cam['Y_proj'] = pts_cam.Y_cam.values / pts_cam.Z_cam.values
        pts_cam['X_img'] = - pts_cam['X_proj'] * Foc + DCx
        pts_cam['Y_img'] = - pts_cam['Y_proj'] * Foc + DCy
		
        # Removing points outside of the camera field of view
        pts_cam.X_img.loc[(pts_cam.X_img < 0) | (pts_cam.X_img > imsize[1])] = np.nan
        pts_cam.Y_img.loc[(pts_cam.Y_img < 0) | (pts_cam.Y_img > imsize[0])] = np.nan
        pts_cam = pts_cam.dropna(axis=0)
		
        		
        pts_cam['X_distort'] = self.X_undistort[pts_cam.Y_img.astype(int), pts_cam.X_img.astype(int)] - self.X_undistort.astype(int).min()
        pts_cam['Y_distort'] = self.Y_undistort[pts_cam.Y_img.astype(int), pts_cam.X_img.astype(int)] - self.Y_undistort.astype(int).min()

        pts_cam.X_distort.loc[(pts_cam.X_distort < 0) | (pts_cam.X_distort > self.image_undistort.shape[1])] = np.nan
        pts_cam.Y_distort.loc[(pts_cam.Y_distort < 0) | (pts_cam.Y_distort > self.image_undistort.shape[0])] = np.nan
        pts_cam = pts_cam.dropna(axis=0)

        return pts_cam

    def img_correct_distortion(self, return_image=True):
        """
        Function to correction radial and tangential distortion of an image
        :param img: original image
        :return 
        """
        img = self.image
        Foc = self.cam_param[2]
        [DCx, DCy] = self.cam_param[3]
        [K1, K2, K3, K4, K5, K6, P1, P2, P3, P4, P5, P6, P7] = self.cam_param[4]

        Xs_dis, Ys_dis = np.meshgrid(np.arange(0, img.shape[1]),
                                     np.arange(0, img.shape[0]))
        # Distortion model        
        X_centered=(Xs_dis - DCx) / Foc
        Y_centered=(Ys_dis - DCy) / Foc
        R = np.sqrt(pow(X_centered, 2) + pow(Y_centered, 2))

        x_im_nodist = Foc * X_centered * (
                1 + K1 * pow(R, 2) + K2 * pow(R, 4) + K3 * pow(R, 6)) / (
                1 + K4 * pow(R, 2) + K5 * pow(R, 4) + K6 * pow(R, 6)) + (
                              P1 * (pow(R, 2) + 2 * pow(X_centered, 2)) + 2 * P2 * X_centered) * Y_centered * (
                              1 + P3 * pow(R, 2) + P4 * pow(R, 4) + P5 * pow(R, 6) + P6 * pow(R, 8) + P7 * pow(R, 10))
                                  
        y_im_nodist = Foc * Y_centered * (
                1 + K1 * pow(R, 2) + K2 * pow(R, 4) + K3 * pow(R, 6)) / (
                1 + K4 * pow(R, 2) + K5 * pow(R, 4) + K6 * pow(R, 6)) + (
                              P1 * (pow(R, 2) + 2 * pow(Y_centered, 2)) + 2 * P2 * Y_centered) * X_centered * (
                              1 + P3 * pow(R, 2) + P4 * pow(R, 4) + P5 * pow(R, 6) + P6 * pow(R, 8) + P7 * pow(R, 10))

        if return_image:

            # Project distorted image
            img_cor = np.empty([y_im_nodist.astype(int).max() - y_im_nodist.astype(int).min() + 1,
                                x_im_nodist.astype(int).max() - x_im_nodist.astype(int).min() + 1, 3])
            img_cor[(y_im_nodist.astype(int) - y_im_nodist.astype(int).min()),
            (x_im_nodist.astype(int) - x_im_nodist.astype(int).min()), :] = img
            img_cor = img_cor.astype(int)

            img = None
            return x_im_nodist, y_im_nodist, img_cor
        else:
            return x_im_nodist, y_im_nodist

    def project_img_to_DEM(self, return_raster=True, epsg=None):
        '''
        Function to project an image to a DEM
        '''

        # Project all DEM points in the viewshed to the image
        tic = time.perf_counter()

        x_mesh, y_mesh = np.meshgrid(self.Xs, self.Ys)
        pt_world = pd.DataFrame()
        if self.viewshed_file is not None:
            pt_world['X_world'] = x_mesh[self.viewshed_data == 255].flatten()
            pt_world['Y_world'] = y_mesh[self.viewshed_data == 255].flatten()
            pt_world['Z_world'] = self.dem_data[self.viewshed_data == 255].flatten()
            pt_world['dist_cam'] = np.sqrt((pt_world.X_world - self.cam_param[0][0]) ** 2 +
                                           (pt_world.Y_world - self.cam_param[0][1]) ** 2 +
                                           (pt_world.Z_world - self.cam_param[0][2]) ** 2)
            pt_world['ind_xdem'] = np.argwhere(self.viewshed_data == 255)[:, 1]
            pt_world['ind_ydem'] = np.argwhere(self.viewshed_data == 255)[:, 0]
        else:
            pt_world['X_world'] = x_mesh.flatten()
            pt_world['Y_world'] = y_mesh.flatten()
            pt_world['Z_world'] = self.dem_data.flatten()
            pt_world['dist_cam'] = np.sqrt((pt_world.X_world - self.cam_param[0][0]) ** 2 +
                                           (pt_world.Y_world - self.cam_param[0][1]) ** 2 +
                                           (pt_world.Z_world - self.cam_param[0][2]) ** 2)
            pt_world['ind_xdem'] = np.argwhere(~np.isnan(self.dem_data))[:, 1]
            pt_world['ind_ydem'] = np.argwhere(~np.isnan(self.dem_data))[:, 0]

        self.pt_proj = self.XYZ_to_img(pt_world)

        # search for each DEM point in the viewshed its RGB value in the undistorted image to compute orthophoto
        for i, row in self.pt_proj.iterrows():
            self.ortho[row.ind_ydem.astype(int), row.ind_xdem.astype(int), :] = self.image_undistort[
                                                                                row.Y_distort.astype(int),
                                                                                row.X_distort.astype(int), :]

        toc = time.perf_counter()

        print(f"Ortho computed in {toc - tic:0.4f} seconds")

        if return_raster:
            fileformat = "GTiff"
            driver = gdal.GetDriverByName(fileformat)
            metadata = driver.GetMetadata()
            if metadata.get(gdal.DCAP_CREATE) == "YES":
                print("Driver {} supports Create() method.".format(fileformat))

            if metadata.get(gdal.DCAP_CREATECOPY) == "YES":
                print("Driver {} supports CreateCopy() method.".format(fileformat))

            dst_ds = driver.Create(self.output_file, xsize=self.Xsize, ysize=self.Ysize,
                                   bands=3, eType=gdal.GDT_Byte)

            dst_ds.SetGeoTransform(self.geot)
            srs = osr.SpatialReference()
            if epsg is None:
                srs.ImportFromWkt(self.dem_raster.GetProjectionRef())
            else:
                srs.ImportFromEPSG(epsg)
            dst_ds.SetProjection(srs.ExportToWkt())
            dst_ds.GetRasterBand(1).WriteArray(self.ortho[:, :, 0])
            dst_ds.GetRasterBand(2).WriteArray(self.ortho[:, :, 1])
            dst_ds.GetRasterBand(3).WriteArray(self.ortho[:, :, 2])
            # Once we're done, close properly the dataset
            dst_ds = None
            print('Ortho computed')

        return 0

    def plot_DEM_on_img(self, color_variable='Z_world', **kwargs):
        pyplot.figure()
        pyplot.imshow(self.image_undistort)
        pyplot.scatter(self.pt_proj.X_distort,
                       self.pt_proj.Y_distort,
                       c=self.pt_proj[color_variable],
                       **kwargs)
        pyplot.colorbar()


# =================================================================================
# --------------------------  Practical Functions -------------------------
# =================================================================================


def img_correct_distortion(img, cam_param):
    """
    Function to correction radial and tangential distortion of an image
    :param img: original image
    :param cam_param: Object containing camera calibration coefficients. extracted from the resection.
    :return 
    """

    Xs_dis, Ys_dis = np.meshgrid(np.arange(-img.shape[1] / 2, img.shape[1] / 2),
                                 np.arange(-img.shape[0] / 2, img.shape[0] / 2))
    
    Foc=cam_param[2]
    [DCx, DCy] = cam_param[3]
    [K1, K2, K3, K4, K5, K6, P1, P2, P3, P4, P5, P6, P7] = cam_param[4]
    
    X_centered=(Xs_dis - DCx) / Foc
    Y_centered=(Ys_dis - DCy) / Foc
    R = np.sqrt(pow(X_centered, 2) + pow(Y_centered, 2))

    x_im_nodist = Foc * X_centered * (
            1 + K1 * pow(R, 2) + K2 * pow(R, 4) + K3 * pow(R, 6)) / (
            1 + K4 * pow(R, 2) + K5 * pow(R, 4) + K6 * pow(R, 6)) + (
                          P1 * (pow(R, 2) + 2 * pow(X_centered, 2)) + 2 * P2 * X_centered) * Y_centered * (
                          1 + P3 * pow(R, 2) + P4 * pow(R, 4) + P5 * pow(R, 6) + P6 * pow(R, 8) + P7 * pow(R, 10))
                              
    y_im_nodist = Foc * Y_centered * (
            1 + K1 * pow(R, 2) + K2 * pow(R, 4) + K3 * pow(R, 6)) / (
            1 + K4 * pow(R, 2) + K5 * pow(R, 4) + K6 * pow(R, 6)) + (
                          P1 * (pow(R, 2) + 2 * pow(Y_centered, 2)) + 2 * P2 * Y_centered) * X_centered * (
                          1 + P3 * pow(R, 2) + P4 * pow(R, 4) + P5 * pow(R, 6) + P6 * pow(R, 8) + P7 * pow(R, 10))

    img_cor = np.empty([y_im_nodist.astype(int).max() - y_im_nodist.astype(int).min() + 1,
                        x_im_nodist.astype(int).max() - x_im_nodist.astype(int).min() + 1, 3])
    img_cor[(y_im_nodist.astype(int) - y_im_nodist.astype(int).min()),
    (x_im_nodist.astype(int) - x_im_nodist.astype(int).min()), :] = img
    img_cor = img_cor.astype(int)

    return x_im_nodist, y_im_nodist, img_cor
