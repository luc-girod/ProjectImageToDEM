#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Adaptation by S. Filhol from:

> Author:  Jeffrey T. Walton, Paul Smith's College, New York
>
>   Single-photo resection - calculates the camera orientation and location
>       given camera calibration parameters, control point photo and world
>        coordinates and initial guesses for camera exterior orientation.
>
>   based on MATLAB code from:
>   Introduction to Modern Photogrammetry by Mikhail, Bethel, McGlone
>   John Wiley & Sons, Inc. 2001
>
> https://github.com/jeffwalton/photogrammetry-resection



'''


import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize 
import json
import pandas as pd
from types import SimpleNamespace
from math import sin, cos

# ignore pandas warnings
import warnings
warnings.filterwarnings("ignore")

class Resection():
    '''
    Class to compute resection parameter from camera parameters and GCPs 
    '''
    
    def __init__(self, camera_file, GCP_file, image_file, delimiter_GCP=' ', x_offset=None, y_offset=None, z_offset=None, free_param=['omega', 'phi', 'kappa'], param_bounds=([-3.15, -3.15, -3.15], [3.15,3.15,3.15])):
        '''
        camera_file: json file with camera parameters
        GCP_file: csv file with GCPs names and coordinates
        image_file: image file (jpeg)
        delimiter_GCP: delimiter of the csv file GCP
        x_offset: offset in X to center the least square
        y_offset: offset in Y to center the least square
        z_offset: offset in Z to center the least square
        free_param: list of choice of free parameters to fit the least square. Can be 'omega', 'phi', 'kappa', 'X_ini', 'Y_ini', 'Z_ini', 'Foc' and radial distortion parameters 'DCx', 'DCy', 'K1', 'K2', 'K3', 'K4' and 'K5', and tangential parameters 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7'
        param_bounds: free parameter min and max value to bound the least square. Must be ([min1, min2, min3, ...],[max1, max2, max3, ...])
        
        '''
        
        #Load camera parameters
        with open(camera_file, 'r') as myfile:
            self.cam = json.loads(myfile.read(), object_hook=lambda d: SimpleNamespace(**d))
        
        # Load GCP coordinates
        # headers must be: name x_img y_img x_world y_world z_world
        self.GCPs = pd.read_csv(GCP_file, delimiter=delimiter_GCP)
        
        # define offesets if not provided
        if x_offset is None:
            self.x_offset = self.cam.eop.X_ini
        else:
            self.x_offset = x_offset
            
        if y_offset is None:
            self.y_offset = self.cam.eop.Y_ini
        else:
            self.y_offset = y_offset
        
        if z_offset is None:
            self.z_offset = self.cam.eop.Z_ini
        else:
            self.z_offset = z_offset
        
        self.free_param = free_param
        
        # Build the x0 vector including the tunning parameter for the least square
        self.x0_dict = {}
        for param in free_param:
            if param in ['omega', 'phi', 'kappa', 'X_ini', 'Y_ini', 'Z_ini']:
                p = self.cam.eop.__getattribute__(param)
                self.x0_dict[param] = p
                if param =='X_ini':
                    self.x0_dict[param] = 0#p - self.x_offset
                    # self.x_offset=self.x0_dict[param]
                elif param =='Y_ini':
                    self.x0_dict[param] = 0#p - self.y_offset
                    # self.y_offset=self.x0_dict[param]
                elif param =='Z_ini':
                    self.x0_dict[param] = 0#p - self.z_offset
                    # self.z_offset=self.x0_dict[param]
            #if param in ['Foc', 'DCx', 'DCy','R1', 'R3', 'R5']:
            if param in ['Foc', 'DCx', 'DCy','K1', 'K2', 'K3', 'K4', 'K5', 'K6', 'P1','P2','P3','P4', 'P5', 'P6', 'P7']:
                p = self.cam.iop.__getattribute__(param)
                self.x0_dict[param] = p
        self.x0 = list(self.x0_dict.values())
        self.param_bounds = (param_bounds)
        print(self.x0)
        # class to store camera parameters after least square
        class new_cam:
            def __init__(self):
                self.center = None
                self.rotation = None
                self.proj_param = None
                self.omega = None
                self.kappa = None
                self.phi = None
        self.new_cam = new_cam
        
        
        self.GCPs['x_world_offset'] = self.GCPs.x_world - self.x_offset
        self.GCPs['y_world_offset'] = self.GCPs.y_world - self.y_offset
        self.GCPs['z_world_offset'] = self.GCPs.z_world - self.z_offset
        
        # initialize least square with direct projection
        res_ini = self.collinearity_func(self.x0)
        
        idx = np.arange(0,res_ini.__len__(),2)
        self.GCPs['residual_x_lstsq'] = np.nan
        self.GCPs['residual_y_lstsq'] = np.nan
        self.GCPs['residual_x_ini'] = np.nan
        self.GCPs['residual_y_ini'] = np.nan
        self.GCPs['residual_x_ini'].loc[self.GCPs.lstsq_IO.astype(bool)] = res_ini[idx]
        self.GCPs['residual_y_ini'].loc[self.GCPs.lstsq_IO.astype(bool)] = res_ini[idx+1]
        
        # load image
        self.image = plt.imread(image_file)

    
    def change_free_params(self, free_param=['omega', 'phi', 'kappa'], param_bounds=([-3.15, -3.15, -3.15], [3.15,3.15,3.15])):
        self.cam.eop.omega = self.new_cam.omega
        self.cam.eop.kappa = self.new_cam.kappa
        self.cam.eop.phi = self.new_cam.phi
        self.cam.iop.Foc = self.new_cam.Foc
        self.cam.eop.X_ini = self.new_cam.center[0]
        self.cam.eop.Y_ini = self.new_cam.center[1]
        self.cam.eop.Z_ini = self.new_cam.center[2]
        self.cam.iop.K1 = self.new_cam.K1
        self.cam.iop.K2 = self.new_cam.K2
        self.cam.iop.K3 = self.new_cam.K3
        self.cam.iop.K4 = self.new_cam.K4
        self.cam.iop.K5 = self.new_cam.K5
        self.cam.iop.K6 = self.new_cam.K6
        self.cam.iop.P1 = self.new_cam.P1
        self.cam.iop.P2 = self.new_cam.P2
        self.cam.iop.P3 = self.new_cam.P3
        self.cam.iop.P4 = self.new_cam.P4
        self.cam.iop.P5 = self.new_cam.P5
        self.cam.iop.P6 = self.new_cam.P6
        self.cam.iop.P7 = self.new_cam.P7
        self.cam.iop.DCx = self.new_cam.DCx
        self.cam.iop.DCy = self.new_cam.DCy

        self.x0_dict = {}
        for param in free_param:
            if param == 'omega':
                self.x0_dict['omega'] = self.new_cam.omega
            elif param == 'kappa':
                self.x0_dict['kappa'] = self.new_cam.kappa
            elif param == 'phi':
                self.x0_dict['phi'] = self.new_cam.phi
            elif param == 'Foc':
                self.x0_dict['Foc'] = self.new_cam.Foc
            elif param =='X_ini':
                self.x0_dict['X_ini'] = self.new_cam.center[0] - self.x_offset
            elif param =='Y_ini':  
                self.x0_dict['Y_ini'] = self.new_cam.center[1] - self.y_offset
            elif param =='Z_ini':
                self.x0_dict['Z_ini'] = self.new_cam.center[2] - self.z_offset
            elif param == 'K1':
                self.x0_dict['K1'] = self.new_cam.K1
            elif param == 'K2':
                self.x0_dict['K2'] = self.new_cam.K2
            elif param == 'K3':
                self.x0_dict['K3'] = self.new_cam.K3
            elif param == 'K4':
                self.x0_dict['K4'] = self.new_cam.K4
            elif param == 'K5':
                self.x0_dict['K5'] = self.new_cam.K5
            elif param == 'P1':
                self.x0_dict['P1'] = self.new_cam.P1
            elif param == 'P2':
                self.x0_dict['P2'] = self.new_cam.P2
            elif param == 'P3':
                self.x0_dict['P3'] = self.new_cam.P3
            elif param == 'P4':
                self.x0_dict['P4'] = self.new_cam.P4
            elif param == 'P5':
                self.x0_dict['P5'] = self.new_cam.P5
            elif param == 'P6':
                self.x0_dict['P6'] = self.new_cam.P6
            elif param == 'P7':
                self.x0_dict['P7'] = self.new_cam.P7
            elif param == 'DCx':
                self.x0_dict['DCx'] = self.new_cam.DCx
            elif param == 'DCy':
                self.x0_dict['DCy'] = self.new_cam.DCy
                    
        for param in free_param:
            self.x0 = [self.x0_dict.get(x) for x in free_param]
            self.param_bounds = (param_bounds)
        print(self.x0)
        
        
        
        
    def rot_matrix_from_angles(self, omega, phi, kappa):
        '''
        Rotation matrix from angles following Micmac convention
        '''
        RX = np.array([[1,0,0],
                 [0, cos(omega), -sin(omega)],
                 [0, sin(omega), cos(omega)]])    
        RY = np.array([[cos(phi), 0, sin(phi)],
                 [0,1,0],    
                 [-sin(phi), 0, cos(phi)]])
        RZ = np.array([[cos(kappa),-sin(kappa),0],
                 [sin(kappa), cos(kappa),0],
                 [0,0,1]])
        M = RX.dot(RY.dot(RZ)).dot(np.array([[1,0,0],[0,-1,0],[0,0,1]]))
        return M
    
    def collinearity_func(self, indep_vars):
        """
        collinearity function calculates a sum of the squared residuals of the
            collinearity equations for all of the control points
        This function is passed to scipy.optimize.minimize()

        Inputs:
            indep_vars (passed) are the exterior orientation parameters of the camera

        Returns:
            sum of squared residuals of collinearity eqns
        """
        
        omega = self.cam.eop.omega
        phi = self.cam.eop.phi
        kappa = self.cam.eop.kappa
        X_ini = self.cam.eop.X_ini - self.x_offset
        Y_ini = self.cam.eop.Y_ini - self.y_offset
        Z_ini = self.cam.eop.Z_ini - self.z_offset
        Foc = self.cam.iop.Foc
        DCx = self.cam.iop.DCx
        DCy = self.cam.iop.DCy
        K1 = self.cam.iop.K1
        K2 = self.cam.iop.K2
        K3 = self.cam.iop.K3
        K4 = self.cam.iop.K4
        K5 = self.cam.iop.K5
        K6 = self.cam.iop.K6
        P1 = self.cam.iop.P1
        P2 = self.cam.iop.P2
        P3 = self.cam.iop.P3
        P4 = self.cam.iop.P4
        P5 = self.cam.iop.P5
        P6 = self.cam.iop.P6
        P7 = self.cam.iop.P7
        
        # logic to grab value from x0 no matter the order indicated
        if 'omega' in self.x0_dict.keys():
            omega = indep_vars[list(self.x0_dict.keys()).index('omega')]
        if 'phi' in self.x0_dict.keys():
            phi = indep_vars[list(self.x0_dict.keys()).index('phi')]
        if 'kappa' in self.x0_dict.keys():
            kappa = indep_vars[list(self.x0_dict.keys()).index('kappa')]
        if 'X_ini' in self.x0_dict.keys():
            X_ini = indep_vars[list(self.x0_dict.keys()).index('X_ini')]
        if 'Y_ini' in self.x0_dict.keys():
            Y_ini = indep_vars[list(self.x0_dict.keys()).index('Y_ini')]
        if 'Z_ini' in self.x0_dict.keys():
            Z_ini = indep_vars[list(self.x0_dict.keys()).index('Z_ini')]
        if 'Foc' in self.x0_dict.keys():
            Foc = indep_vars[list(self.x0_dict.keys()).index('Foc')]
        if 'DCx' in self.x0_dict.keys():
            DCx = indep_vars[list(self.x0_dict.keys()).index('DCx')]
        if 'DCy' in self.x0_dict.keys():
            DCy = indep_vars[list(self.x0_dict.keys()).index('DCy')]
        if 'K1' in self.x0_dict.keys():
            K1 = indep_vars[list(self.x0_dict.keys()).index('K1')]
        if 'K2' in self.x0_dict.keys():
            K2 = indep_vars[list(self.x0_dict.keys()).index('K2')]
        if 'K3' in self.x0_dict.keys():
            K3 = indep_vars[list(self.x0_dict.keys()).index('K3')]
        if 'K4' in self.x0_dict.keys():
            K4 = indep_vars[list(self.x0_dict.keys()).index('K4')]
        if 'K5' in self.x0_dict.keys():
            K5 = indep_vars[list(self.x0_dict.keys()).index('K5')]
        if 'K6' in self.x0_dict.keys():
            K6 = indep_vars[list(self.x0_dict.keys()).index('K6')]
        if 'P1' in self.x0_dict.keys():
            P1 = indep_vars[list(self.x0_dict.keys()).index('P1')]
        if 'P2' in self.x0_dict.keys():
            P2 = indep_vars[list(self.x0_dict.keys()).index('P2')]
        if 'P3' in self.x0_dict.keys():
            P3 = indep_vars[list(self.x0_dict.keys()).index('P3')]
        if 'P4' in self.x0_dict.keys():
            P4 = indep_vars[list(self.x0_dict.keys()).index('P4')]
        if 'P5' in self.x0_dict.keys():
            P5 = indep_vars[list(self.x0_dict.keys()).index('P5')]
        if 'P6' in self.x0_dict.keys():
            P6 = indep_vars[list(self.x0_dict.keys()).index('P6')]
        if 'P7' in self.x0_dict.keys():
            P7 = indep_vars[list(self.x0_dict.keys()).index('P7')]

        
        M = self.rot_matrix_from_angles(omega, phi, kappa)
        tmp = self.GCPs.loc[self.GCPs.lstsq_IO.astype(bool)].reset_index(drop=True)
        
        F = np.zeros(2*tmp.shape[0])
        
        
        for i, row in tmp.iterrows():
            uvw = M * np.matrix([[row.x_world_offset - X_ini], [row.y_world_offset - Y_ini], [row.z_world_offset - Z_ini]])
            xproj_nodist = -Foc * uvw[0,0] / uvw[2,0]
            yproj_nodist = -Foc * uvw[1,0] / uvw[2,0]
            
            # Brown-Conrady distortion model
            X_centered=(row.x_img - DCx) / Foc
            Y_centered=(row.y_img - DCy) / Foc
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
            
            resx = x_im_nodist - xproj_nodist
            resy = y_im_nodist - yproj_nodist
            F[2*i], F[2*i+1] = resx, resy
            
        return F
    
    def project_GCPs_to_img(self, plot=True):
        self.GCPs['x_img_repoj']=self.GCPs['x_img']-self.GCPs['residual_x_lstsq']
        self.GCPs['y_img_repoj']=self.GCPs['y_img']-self.GCPs['residual_y_lstsq']
        if plot:
            fig, ax = plt.subplots(1,1)
            ax.imshow(self.image)
            ax.scatter(self.GCPs['x_img'],(self.GCPs['y_img']), label='Original positions')
            for i, txt in enumerate(self.GCPs['name']):
                ax.annotate(self.GCPs['name'][i], (self.GCPs['x_img'][i],(self.GCPs['y_img'][i])),color='blue', fontsize=8)


            ax.scatter(self.GCPs['x_img_repoj'],self.GCPs['y_img_repoj'], label='Reprojected positions')
            for i, txt in enumerate(self.GCPs['name']):
                ax.annotate(self.GCPs['name'][i], (self.GCPs['x_img_repoj'][i],self.GCPs['y_img_repoj'][i]),color='orange', fontsize=8)
            ax.legend()                   
        
    
    def estimate_cam(self, method='dogbox', loss='cauchy', verbose=2, f_scale=1, xtol=1e-8, ftol=1e-08, gtol=1e-09): 
        '''
        Function to perform camera estimation of the free parameters using least square.
        
        see: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html
        '''
        
        # perform least square
        res = optimize.least_squares(self.collinearity_func, self.x0, 
                                     loss=loss, method=method, verbose=verbose,
                                    bounds=self.param_bounds, f_scale=f_scale, jac='3-point', ftol=ftol, xtol=xtol, gtol=gtol)
        self.new_cam.RMSE = np.sqrt(np.sum(res.fun**2)/res.fun.__len__())
        self.new_cam.lstsq_results = res
        
        # extract new camera info from the least square results
        if 'omega' in self.x0_dict.keys():
            self.new_cam.omega = res.x[list(self.x0_dict.keys()).index('omega')]
        else:
            self.new_cam.omega = self.cam.eop.omega
        if 'phi' in self.x0_dict.keys():
            self.new_cam.phi = res.x[list(self.x0_dict.keys()).index('phi')]
        else:
            self.new_cam.phi = self.cam.eop.phi
        if 'kappa' in self.x0_dict.keys():
            self.new_cam.kappa = res.x[list(self.x0_dict.keys()).index('kappa')]
        else:
            self.new_cam.kappa = self.cam.eop.kappa
        if 'X_ini' in self.x0_dict.keys():
            self.new_cam.X_ini = res.x[list(self.x0_dict.keys()).index('X_ini')]
        else:
            self.new_cam.X_ini = self.cam.eop.X_ini - self.x_offset
        if 'Z_ini' in self.x0_dict.keys():
            self.new_cam.Z_ini = res.x[list(self.x0_dict.keys()).index('Z_ini')]
        else:
            self.new_cam.Z_ini = self.cam.eop.Z_ini - self.z_offset
        if 'Y_ini' in self.x0_dict.keys():
            self.new_cam.Y_ini = res.x[list(self.x0_dict.keys()).index('Y_ini')]
        else:
            self.new_cam.Y_ini = self.cam.eop.Y_ini - self.y_offset
        if 'Foc' in self.x0_dict.keys():
            self.new_cam.Foc = res.x[list(self.x0_dict.keys()).index('Foc')]
        else:
            self.new_cam.Foc = self.cam.iop.Foc     
        if 'DCx' in self.x0_dict.keys():
            self.new_cam.DCx = res.x[list(self.x0_dict.keys()).index('DCx')]
        else:
            self.new_cam.DCx = self.cam.iop.DCx
        if 'DCy' in self.x0_dict.keys():
            self.new_cam.DCy = res.x[list(self.x0_dict.keys()).index('DCy')]
        else:
            self.new_cam.DCy = self.cam.iop.DCy    
        if 'K1' in self.x0_dict.keys():
            self.new_cam.K1 = res.x[list(self.x0_dict.keys()).index('K1')]
        else:
            self.new_cam.K1 = self.cam.iop.K1
        if 'K2' in self.x0_dict.keys():
            self.new_cam.K2 = res.x[list(self.x0_dict.keys()).index('K2')]
        else:
            self.new_cam.K2 = self.cam.iop.K2
        if 'K3' in self.x0_dict.keys():
            self.new_cam.K3 = res.x[list(self.x0_dict.keys()).index('K3')]
        else:
            self.new_cam.K3 = self.cam.iop.K3
        if 'K4' in self.x0_dict.keys():
            self.new_cam.K4 = res.x[list(self.x0_dict.keys()).index('K4')]
        else:
            self.new_cam.K4 = self.cam.iop.K4
        if 'K5' in self.x0_dict.keys():
            self.new_cam.K5 = res.x[list(self.x0_dict.keys()).index('K5')]
        else:
            self.new_cam.K5 = self.cam.iop.K5
        if 'K6' in self.x0_dict.keys():
            self.new_cam.K6 = res.x[list(self.x0_dict.keys()).index('K6')]
        else:
            self.new_cam.K6 = self.cam.iop.K6
        if 'P1' in self.x0_dict.keys():
            self.new_cam.P1 = res.x[list(self.x0_dict.keys()).index('P1')]
        else:
            self.new_cam.P1 = self.cam.iop.P1
        if 'P2' in self.x0_dict.keys():
            self.new_cam.P2 = res.x[list(self.x0_dict.keys()).index('P2')]
        else:
            self.new_cam.P2 = self.cam.iop.P2
        if 'P3' in self.x0_dict.keys():
            self.new_cam.P3 = res.x[list(self.x0_dict.keys()).index('P3')]
        else:
            self.new_cam.P3 = self.cam.iop.P3
        if 'P4' in self.x0_dict.keys():
            self.new_cam.P4 = res.x[list(self.x0_dict.keys()).index('P4')]
        else:
            self.new_cam.P4 = self.cam.iop.P4
        if 'P5' in self.x0_dict.keys():
            self.new_cam.P5 = res.x[list(self.x0_dict.keys()).index('P5')]
        else:
            self.new_cam.P5 = self.cam.iop.P5
        if 'P6' in self.x0_dict.keys():
            self.new_cam.P6 = res.x[list(self.x0_dict.keys()).index('P6')]
        else:
            self.new_cam.P6 = self.cam.iop.P6
        if 'P7' in self.x0_dict.keys():
            self.new_cam.P7 = res.x[list(self.x0_dict.keys()).index('P7')]
        else:
            self.new_cam.P7 = self.cam.iop.P7

            
        self.new_cam.center = [self.new_cam.X_ini + self.x_offset, self.new_cam.Y_ini + self.y_offset, self.new_cam.Z_ini + self.z_offset]
        self.new_cam.rotation = self.rot_matrix_from_angles(self.new_cam.omega, self.new_cam.phi, self.new_cam.kappa)
        self.new_cam.distortion_center = [self.new_cam.DCx, self.new_cam.DCy]
        self.new_cam.distortion_params = [self.new_cam.K1, self.new_cam.K2, self.new_cam.K3, self.new_cam.K4, self.new_cam.K5, self.new_cam.K6, self.new_cam.P1, self.new_cam.P2, self.new_cam.P3, self.new_cam.P4, self.new_cam.P5, self.new_cam.P6, self.new_cam.P7]
        self.new_cam.proj_param = [self.new_cam.center, self.new_cam.rotation, self.new_cam.Foc, self.new_cam.distortion_center, self.new_cam.distortion_params]
        
        idx = np.arange(0,res.fun.__len__(),2)
        self.GCPs['residual_x_lstsq'] = np.nan
        self.GCPs['residual_y_lstsq'] = np.nan
        self.GCPs['residual_x_lstsq'].loc[self.GCPs.lstsq_IO.astype(bool)] = res.fun[idx]
        self.GCPs['residual_y_lstsq'].loc[self.GCPs.lstsq_IO.astype(bool)] = res.fun[idx+1]
        
        # Compute RMSE of the final fit
        RMSE_ini = np.sqrt(np.sum(self.GCPs.residual_x_ini**2 + self.GCPs.residual_y_ini**2)/res.fun.__len__())
        print('RMSE initial = ', RMSE_ini)
        print('RMSE lstsq = ', self.new_cam.RMSE)
        print('Camera = ', self.new_cam.proj_param)
        return self.new_cam.proj_param
    
    
    def print_residuals(self):
        '''
        Function to print to screen the least square residual before and after
        '''
        
        print(self.GCPs[['name','residual_x_ini', 'residual_y_ini', 'residual_x_lstsq',
       'residual_y_lstsq']].to_string())

    def plot_residuals(self):
        '''
        Visualize residuals
        '''
        fig, ax = plt.subplots(3,2,sharex=True, sharey=True)
        sc1 = ax[0,0].scatter(self.GCPs.x_img, self.GCPs.y_img, 
                              c=self.GCPs.residual_x_ini,
                              cmap=plt.cm.RdBu,
                              vmin= -self.GCPs[['residual_x_ini', 'residual_y_ini']].abs().max().max(),
                              vmax= self.GCPs[['residual_x_ini', 'residual_y_ini']].abs().max().max())
        ax[0,0].set_title('Initial Residuals in X')
        plt.colorbar(sc1, ax=ax[0,0])
        sc2 = ax[0,1].scatter(self.GCPs.x_img, self.GCPs.y_img,
                              c=self.GCPs.residual_y_ini,
                              cmap=plt.cm.RdBu,
                              vmin= -self.GCPs[['residual_x_ini', 'residual_y_ini']].abs().max().max(),
                              vmax= self.GCPs[['residual_x_ini', 'residual_y_ini']].abs().max().max())
        ax[0,1].set_title('Initial Residuals in Y')
        plt.colorbar(sc2, ax=ax[0,1])
        sc3 = ax[1,0].scatter(self.GCPs.x_img, self.GCPs.y_img,
                              c=self.GCPs.residual_x_lstsq,
                              cmap=plt.cm.RdBu,
                              vmin= -self.GCPs[['residual_x_lstsq', 'residual_y_lstsq']].abs().max().max(),
                              vmax= self.GCPs[['residual_x_lstsq', 'residual_y_lstsq']].abs().max().max())
        ax[1,0].set_title('LstSq Residuals in X')
        plt.colorbar(sc3, ax=ax[1,0])
        sc4 = ax[1,1].scatter(self.GCPs.x_img, self.GCPs.y_img,
                              c=self.GCPs.residual_y_lstsq,
                              cmap=plt.cm.RdBu,
                              vmin= -self.GCPs[['residual_x_lstsq', 'residual_y_lstsq']].abs().max().max(),
                              vmax= self.GCPs[['residual_x_lstsq', 'residual_y_lstsq']].abs().max().max())
        ax[1,1].set_title('LstSq Residuals in Y')
        plt.colorbar(sc4, ax=ax[1,1])

        sc5 = ax[2,0].scatter(self.GCPs.x_img, 
                              self.GCPs.y_img, 
                              c=self.GCPs.residual_x_ini - self.GCPs.residual_x_lstsq)
        ax[2,0].set_title('Diff in X')
        plt.colorbar(sc5, ax=ax[2,0])
        sc6 = ax[2,1].scatter(self.GCPs.x_img, 
                              self.GCPs.y_img, 
                              c=self.GCPs.residual_y_ini - self.GCPs.residual_y_lstsq)
        ax[2,1].set_title('Diff in Y')
        plt.colorbar(sc6, ax=ax[2,1])