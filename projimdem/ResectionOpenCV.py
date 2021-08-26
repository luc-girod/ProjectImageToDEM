# -*- coding: utf-8 -*-
"""
Created on Thu Aug 26 14:52:25 2021

@author: lucg
"""
import numpy as np
import cv2 as cv
import glob
from matplotlib import pyplot as plt
from scipy import optimize 
import json
import pandas as pd
from types import SimpleNamespace
from math import sin, cos

# termination criteria
criteria = (cv.TERM_CRITERIA_EPS + cv.TERM_CRITERIA_MAX_ITER, 30, 0.001)
# prepare object points, like (0,0,0), (1,0,0), (2,0,0) ....,(6,5,0)
objp = np.zeros((6*7,3), np.float32)
objp[:,:2] = np.mgrid[0:7,0:6].T.reshape(-1,2)
# Arrays to store object points and image points from all the images.
objpoints = [] # 3d point in real world space
imgpoints = [] # 2d points in image plane.
images = glob.glob('*.jpg')
for fname in images:
    img = cv.imread(fname)
    gray = cv.cvtColor(img, cv.COLOR_BGR2GRAY)
    # Find the chess board corners
    ret, corners = cv.findChessboardCorners(gray, (7,6), None)
    # If found, add object points, image points (after refining them)
    if ret == True:
        objpoints.append(objp)
        corners2 = cv.cornerSubPix(gray,corners, (11,11), (-1,-1), criteria)
        imgpoints.append(corners)
        # Draw and display the corners
        cv.drawChessboardCorners(img, (7,6), corners2, ret)
        cv.imshow('img', img)
        cv.waitKey(500)
cv.destroyAllWindows()


GCP_file = '../FinseDemoData/GCPs_pointagev4.csv'
GCPs = pd.read_csv(GCP_file, delimiter=",")
image_file = '../FinseDemoData/2019-05-24_12-00_ori.jpg'
img = cv.imread(image_file)
gray = cv.cvtColor(img, cv.COLOR_BGR2GRAY)

objpoints = [np.concatenate((np.array([GCPs.x_world.loc[:]]).T,np.array([GCPs.y_world.loc[:]]).T,np.array([GCPs.z_world.loc[:]]).T), axis=1).astype('float32')]
imgpoints = [np.concatenate((np.array([GCPs.x_img_prev.loc[:]]).T,np.array([GCPs.y_img_prev.loc[:]]).T), axis=1).astype('float32')]
matrix_init=np.array([[1484,0,0],[0,1484,0],[0,0,1]])
dist_init=np.array([[0,0,0,0,0,0]])
rvecs_init=[np.array([[1.49],[-0.39],[-0.39]])]
tvecs_init=[np.array([[419169.2],[6718421.3],[1212]])]


ret, mtx, dist, rvecs, tvecs = cv.calibrateCamera(objpoints, imgpoints, gray.shape[::-1], int_matrix, dist_init, rvecs_init, tvecs_init, cv.CALIB_USE_INTRINSIC_GUESS)