import os
directory="/mnt/e/Boternosa"
for filename in os.listdir(directory):
    f = os.path.join(directory, filename)
    # checking if it is a file
    if (os.path.isfile(f) and f.endswith("jpg")):
        print(f)
        image_file     = f
        output_file    = '/mnt/e/Boternosa/Orthos/' + filename.removesuffix(".jpg") + '_Ortho.tif'
        print(output_file)
        Bolternosa_proj = pi.Projection(dem_file=dem_file,
                          viewshed_file=viewshed_file,
                          image_file=image_file,
                          cam_param=Bolternosa.new_cam.proj_param,
                           output_file=output_file
                          )
        Bolternosa_proj.project_img_to_DEM(return_raster=True, epsg=32633)



import os
directory="/mnt/e/FinseWebcam"
os.mkdir('/mnt/e/FinseWebcam/Orthos/')
for filename in os.listdir(directory):
    f = os.path.join(directory, filename)
    # checking if it is a file
    if (os.path.isfile(f) and f.endswith("jpg")):
        print(f)
        image_file     = f
        output_file    = '/mnt/e/FinseWebcam/Orthos/' + filename.removesuffix(".jpg") + '_Ortho.tif'
        print(output_file)
        Finse_proj = pi.Projection(dem_file=dem_file,
                          viewshed_file=viewshed_file,
                          image_file=image_file,
                          cam_param=finse.new_cam.proj_param,
                           output_file=output_file
                          )
        Finse_proj.project_img_to_DEM(return_raster=True, epsg=32632)

