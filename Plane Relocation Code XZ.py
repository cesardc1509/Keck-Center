# -*- coding: utf-8 -*-
"""
Created on Mon Nov  8 10:45:45 2021

@authors: Cesar Diaz, Jaime Mata
cdiazcarav@miners.utep.edu, jrmata2@miners.utep.edu

W.M. Keck Center for 3D Innovation
The University of Texas at El Paso

Code for plane relocation of artifacts, as part of the Keck Center efforts on 
the Global Test Artifact Data Exchange Program (GTADEXP). The code takes as 
input the raw coordinates as read in MIPAR using the geometry features the 
artifact includes, and the pore information as exported from MIPAR or ImageJ.
In an original picture as returned from the Keyance VHX 7000 microscope, the 
origin is located in the top left coordinate of the image. This code changes
this origin to a position near the central lattice that is standard for each 
artifact, and changes the coordinates' units from pixels to micrometers.

The original code was written in MATLAB and was separated in two scripts, one 
containing the process for calculating the equation of the plane, and one for 
converting the pixel coordinates using the transformation factors obtained from 
the plane equation code (called DZPlaneCalc.m and PoreLocations.m respectively).
Functions from coordinate arragement to STEP 9B correspond to the equation of 
the plane calculation process. Starting from STEP 10 they correspond to the 
process of converting the coordinates. 

"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import numpy as np
import math

# %% Data input section - Only section that needs to be modified by the user

# Please introduce below the name of the file containing the data obtained from
# MIPAR or ImageJ, the name of the input coordinates file obtained for the
# pore relocation process, the column number on the pores Excel file containing
# the x and y coordinates for each pore, the total number of columns on the
# pores Excel file, and the name of the artifact being analyzed.

# Both input coordinates and pore information must be in the same directory
# as the code for this to work.

pores_data_name = 'QTA 10 Results - From ROI.xlsx'
triangles_data_name = 'InputCorQTA0110.xlsx'
x_coords_index = 2  # For column numbers starting at 0
y_coords_index = 3
total_columns = 16  # Total number of columns on the pores Excel file
name_artifact = 'QTA0110'

# %% Coordinate arrangement

def arrange_coords(df_coords):
    # Function takes as imput a pandas DataFrame with the set of coordinates
    # obtained from MIPAR for an specific artifact and returns those
    # coordinates arranged in two 1D arrays, for lower and upper coordinates.
    # Those 1D arrays contain the coordinates on x and z points for each point.

    coords = np.array(df_coords[['X points', 'Z Points']])

    upper_coords = coords[0:4, :]
    lower_coords = coords[4:8, :]

    upper = np.zeros([8])
    lower = np.zeros([8])

    i = 0
    j = 0

    # Rearangement of upper coordintates
    for i in range(upper_coords.shape[0]):
        upper[j] = upper_coords[i, 0]
        upper[j + 1] = upper_coords[i, 1]
        j = j + 2

    i = 0
    j = 0

    # Rearangement of lower coordinates
    for i in range(lower_coords.shape[0]):
        lower[j] = lower_coords[i, 0]
        lower[j + 1] = lower_coords[i, 1]
        j = j + 2

    return upper, lower

# %% Conversion to microns

def convert_to_microns(upper_pixels, lower_pixels):
    # Necessary conversion to microns, function just multiplies the pixels
    # array times the conversion factor micometers/pixel, which is obtained
    # from the raw image taken in the Keyance VHX 7000 Digital Microscope

    conversion_factor = 2.06  # Micrometers/pixel

    upper_microns = upper_pixels * conversion_factor
    lower_microns = lower_pixels * conversion_factor

    return upper_microns, lower_microns

# %% STEP 1 - Calculation of fiducial parameters (a', D'...)

def get_a_and_d(upper_points, lower_points):
    # T stands for top and B stands for bottom. This is linked with the location of the
    # triangles
    # A is the distance between the points gathered in the x direction
    # R stands for right and L for left
    # D is the distance between the triangles

    top_a_left = upper_points[0] - upper_points[2]
    top_a_right = upper_points[4] - upper_points[6]
    bottom_a_left = lower_points[0] - lower_points[2]
    bottom_a_right = lower_points[4] - lower_points[6]

    a = np.array([top_a_left, top_a_right, bottom_a_left, bottom_a_right])

    dp_top = upper_points[2] - upper_points[4]
    dp_bottom = lower_points[2] - lower_points[4]
    d_prime = [dp_top, dp_bottom]

    return a, d_prime

# %% STEP 2 - Calculation of Theta for both the upper and lower section

def get_theta(distance):
    global gamma
    global d_top
    global d_bottom

    top_th = math.asin((d_top / distance[0]) * math.sin(gamma)) - gamma
    bottom_th = math.asin((d_bottom / distance[1]) * math.sin(gamma)) - gamma

    return top_th, bottom_th

# %% STEP 3 - Calculation of A Prime

def cal_a_prime(a, top_th, bottom_th):
    global gamma

    a_prime_l_top = a[0] * (math.sin(gamma + top_th)) / math.sin(gamma)
    a_prime_r_top = a[1] * (math.sin(gamma + top_th)) / math.sin(gamma)
    a_prime_l_bottom = a[2] * (math.sin(gamma + bottom_th)) / math.sin(gamma)
    a_prime_r_bottom = a[3] * (math.sin(gamma + bottom_th)) / math.sin(gamma)

    a_prime_array = np.array([a_prime_l_top, a_prime_r_top, a_prime_l_bottom,
                              a_prime_r_bottom])

    return a_prime_array

# %% STEP 4 - Calculation of H

def get_h(a_prime_array):
    global gamma

    h_top_l = (a_prime_array[0] / 2) * math.tan(gamma)
    h_top_r = (a_prime_array[1] / 2) * math.tan(gamma)
    h_bottom_l = (a_prime_array[2] / 2) * math.tan(gamma)
    h_bottom_r = (a_prime_array[3] / 2) * math.tan(gamma)

    h_array = np.array([h_top_l, h_top_r, h_bottom_l, h_bottom_r])

    return h_array

# %% STEP 5 - Calculation of Y Coordinates for the four points of interest

def get_y(upper, lower):
    # The actual Y coordinates of the artifact after the cut are calculated. The
    # origin and the zero value are located in the middle of the artifact.
    # Therefore, since the cut made is not perfectly aligned nor is it
    # completely straight, the actual y coordinate in the image needs to be
    # acquired.
    # Fixed Values are the following:

    global d_top
    d_top = 8500  # Distance between upper traingles in microns

    global d_bottom
    d_bottom = 22500  # Distance between lower triangles in microns

    global gamma
    gamma = 1.106538746  # This value is in radians

    base_triangle = 3000  # Triangle base value in microns

    # First Step: Convert points provided into microns using the 2.06 micron/pixel
    # conversion factor.
    upper_points, lower_points = convert_to_microns(upper, lower)
    As, DPs = get_a_and_d(upper_points, lower_points)
    top_th, bottom_th = get_theta(DPs)
    APs = cal_a_prime(As, top_th, bottom_th)
    h = get_h(APs)

    y_left_top = -1500 + h[0]
    y_right_top = 1500 - h[1]
    y_left_bottom = -1500 + h[2]
    y_right_bottom = 1500 - h[3]

    y = np.array([y_left_top, y_right_top, y_left_bottom, y_right_bottom])

    return y

# %% STEP 6 - Vector Calculation from Bottom Left for Plane Calculation

def vector_calculation(upper, lower):
    # Fixed Values of Location in CAD. The order is from left to right and from
    # top to bottom. Values are given in microns

    x_values = np.array([-5000, 5000, -12000, 12000])
    z_values = np.array([22000, 22000, -11500, -11500])

    y = get_y(upper, lower)

    global top_right
    top_right = np.array([x_values[1], y[1], z_values[1]])
    global bottom_left
    bottom_left = np.array([x_values[2], y[2], z_values[2]])
    global bottom_right
    bottom_right = np.array([x_values[3], y[3], z_values[3]])

    bottom_left_to_top_right = np.array([(top_right[0] - bottom_left[0]),
                                         (top_right[1] - bottom_left[1]),
                                         (top_right[2] - bottom_left[2])])
    bottom_left_to_bottom_right = np.array([(bottom_right[0] - bottom_left[0]),
                                            (bottom_right[1] - bottom_left[1]),
                                            (bottom_right[2] - bottom_left[2])])

    return bottom_left_to_top_right, bottom_left_to_bottom_right

# %% STEP 7A - Plane Calculation using unit vectors

def unit_plane_calculation(upper_points, lower_points):
    # Normalized plane calculation is done which takes into account the
    # previously calulated actual y-coordinates

    u, v = vector_calculation(upper_points, lower_points)
    global top_right
    normal_dis_u = math.sqrt(u[0] ** 2 + u[1] ** 2 + u[2] ** 2)
    normal_dis_v = math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)

    unit_vector_u = u / normal_dis_u
    unit_vector_v = v / normal_dis_v

    a_u = (unit_vector_u[1] * unit_vector_v[2]) - (unit_vector_u[2] *
                                                   unit_vector_v[1])
    b_u = (unit_vector_u[2] * unit_vector_v[0]) - (unit_vector_v[2] *
                                                   unit_vector_u[0])
    c_u = (unit_vector_u[0] * unit_vector_v[1]) - (unit_vector_u[1] *
                                                   unit_vector_v[0])
    d_u = (a_u * top_right[0] + b_u * top_right[1] + c_u * top_right) * -1

    d_u_1 = d_u[2]

    abcd_u = np.array([a_u, b_u, c_u, d_u_1])

    return abcd_u

# %% STEP 7B - Plane Calculaion without unit vectors

def plane_calculator(upper_points, lower_points):
    # This function is to calculate the plane without using unit vectors

    global top_right
    bottom_left_top_right, bottom_left_bottom_right = vector_calculation(upper_points,
                                                                         lower_points)
    bltr = bottom_left_top_right
    blbr = bottom_left_top_right

    a = (bltr[1] * blbr[2]) - (bltr[2] * blbr[1])
    b = (bltr[2] * blbr[0]) - (blbr[2] * bltr[0])
    c = (bltr[0] * blbr[1]) - (bltr[1] * blbr[0])
    dp = (a * top_right[0] + b * top_right[1] + c * top_right)

    d = dp[2]

    abcd = np.array([a, b, c, d])

    return abcd

# %% STEP 8 - Calculation of Angle of Rotation

def get_angle(uc, lc):
    # UC stands for upper coordinates and LC stands for lower coordinates
    # Angle is given in radians. Using a line fitting function, the average of
    # the slopes calculated is acquired to determine the angle of tilt of the
    # image

    # Triple Indexing in Python: [start_point:end_point(Not Inclusive):increment]
    # Triple Indexing in MATLAB: [start_point:incremet:end_point(Inclusive)]
    upper = np.polyfit(uc[0:8:2], uc[1:8:2], 1)
    lower = np.polyfit(lc[0:8:2], lc[1:8:2], 1)

    avg = (upper[0] + lower[0]) / 2
    angle = math.atan(avg)

    return angle

# %% STEP 9 - Calculation of Transfomation Factors

def get_trans_factors(upper, lower):
    # Transformation factors needed to translate coordinates within the
    # previously calculated plane are acquired. The center of the left lower
    # traingle is calculated in microns

    upper_microns, lower_microns = convert_to_microns(upper, lower)
    y_coordinates = get_y(upper, lower)
    abcd = unit_plane_calculation(upper, lower)
    angle1 = get_angle(upper, lower)
    x = (lower_microns[0] + lower_microns[2]) / 2
    y1 = (-abcd[0] * lower_microns[0] - abcd[2]
          * lower_microns[1] - abcd[3]) / abcd[1]
    y2 = (-abcd[0] * lower_microns[2] - abcd[2]
          * lower_microns[3] - abcd[3]) / abcd[1]

    y = (y1 + y2) / 2
    z = (lower_microns[1] + lower_microns[3]) / 2

    x1 = x * math.cos(angle1) + z * math.sin(angle1)
    z1 = -x * math.sin(angle1) + z * math.cos(angle1)
    x2 = -12000 - x1
    y2 = y_coordinates[2] - y
    z3 = -11500 - z1

    trans_factors = np.array([x2, y2, z3])

    return trans_factors

# %% STEP 9B - Calculation of Transformation Factor Without Unit Vectors

def get_trans_factors_without_unit_vector(upper, lower):

    upper_microns, lower_microns = convert_to_microns(upper, lower)
    y_coordinates = get_y(upper, lower)
    abcd = plane_calculator(upper, lower)
    angle1 = get_angle(upper, lower)

    x = (lower_microns[0] + lower_microns[2]) / 2
    y1 = (-abcd[0] * lower_microns[0] - abcd[2]
          * lower_microns[1] - abcd[3]) / abcd[1]
    y2 = (-abcd[0] * lower_microns[2] - abcd[2]
          * lower_microns[3] - abcd[3]) / abcd[1]
    y = (y1 + y2) / 2
    z = (lower_microns[1] + lower_microns[3]) / 2
    x1 = x * math.cos(angle1) + z * math.sin(angle1)
    z1 = -x * math.sin(angle1) + z * math.cos(angle1)
    x2 = -12000 - x1
    y2 = y_coordinates[2] - y
    z3 = -11500 - z1

    trans_factors_no_unit = np.array([x2, y2, z3])

    return trans_factors_no_unit

# %% STEP 10 - Calculation of Y for each Pore Location

def get_y_pores(x_z, abcd):

    length_x_z_rows, length_x_z_cols = np.shape(x_z)

    y = np.ones([length_x_z_rows, 1])
    i = 0

    for i in range(length_x_z_rows):
        y[i] = (-abcd[0] * x_z[i, 0] - abcd[2] * x_z[i, 1] - abcd[3]) / abcd[1]

    return y

# %% STEP 11 - Rotation of X and Z

def get_rotation(x_z, angle):

    length_x_z_rows, length_x_z_cols = np.shape(x_z)

    x_z_rot = np.ones([length_x_z_rows, 2])
    i = 0

    for i in range(length_x_z_rows):
        x_z_rot[i, 0] = x_z[i, 0] * \
            math.cos(angle) + x_z[i, 1] * math.sin(angle)
        x_z_rot[i, 1] = -x_z[i, 0] * \
            math.sin(angle) + x_z[i, 1] * math.cos(angle)

    return x_z_rot

# %% STEP 12 - Translating Coordinates

def get_translation(coords, translation_factor):

    length_coords_rows, length_coords_cols = np.shape(coords)

    coords_translated = np.ones([length_coords_rows, 3])
    i = 0

    for i in range(length_coords_rows):
        coords_translated[i, 0] = coords[i, 0] + translation_factor[0]
        coords_translated[i, 1] = coords[i, 1] + translation_factor[1]
        coords_translated[i, 2] = coords[i, 2] + translation_factor[2]

    return coords_translated

# %% Variable Output Printing

def output_printing(name, variable):
    # Simple function for printing the output of a variable. Performs a similar
    # task to MATLAB printing. Function used mainly for testing outputs.

    print('{} = {} \n'.format(name, variable))

# %% DZ Plane Calculation (Section calling all the functions for the DZ Plane
#   Calculation)

data = pd.read_excel(triangles_data_name)
data = pd.DataFrame(data)

upper, lower = arrange_coords(data)

bltr, blbr = vector_calculation(upper, lower)

transformation_factors = get_trans_factors(upper, lower)
angle1 = get_angle(upper, lower)
abcd = unit_plane_calculation(upper, lower)

output_printing('Upper', upper)
output_printing('Lower', lower)
output_printing('abcd', abcd)
output_printing('transformation_factors', transformation_factors)
output_printing('angle1', angle1)

# %% Pore Locations (Section calling all the functions and performing operations
#   for obtaining the transformed pore locations)

# Getting input data
pores = pd.read_excel(pores_data_name)
pores = pd.DataFrame(pores)

column_names = pores.columns
pores = np.array(pores)

# Getting x and y coordinates on pores info
x_z = pores[:, [x_coords_index, y_coords_index]]

# Performing Operations
x_z[:, 1] = x_z[:, 1] * -1
output_printing('x_z', x_z)

y = get_y_pores(x_z, abcd)
output_printing('y', y)

x_z_2 = get_rotation(x_z, angle1)

length_xz2_rows, length_xz2_cols = np.shape(x_z_2)
length_y_rows, length_y_cols = np.shape(y)

output_printing('xz2 rows', length_xz2_rows)
output_printing('xz2 columns', length_xz2_cols)
output_printing('length y rows', length_y_rows)
output_printing('length y columns', length_y_cols)

# The reshape is added to the columns of x_z_2 because hstack can only add
# arrays of the same number of dimensions. Indexing a single column on x_z_2
# originally creates a one dimensional array. The reshape function turns it on
# a 2D array
coords = np.hstack((x_z_2[:, 0].reshape((length_xz2_rows, 1)), y,
                    x_z_2[:, 1].reshape((length_xz2_rows, 1))))

# Final translation of coordinates
trans_coords_mic = get_translation(coords, transformation_factors)
output_printing('Transformed Coordinates', trans_coords_mic)

# Table with area and X, Y, Z Coordinates
axyz = np.hstack((pores[:, 0:x_coords_index], trans_coords_mic / 1000,
                  pores[:, y_coords_index + 1:total_columns]))

# Getting the column arrays for the outputing file
column_names = np.hstack((column_names[0:y_coords_index + 1], 'Z',
                          column_names[y_coords_index + 1:total_columns]))

# Output Excel Sheet with Pore XYZ Coordinates
df_axyz = pd.DataFrame(axyz, columns=column_names)
df_axyz.to_csv('{} Changed Coords.csv'.format(name_artifact))

# %% Plotting Values
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(axyz[:, x_coords_index], axyz[:, y_coords_index],
           axyz[:, y_coords_index + 1], color='b')
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')
plt.show()

plt.figure(figsize=(11, 10))
plt.scatter(axyz[:, x_coords_index], axyz[:, y_coords_index + 1], s=10)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlim([-22, 22])
plt.ylim([-15, 25])
plt.show()

