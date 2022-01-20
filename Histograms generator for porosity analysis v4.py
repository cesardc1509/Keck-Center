# -*- coding: utf-8 -*-
"""
Created on Fri Oct  8 11:31:12 2021

@author: Cesar Diaz
cdiazcarav@miners.utep.edu

Histograms generator for porosity analysis v4
W.M. Keck Center for 3D Innovation

Current version of the code created for automatically processing the porosity 
analysis data of multiple artifacts, as part of the Keck Center efforts on 
the Global Test Artifact Data Exchange Project (GTADExP)

Code takes as inputs the csv file for the pores generated either by MIPAR or 
ImageJ analysis and returns the histograms and statistics for a set of selected
features of interest

Features must be added manually inside the code for now in case more are needed. 
Calling section for functions can be found at the bottom of the code. 

Note: area of the pores must always be the feature at index 1 (starting index 
with 0), otherwise the automatic feature rejection will not work. The order 
in which the other features are placed will not affect the functioning of the 
code.

"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math
import copy


#%% Function for creating the histogram of an specific dataset and feature

def plotting(dataset, bin_size, color, name, parameter_title, xlabel, mean, median, 
             std, limits, xlim, std_xlocation):
    
    # Bin width is provided as parameter for the function. Next line calculates
    # the total number of beans the histogram needs to get that specific width/delta x
    num_bins = math.ceil((np.max(dataset) - np.min(dataset)) / bin_size)
    
    # Props for bbox parameter in the std text box. Matplotlib related stuff.
    props = dict(boxstyle='round', facecolor='white', edgecolor='gray', alpha=0.5)
    
    # Plotting function for the dataset
    ax = plt.subplot()
    plt.hist(dataset, bins = num_bins, color = color)
    plt.title('{} {} Statistical Distrubution'.format(name, parameter_title))
    plt.xlabel(xlabel)
    plt.ylabel('Amount')
    plt.axvline(mean, color = 'black', linestyle = '--', 
                label = 'Mean = {}'.format(round(mean, 2)))
    plt.axvline(median, color = 'blue', linestyle = '--', 
                label = 'Median = {}'.format(round(median, 2)))
    plt.grid(True)
    ax.text(std_xlocation, 0.80,'Standard Deviation: \n {}'.format(round(std, 2)),
                 transform=ax.transAxes, verticalalignment='top', bbox=props)
    
    # Boolean for determining if limits will be included for the histogram. 
    if(limits): 
        plt.xlim(xlim)
    
    plt.legend()
    plt.show()
    
#%% Plotting function for just one feature. 
# Able to plot 3 datasets overlaping over transparent way. Same as the plotting 
# function but without the plt.show and plt.legend commands

def multiple_plotting(dataset, bin_size, color, name, parameter_title, xlabel,
                      limits, xlim, alpha_value):
    
    bins = math.ceil((np.max(dataset) - np.min(dataset)) / bin_size)
    
    # Plotting function for the dataset
    ax = plt.subplot()
    plt.hist(dataset, bins = bins, color = color, label = name, alpha = alpha_value)
    plt.title('{} Statistical Distrubution'.format(parameter_title))
    plt.xlabel(xlabel)
    plt.ylabel('Amount')
    plt.grid(True)
    
    # Boolean for determining if limits will be included for the histogram. 
    if(limits): 
        plt.xlim(xlim)

#%% Function for returning the means, medians, and standard deviations of a 
# dataset. Returns arrays with the statistics of each column in the array given
# as parameter.

def statistics(dataset):
    
    means = dataset.mean(axis = 0)
    medians = np.median(dataset, axis = 0)
    stds = dataset.std(axis = 0)
    
    statistics = np.full((3, len(means)), [means, medians, stds])
    
    return statistics

#%% Function for normalizing the data of a column

def normalize_dataset(dataset):
    
    columns = dataset.shape[1]
    dataset_normalized = np.empty((dataset.shape[0], 0))
    
    for index in range(columns):
        
        # Formula for data normalization: (sample value - mean) / standard deviation
        normalized_stack = (dataset[:, index] - np.mean(dataset[:, index])) / np.std(dataset[:, index])
        
        normalized_stack = np.reshape(normalized_stack, (len(normalized_stack), 1))
        
        dataset_normalized = np.hstack([dataset_normalized, normalized_stack])
        
    return dataset_normalized

#%% Function for eliminating the values of a dataset greater than a certain 
# threshold. Intended for use in histogram plotting.

def cutoff_ratio(dataset, index, cutoff_value, name):
    
    indexes_to_delete = []
    
    for i in range(dataset.shape[0]):
        
        if(dataset[i, index] < cutoff_value):
            
            indexes_to_delete.append(i)
            print('Deleated Index {} from artifact {}'.format(i, name))
            
    dataset = np.delete(dataset, indexes_to_delete, 0)
    
    return dataset
    
#%% Function that uses previously defined functions to return the statistics 
# array and the nomalized dataset of a dataset

def get_and_save_stats(dataset, name):
    
    columns = list(dataset.columns)
    dataset = np.array(dataset)
    
    dataset = cutoff_ratio(dataset, 1, 177, name)# Being 177 the area of a circle of radius 7.5
    df_cutoff = pd.DataFrame(dataset)
    df_cutoff.to_csv('{} Results above 15 microns diameter.csv'.format(name))
    
    statistics_array = statistics(dataset)
    df_statistics_array = pd.DataFrame(statistics_array, index = ['Mean', 'Median', 
                                                              'Standard Deviation'],
                                                               columns = columns)
    df_statistics_array.to_csv('{} statistics array.csv'.format(name))
    
    return dataset, statistics_array

#%% Function for plotting several axes using a for loop

def plotting_file(name_csv, name_artifact):
    
    df_dataset = pd.read_csv(name_csv)
    df_dataset = pd.DataFrame(df_dataset)
    
    columns = df_dataset.columns
    
    dataset, statistics_array = get_and_save_stats(df_dataset, name_artifact)
    columns_number = dataset.shape[1]
    
    i = 1
    
    for i in range(columns_number):
        
        column_name = columns[i]
        
        plot_feature = False
        
        # Matplotlib features definition section 
        if(column_name == 'Area' or column_name == 'Area (um^2)'):
            
            color = 'b'
            parameter = 'Area'
            xlabel = 'Area (μm$^{2}$)'
            limits = True
            xlimits = [0, 5000]
            std_box_location_x = 0.67
            bin_size = 100
            
            plot_feature = True
        
        if(column_name == 'AR' or column_name == 'Aspect Ratio'):
            
            color = 'b'
            parameter = 'Aspect Ratio'
            xlabel = 'Aspect Ratio'
            limits = True
            xlimits = [0, 10]
            std_box_location_x = 0.67
            bin_size = 0.1
            
            plot_feature = True
        
        if(column_name == 'Round' or column_name == 'Roundness'):
            
            color = 'r'
            parameter = 'Roundness'
            xlabel = 'Roundness'
            limits = True
            xlimits = [0, 1]
            std_box_location_x = 0.023
            bin_size = 0.02
            
            plot_feature = True
            
        if(column_name == 'Perim.' or column_name == 'Perimeter'):
            
            color = 'b'
            parameter = 'Perimeter'
            xlabel = 'Perimeter'
            limits = True
            xlimits = [0, 500]
            std_box_location_x = 0.67
            bin_size = 10
            
            plot_feature = True
            
        if(column_name == 'Circularity'):
            
            color = 'g'
            parameter = 'Circularity'
            xlabel = 'Circularity'
            limits = True
            xlimits = [0, 1]
            std_box_location_x = 0.67
            bin_size = 0.02
            
            plot_feature = True
            
        if(column_name == 'Feret'):
            
            color = 'g'
            parameter = 'Maximum Feret Diameter'
            xlabel = 'Maximum Feret Diameter (μm)'
            limits = True
            xlimits = [0, 200]
            std_box_location_x = 0.67
            bin_size = 5
            
            plot_feature = True
            
        if(column_name == 'FeretAngle'):
            
            color = 'r'
            parameter = 'Feret Diameter Angle'
            xlabel = 'Feret Diameter Angle (degrees)'
            limits = True 
            xlimits = [0, 180]
            std_box_location_x = 0.67
            bin_size = 5
            
            plot_feature = True
            
        if(column_name == 'MinFeret'):
            
            color = 'g'
            parameter = 'Minimum Feret Diameter'
            xlabel = 'Minimum Feret Diameter (μm)'
            limits = True
            xlimits = [0, 150]
            std_box_location_x = 0.67
            bin_size = 2
            
            plot_feature = True
            
        if(column_name == 'Solidity'):
            
            color = 'r'
            parameter = 'Solidity'
            xlabel = 'Solidity'
            limits = True
            xlimits = [0, 1]
            std_box_location_x = 0.67
            bin_size = 0.02
            
            plot_feature = True
            
        if(plot_feature):
            
            mean = statistics_array[0, i]
            median = statistics_array[1, i]
            std = statistics_array[2, i]
            
            plotting(dataset[:, i], bin_size = bin_size, color = color, name = name_artifact, 
                     parameter_title = parameter, xlabel = xlabel, mean = mean, 
                     median = median, std = std, limits = limits, xlim = xlimits, 
                     std_xlocation = std_box_location_x)


#%% Function currently not being used

def multiple_plotting_feature(names_files, names_artifacts, feature):
    
    i = 0
    alpha_value = 0.3
    
    for i in range(len(names_files)):
        
        df_dataset = pd.read_csv(names_files[i])
        df_dataset = pd.DataFrame(df_dataset)
        
        columns = df_dataset.columns
        name_artifact = names_artifacts[i]
        
        dataset, statistics_array = get_and_save_stats(df_dataset, name_artifact)
        columns_number = dataset.shape[1]
        
        i = 1
        
        for i in range(columns_number):
            
            column_name = columns[i]
            plot_feature = False
            
            # Feature selection for matplotlib
            if(column_name == 'Area' or column_name == 'Area (um^2)'):
                
                color = 'b'
                parameter = 'Area'
                xlabel = 'Area (μm$^{2}$)'
                limits = True
                xlimits = [0, 5000]
                bin_size = 100
                
                plot_feature = True
                
            if(plot_feature):
                
                multiple_plotting(dataset[:, i], bin_size = bin_size, color = color, 
                                  name = name_artifact, parameter_title = parameter, 
                                  xlabel = xlabel, limits = limits, xlim = xlimits, 
                                  alpha_value = alpha_value)
    
    return -1 #For people not familiar with software engineering, -1 is a 
              # commonly returned value when a function is not being used

#%% Function for plotting the features of multiple data sets overlapping and 
# in the same plot (with trasnparency)

def multiple_plotting_file(names_files, names_artifacts):
    
    i = 0
    alpha_value = 0.5
    
    color_counter = 0
    
    for i in range(len(names_files)):
        
        df_dataset = pd.read_csv(names_files[i])
        df_dataset = pd.DataFrame(df_dataset)
        
        columns = df_dataset.columns
        name_artifact = names_artifacts[i]
        
        dataset, statistics_array = get_and_save_stats(df_dataset, name_artifact)
        columns_number = dataset.shape[1]
        
        i = 1
        
        for i in range(columns_number):
            
            column_name = columns[i]
            #print(i)
            
            plot_feature = False
            
            # Feature selection for matplotlib
            if(column_name == 'Area' or column_name == 'Area (um^2)'):
                
                if(color_counter == 0):
                    
                    color = 'b'
                
                elif(color_counter == 1):
                    
                    color = 'r'
                
                elif(color_counter == 2):
                    
                    color = 'g'
                    
                else:
                    
                    color = 'b'
                    color_counter = 1
                
                parameter = 'Area'
                xlabel = 'Area (μm$^{2}$)'
                limits = True
                xlimits = [0, 5000]
                bin_size = 100
                
                plot_feature = True
                
                color_counter += 1
                
            if(plot_feature):
                
                multiple_plotting(dataset[:, i], bin_size = bin_size, color = color, 
                                  name = name_artifact, parameter_title = parameter, 
                                  xlabel = xlabel, limits = limits, xlim = xlimits, 
                                  alpha_value = alpha_value)
    
    plt.legend()
    plt.show()
    
    color_counter = 0
    
    for i in range(len(names_files)):
        
        df_dataset = pd.read_csv(names_files[i])
        df_dataset = pd.DataFrame(df_dataset)
        
        columns = df_dataset.columns
        name_artifact = names_artifacts[i]
        
        dataset, statistics_array = get_and_save_stats(df_dataset, name_artifact)
        columns_number = dataset.shape[1]
        
        i = 1
        
        for i in range(columns_number):
            
            column_name = columns[i]
            #print(i)
            
            plot_feature = False
            
            # Feature selection for matplotlib
            if(column_name == 'AR' or column_name == 'Aspect Ratio'):
            
                if(color_counter == 0):
                    
                    color = 'b'
                
                elif(color_counter == 1):
                    
                    color = 'r'
                
                elif(color_counter == 2):
                    
                    color = 'g'
                    
                else:
                    
                    color = 'b'
                    color_counter = 1
            
                parameter = 'Aspect Ratio'
                xlabel = 'Aspect Ratio'
                limits = True
                xlimits = [0, 10]
                bin_size = 0.1
                
                plot_feature = True
                
                color_counter += 1
                
            if(plot_feature):
                
                multiple_plotting(dataset[:, i], bin_size = bin_size, color = color, 
                                  name = name_artifact, parameter_title = parameter, 
                                  xlabel = xlabel, limits = limits, xlim = xlimits, 
                                  alpha_value = alpha_value)
        
    plt.legend()
    plt.show()
    
    color_counter = 0
    
    for i in range(len(names_files)):
        
        df_dataset = pd.read_csv(names_files[i])
        df_dataset = pd.DataFrame(df_dataset)
        
        columns = df_dataset.columns
        name_artifact = names_artifacts[i]
        
        dataset, statistics_array = get_and_save_stats(df_dataset, name_artifact)
        columns_number = dataset.shape[1]
        
        i = 1
        
        for i in range(columns_number):
            
            column_name = columns[i]
            #print(i)
            
            plot_feature = False
            
            # Feature selection for matplotlib
            if(column_name == 'Round' or column_name == 'Roundness'):
                
                if(color_counter == 0):
                    
                    color = 'b'
                
                elif(color_counter == 1):
                    
                    color = 'r'
                
                elif(color_counter == 2):
                
                    color = 'g'
                    
                else:
                    
                    color = 'b'
                    color_counter = 1    
                
                parameter = 'Roundness'
                xlabel = 'Roundness'
                limits = True
                xlimits = [0, 1]
                bin_size = 0.02
            
                plot_feature = True
                
                color_counter += 1
            
            if(plot_feature):
                
                multiple_plotting(dataset[:, i], bin_size = bin_size, color = color, 
                                  name = name_artifact, parameter_title = parameter, 
                                  xlabel = xlabel, limits = limits, xlim = xlimits, 
                                  alpha_value = alpha_value)
    
    plt.legend()
    plt.show()
    
    color_counter = 0
    
    for i in range(len(names_files)):
        
        df_dataset = pd.read_csv(names_files[i])
        df_dataset = pd.DataFrame(df_dataset)
        
        columns = df_dataset.columns
        name_artifact = names_artifacts[i]
        
        dataset, statistics_array = get_and_save_stats(df_dataset, name_artifact)
        columns_number = dataset.shape[1]
        
        i = 1
        
        for i in range(columns_number):
            
            column_name = columns[i]
            #print(i)
            
            plot_feature = False
            
            # Feature selection for matplotlib
            if(column_name == 'Perim.' or column_name == 'Perimeter (um)' or 
               column_name == 'Perimeter'):
                
                if(color_counter == 0):
                    
                    color = 'b'
                
                elif(color_counter == 1):
                    
                    color = 'r'
                
                elif(color_counter == 2):
                    
                    color = 'g'
                    
                else:
                    
                    color = 'b'
                    color_counter = 1    
            
                parameter = 'Perimeter'
                xlabel = 'Perimeter'
                limits = True
                xlimits = [0, 500]
                bin_size = 10
                
                plot_feature = True
                
                color_counter += 1
                
            if(plot_feature):
                
                multiple_plotting(dataset[:, i], bin_size = bin_size, color = color, 
                                  name = name_artifact, parameter_title = parameter, 
                                  xlabel = xlabel, limits = limits, xlim = xlimits, 
                                  alpha_value = alpha_value)
    
    plt.legend()
    plt.show()
    
    color_counter = 0
    
    for i in range(len(names_files)):
        
        df_dataset = pd.read_csv(names_files[i])
        df_dataset = pd.DataFrame(df_dataset)
        
        columns = df_dataset.columns
        name_artifact = names_artifacts[i]
        
        dataset, statistics_array = get_and_save_stats(df_dataset, name_artifact)
        columns_number = dataset.shape[1]
        
        i = 1
        
        for i in range(columns_number):
            
            column_name = columns[i]
            #print(i)
            
            plot_feature = False
            
            # Feature selection for matplotlib
            if(column_name == 'Circularity'):
                
                if(color_counter == 0):
                    
                    color = 'b'
                
                elif(color_counter == 1):
                    
                    color = 'r'
                
                elif(color_counter == 2):
                
                    color = 'g'
                    
                else:
                    
                    color = 'b'
                    color_counter = 1    
                
                parameter = 'Circularity'
                xlabel = 'Circularity'
                limits = True
                xlimits = [0, 1]
                std_box_location_x = 0.023
                bin_size = 0.02
                
                plot_feature = True
                
                color_counter += 1
            
            if(plot_feature):
                
                multiple_plotting(dataset[:, i], bin_size = bin_size, color = color, 
                                  name = name_artifact, parameter_title = parameter, 
                                  xlabel = xlabel, limits = limits, xlim = xlimits, 
                                  alpha_value = alpha_value)
    
    plt.legend()
    plt.show()
    
    color_counter = 0
    
    for i in range(len(names_files)):
        
        df_dataset = pd.read_csv(names_files[i])
        df_dataset = pd.DataFrame(df_dataset)
        
        columns = df_dataset.columns
        name_artifact = names_artifacts[i]
        
        dataset, statistics_array = get_and_save_stats(df_dataset, name_artifact)
        columns_number = dataset.shape[1]
        
        i = 1
        
        for i in range(columns_number):
            
            column_name = columns[i]
            #print(i)
            
            plot_feature = False
            
            # Feature selection for matplotlib
            if(column_name == 'Feret'):
                
                if(color_counter == 0):
                    
                    color = 'b'
                
                elif(color_counter == 1):
                    
                    color = 'r'
                
                elif(color_counter == 2):
                
                    color = 'g'
                    
                else:
                    
                    color = 'b'
                    color_counter = 1    
                
                parameter = 'Maximum Feret Diameter'
                xlabel = 'Maximum Feret Diameter (μm)'
                limits = False
                xlimits = [0, 1]
                bin_size = 5
                
                plot_feature = True
                
                color_counter += 1
            
            if(plot_feature):
                
                multiple_plotting(dataset[:, i], bin_size = bin_size, color = color, 
                                  name = name_artifact, parameter_title = parameter, 
                                  xlabel = xlabel, limits = limits, xlim = xlimits, 
                                  alpha_value = alpha_value)
    
    plt.legend()
    plt.show()
    
    color_counter = 0
    
    for i in range(len(names_files)):
        
        df_dataset = pd.read_csv(names_files[i])
        df_dataset = pd.DataFrame(df_dataset)
        
        columns = df_dataset.columns
        name_artifact = names_artifacts[i]
        
        dataset, statistics_array = get_and_save_stats(df_dataset, name_artifact)
        columns_number = dataset.shape[1]
        
        i = 1
        
        for i in range(columns_number):
            
            column_name = columns[i]
            #print(i)
            
            plot_feature = False
            
            # Feature selection for matplotlib
            if(column_name == 'FeretAngle'):
                
                if(color_counter == 0):
                    
                    color = 'b'
                
                elif(color_counter == 1):
                    
                    color = 'r'
                
                elif(color_counter == 2):
                
                    color = 'g'
                
                else:
                    
                    color = 'b'
                    color_counter = 1    
                
                parameter = 'Maximum Feret Diameter Angle'
                xlabel = 'Maximum Feret Diameter Angle (degrees)'
                limits = False
                xlimits = [0, 200]
                bin_size = 5
                
                plot_feature = True
                
                color_counter += 1
            
            if(plot_feature):
                
                multiple_plotting(dataset[:, i], bin_size = bin_size, color = color, 
                                  name = name_artifact, parameter_title = parameter, 
                                  xlabel = xlabel, limits = limits, xlim = xlimits, 
                                  alpha_value = alpha_value)
    
    plt.legend()
    plt.show()
    
    color_counter = 0
    
    for i in range(len(names_files)):
        
        df_dataset = pd.read_csv(names_files[i])
        df_dataset = pd.DataFrame(df_dataset)
        
        columns = df_dataset.columns
        name_artifact = names_artifacts[i]
        
        dataset, statistics_array = get_and_save_stats(df_dataset, name_artifact)
        columns_number = dataset.shape[1]
        
        i = 1
        
        for i in range(columns_number):
            
            column_name = columns[i]
            #print(i)
            
            plot_feature = False
            
            # Feature selection for matplotlib
            if(column_name == 'MinFeret'):
                
                if(color_counter == 0):
                    
                    color = 'b'
                
                elif(color_counter == 1):
                    
                    color = 'r'
                
                elif(color_counter == 2):
                
                    color = 'g'
                    
                else:
                    
                    color = 'b'
                    color_counter = 1    
                
                parameter = 'Minimum Feret Diameter'
                xlabel = 'Minimum Feret Diameter (μm)'
                limits = True
                xlimits = [0, 300]
                std_box_location_x = 0.023
                bin_size = 5
            
                plot_feature = True
                
                color_counter += 1
            
            if(plot_feature):
                
                multiple_plotting(dataset[:, i], bin_size = bin_size, color = color, 
                                  name = name_artifact, parameter_title = parameter, 
                                  xlabel = xlabel, limits = limits, xlim = xlimits, 
                                  alpha_value = alpha_value)
    
    plt.legend()
    plt.show()
    
    color_counter = 0
    
    for i in range(len(names_files)):
        
        df_dataset = pd.read_csv(names_files[i])
        df_dataset = pd.DataFrame(df_dataset)
        
        columns = df_dataset.columns
        name_artifact = names_artifacts[i]
        
        dataset, statistics_array = get_and_save_stats(df_dataset, name_artifact)
        columns_number = dataset.shape[1]
        
        i = 1
        
        for i in range(columns_number):
            
            column_name = columns[i]
            #print(i)
            
            plot_feature = False
            
            # Feature selection for matplotlib
            if(column_name == 'Solidity'):
                
                if(color_counter == 0):
                    
                    color = 'b'
                
                elif(color_counter == 1):
                    
                    color = 'r'
                
                elif(color_counter == 2):
                
                    color = 'g'
                    
                else:
                    
                    color = 'b'
                    color_counter = 1    
                
                parameter = 'Solidity'
                xlabel = 'Solidity'
                limits = True
                xlimits = [0, 1]
                std_box_location_x = 0.023
                bin_size = 0.02
            
                plot_feature = True
                
                color_counter += 1
            
            if(plot_feature):
                
                multiple_plotting(dataset[:, i], bin_size = bin_size, color = color, 
                                  name = name_artifact, parameter_title = parameter, 
                                  xlabel = xlabel, limits = limits, xlim = xlimits, 
                                  alpha_value = alpha_value)
    
    plt.legend()
    plt.show()
    
    color_counter = 0
    
#%% Function Calling Section 

plotting_file('QTA0110.csv', 'QTA 0110')

'''
files_for_plotting = ['QTA 0109 manual.csv', 'QTA 0110 manual.csv', 'QTA 0111 manual.csv']
names_artifacts = ['QTA 0109 with ImageJ', 'QTA 0110 with ImageJ', 'QTA 0111 with ImageJ']

multiple_plotting_file(files_for_plotting, names_artifacts)

files_for_plotting_10_11 = ['QTA 0110 manual.csv', 'QTA 0111 manual.csv']
names_artifacts_10_11 = ['QTA 0110 with ImageJ', 'QTA 0111 with ImageJ']

multiple_plotting_file(files_for_plotting_10_11, names_artifacts_10_11)'''

