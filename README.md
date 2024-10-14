# Altimetric-Data-Check


This repository contains:
1. River water level altimetric data from multiple sources for two basins (Niger and Ganges-Brahmaputra) compressed as .tar.xz file.
2. A MATLAB script for performing a simple analysis and visualization of this data.
3. Python scripts for downloading altimetric data from each data source.


## Latest MATLAB Scripts

Visualizes the time series of water levels at all virtual stations. Plots the stations' locations on a map in two different ways:
  1. A scatter plot where the color represents the mean water level at each station.
  2. A map that separates stations based on the source of the altimetric data.

V07: user can specify path to each data in the "Specify Paths" section

default paths are mentioned below, where REGION is either Niger or Ganges-Brahmaputra:                                                                     
Dahiti Altimetry Data: './Dahiti/REGION Basin/'                                                                     
Hydroweb Altimetry Data: './Hydroweb/REGION Basin/'                                                                      
CLMS Altimetry data: './CLMS/REGION Basin/'                                                                      


V08: opens a window for each data source to select desired data. a data source can be ignored by closing it's window without selecting any data.
