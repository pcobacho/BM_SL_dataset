% *************************************************************************
% AUTHORS: Francisco J. Martin Vega
% *************************************************************************
% GROUP: Lab 1.3.5., Communications and Signal Processing Lab (ComSP), 
% TELMA, ETSIT, Universiry of Malaga
% *************************************************************************
% DESCRIPTION:
% This script adds paths to code folders. 
% OUTPUTs:
% +) currentDir: It is the set of working paths including the added paths
% +) oldPath: It is the set of working paths before adding those of the
% project. 
% *************************************************************************

function [currentDir, oldPath] = addPaths()
oldPath = path;
currentDir = cd;
path([currentDir '\Lib'], path);
cd(currentDir);