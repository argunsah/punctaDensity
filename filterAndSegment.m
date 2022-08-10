% Ali Ozgur Argunsah, Zurich, 2022

% Part2: Reads .mat Files (data and cube) and creates _info.mat files
% Ch1: Filters and segments boutons (Otsu + Watershed)
% Ch2: Segments the structure using Otsu

clear all
close all
clc

currentDir = pwd;
addpath(genpath(currentDir));

% Select the data to be analyzed
D1          = uipickfiles; % Select previously extracted Files (.mat Files)

% Change the save folder accordingly
saveFolder  = currentDir;
% or write the full folder name
% i.e. : saveFolder  = '/Users/Argunsah/Documents/punctaDensityResults';

for i = 1:size(D1,2)
    
    load(D1{i}); 
    puncta          = bil_wiener_filt(cube{1},3,1);
    tomato          = cube{2};
    
    thr             = graythresh(puncta(:));
    tomato_thr      = graythresh(tomato(:));
    bw_tomato       = tomato>tomato_thr;

    bw              = (puncta>thr);
    D               = bwdist(~bw);
    D               = -D;
    D(~bw)          = Inf;
    
    L               = watershed(D);
    L(~bw)          = 0;
    
    data.L          = L;
    data.puncta     = puncta;
    data.bw_tomato  = tomato;
    data.tomato     = tomato;
    
    save(fullfile(saveFolder,sprintf('%s_info.mat',data.saveDataName)),'-struct', 'data'); 
end