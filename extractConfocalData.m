% Ali Ozgur Argunsah, Zurich, 2022

% Part1: This part is for automatic extraction of Confocal Data
% Reads Confocal Files and creates .mat files
% 2 Channels
% Ch1: Puncta (GFP)
% Ch2: Structure (td-Tomato)

clear all
close all
clc

currentDir = pwd;
addpath(genpath(currentDir));

% Select the data to be analyzed
D1          = uipickfiles; % Select rawConfocal Files (i.e. .oif Files, not Folders)

for i = 1:size(D1,2)
    close all;
    
    % Create Info Structure
    clear data;
    data  = struct;

    [fileDir,fileName,fileExt]  = fileparts(D1{i});
    saveDataName                = fileName;
    saveFolder                  = currentDir;
    
    data.name           = D1{i};
    data.analysisday    = date;
    data.Analyzer       = 'Automatic';  % What is the name of the person who analyzed the data
    data.DataCollected  = 'Lorenzo';    % Who collected the data?
    data.saveDataName   = saveDataName;
    data.saveFolder     = currentDir;   % If you don't want to save in Current Folder, enter a folder name manually

    z_param = 2;
    
    clear rTemp;
    
    rTemp               = bfopen(D1{i});
    reader              = bfGetReader(D1{i});

    omeMetadata         = reader.getMetadataStore;
    voxelSizeX          = omeMetadata.getPixelsPhysicalSizeX(0).value;
    voxelSizeXdouble    = voxelSizeX.doubleValue();
    voxelSizeZ          = omeMetadata.getPixelsPhysicalSizeZ(0).value;
    voxelSizeZdouble    = voxelSizeZ.doubleValue(); 
    data.xy_size        = voxelSizeXdouble;
    data.z_size         = voxelSizeZdouble;

    zSize   = size(rTemp{1,1},1);

    % Create Z-Stack with Interleaved Channels
    clear temp;
    clear Icube;

    Icube = zeros(size(rTemp{1,1}{1,1},1),size(rTemp{1,1}{1,1},2),zSize);
    
    for z = 1:zSize
        temp             = rTemp{1,1}{z,1};
        Icube(:,:,z)     = double(temp);
        if (z == zSize)
            metadata = rTemp{1,2};
        end
    end
    z_num = size(Icube,3);

    data.metadata = metadata;

    I1 = double(Icube(:,:,1:2:z_num));
    I2 = double(Icube(:,:,2:2:z_num));
    
    clear I1norm;
    clear I2norm;
    
    I1norm = I1/max(I1(:));
    I2norm = I2/max(I2(:));
    
    clear cube;
    
    cube{1} = I1norm;
    cube{2} = I2norm;
    
    save(fullfile(saveFolder,sprintf('%s.mat',saveDataName)), 'cube','data'); 
end