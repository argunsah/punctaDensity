% Microglia vs Neuron Interactions/Co-localization Analysis
% Ali Ozgur Argunsah and Lorenzo Gesuita, 2017, Zurich.
% Input: 2 channel RGB z-stack

%% This part is for automatic extraction First Part, to create .mat files
clear all
close all
clc

addpath(genpath('C:\Users\argunsah\Documents\lorenzoPaperCodes\punctaDensityCodes'));

% Change data folder accordingly
dataFolder  = 'W:\Assembly\Lorenzo\Pics_Confocal\Confocal_Microglia Project Images\Microglia_Iontoporations_mSyp_Cx3cl1';

% Select the data to be analyzed
D1          = uipickfiles('FilterSpec',dataFolder); % Select Files with .oif extensions

% Change the save folder accordingly
saveFolder  = 'E:\';

for i = 1:size(D1,2)
    close all;
    
    % Create Info Structure
    clear data;
    data  = struct;

    regexpInd       = regexp(D1{i},'\');
    saveDataName    = D1{i}(regexpInd(end)+1:end);
    saveDatanameInd = regexp(saveDataName,'\.');
    saveDataName    = saveDataName(1:saveDatanameInd(end)-1);
    
    saveDatanameEmptySpaces = regexp(saveDataName,' ');
    for tt = 1:length(saveDatanameEmptySpaces)
        saveDataName(saveDatanameEmptySpaces(tt)) = '_';
    end
    
    data.name           = D1{i};
    data.analysisday    = date;
    data.Analyzer       = 'Automatic';  % What is the name of the person who analyzed the data
    data.DataCollected  = 'Lorenzo';    % Who collected the data?
    data.saveDataName   = saveDataName;

    data_name_all       = fullfile(D1{i});
    temp                = regexp(data_name_all,'/');

    if isempty(temp)
        temp = regexp(data_name_all,'\');
    end

    data_name           = data_name_all(temp(end)+1:end);
    tempDot             = regexp(data_name,'\.');
    data.saveFolder     = saveFolder;

    z_param = 2;
    
    clear rTemp;
    
    rTemp   = bfopen(D1{i});
    zSize   = size(rTemp{1,1},1);

    % Create Z-Stack with Interleaved Channels
    clear temp;
    clear Icube;
    
    for z = 1:zSize
        temp             = rTemp{1,1}{z,1};
        Icube(:,:,z)     = double(temp);
        if (z == zSize)
            metadata = rTemp{1,2};
        end
    end
    z_num = size(Icube,3);

    I1      = double(Icube(:,:,1:2:z_num));
    I2      = double(Icube(:,:,2:2:z_num));
    
    clear I1norm;
    clear I2norm;
    
    I1norm   = I1/max(I1(:));
    I2norm   = I2/max(I2(:));
    
    clear cube;
    
    cube{1} = I1norm;
    cube{2} = I2norm;
    
    save(fullfile(saveFolder,sprintf('%s.mat',saveDataName)), 'cube'); 
end

%% This part is for automatic extraction Second Part, to create .info files
D1 = uipickfiles('FilterSpec',dataFolder); % Select Files in mat Folder

for i = 1:size(D1,2)
    i
    load(D1{i}); 
    puncta      = bil_wiener_filt(cube{1},3,1);
    tomato      = cube{2};
    
    thr         = graythresh(puncta(:));
    tomato_thr  = graythresh(tomato(:));
    bw_tomato   = tomato>tomato_thr;

    bw          = (puncta>thr);
    D           = bwdist(~bw);
    D           = -D;
    D(~bw)      = Inf;
    
    L           = watershed(D);
    L(~bw)      = 0;
    
    data.L          = L;
    data.puncta     = puncta;
    data.bw_tomato  = tomato;
    data.tomato     = tomato;
    
    pos  = regexp(D1{i},'\.');
    save(fullfile(sprintf('%s_info.mat',D1{i}(1:pos(end)-1))), '-struct', 'data');
end
%% This Part is for Manual Tracing 
D1 = uipickfiles('FilterSpec',dataFolder); % Select Files in _info Folder

xy_size = 0.198;
z_size  = 0.5;

for i = 1:size(D1,2)
   i
   infoNum  = struct;
   pos      = regexp(D1{i},'\.');
   
   load(D1{i}); 
   
   smallTomatoBWF = puncta;
   
   figure, imshow(imcomplement(imadjust(adapthisteq(max(smallTomatoBWF,[],3)),[0 1])));
   h = imcontrast(gca);
   waitfor(h);
   
   clear x;
   clear y;

   b      = 1000;
   xx     = 1;
   
   while xx < b
       clear x;
       clear y;
       clear z;
       
       [x,y] = ginputzoom;
       
       try
       if length(x)>1
       
           hA = gca;
           resetplotview(hA,'InitializeCurrentView');
           set(hA,'xlim',[1 1600]);
           set(hA,'ylim',[1 1600]);

           for k = 1:length(x)
               c    = max(squeeze(tomato(round(y(k)),round(x(k)),:)));
               vals = find(squeeze(tomato(round(y(k)),round(x(k)),:))==c);
               z(k) = mean(vals);
           end

           meanSL = [];
           for n = 1:length(x)-1
               SourcePoint    = [round(y(n))   round(x(n))   z(n)]';
               StartPoint     = [round(y(n+1)) round(x(n+1)) z(n+1)]';

               T1_FMM1        = msfm2d(double(max(smallTomatoBWF,[],3)+1e-8), SourcePoint(1:2), true, true);
               T1_FMM2        = msfm2d(double(max(smallTomatoBWF,[],3)+1e-8), StartPoint(1:2) , true, true);

               ShortestLine   = shortestpath(T1_FMM1,StartPoint(1:2),SourcePoint(1:2));
               ShortestLine2  = shortestpath(T1_FMM2,SourcePoint(1:2),StartPoint(1:2));

               meanSL         = [meanSL ;meanOfTwoLines(ShortestLine,ShortestLine2)];
           end
           hold on, 
           hh = scatter(meanSL(:,2),meanSL(:,1),1,'.r');
           drawnow;

            answer = questdlg('Do you like the Tracing?', ...
                'Tracing Quality', ...
                'Yes','No','Done','Yes');
            % Handle response
            switch answer
                case 'Yes'
                    D            = diff(meanSL*xy_size,1,1);
                    len          = double(trace(sqrt(D*D.')));
                    len3D        = sqrt(len.^2+((z(2)-z(1))*z_size).^2);
                    profileC     = improfile(max(double(tomato),[],3),meanSL(:,2),meanSL(:,1));
                    profileL     = improfile(max(double(L)>0,[],3),meanSL(:,2),meanSL(:,1));
                    diffProfileL = diff(profileL);

                    [aaa1,bbb1]  = find(diffProfileL==1);
                    [aaa2,bbb2]  = find(diffProfileL==-1);

                    if aaa1(1) > aaa2(1)
                        individualSizes{i,xx}(1) = aaa2(1)*xy_size;
                        for s = 2:length(aaa2)
                            individualSizes{i,xx}(s) = (aaa2(s) - aaa1(s-1))*xy_size;
                        end
                    else
                        for s = 1:length(aaa2)
                            individualSizes{i,xx}(s) = (aaa2(s) - aaa1(s))*xy_size;
                        end
                    end

                    branchLen{i,xx} = len3D;
                    numDendrites(i) = xx;
                    shortestLinesAll{i,xx} = meanSL;
                    xx = xx + 1;
                case 'No'
                    delete(hh);
                case 'Done'
                    D            = diff(meanSL*xy_size,1,1);
                    len          = double(trace(sqrt(D*D.')));
                    len3D        = sqrt(len.^2+((z(2)-z(1))*z_size).^2);
                    profileC     = improfile(max(double(tomato),[],3),meanSL(:,2),meanSL(:,1));
                    profileL     = improfile(max(double(L)>0,[],3),meanSL(:,2),meanSL(:,1));
                    diffProfileL = diff(profileL);

                    [aaa1,bbb1]  = find(diffProfileL==1);
                    [aaa2,bbb2]  = find(diffProfileL==-1);

                    if aaa1(1) > aaa2(1)
                        individualSizes{i,xx}(1) = aaa2(1)*xy_size;
                        for s = 2:length(aaa2)
                            individualSizes{i,xx}(s) = (aaa2(s) - aaa1(s-1))*xy_size;
                        end
                    else
                        for s = 1:length(aaa2)
                            individualSizes{i,xx}(s) = (aaa2(s) - aaa1(s))*xy_size;
                        end
                    end

                    branchLen{i,xx} = len3D;
                    numDendrites(i) = xx;
                    shortestLinesAll{i,xx} = meanSL;

                    b = xx;
            end         
       end
       catch ME
           ME
           delete(hh);
       end
   end
   
   saveas(gcf,fullfile(sprintf('%s_figure.png',D1{i}(1:pos(end)-1))),'png');
   
   infoNum.individualSizes  = individualSizes;
   infoNum.branchLen        = branchLen;
   infoNum.numDendrites     = numDendrites(i);
   infoNum.shortestLinesAll = shortestLinesAll;
   infoNum.fileName         = D1{i};

   save(fullfile(sprintf('%s_infoNums_20200901.mat',D1{i}(1:pos(end)-1))), '-struct', 'infoNum');
   clear infoNum;
   clear individualSizes;
   clear branchLen;
   clear shortestLinesAll;
end

%% Conversion on InfoNums Files to New Format for Analysis in R (txt - tab separated)
% each column has in row one the lenght of the branch; each puncta size
% is listed below.

D1 = uipickfiles('FilterSpec',dataFolder); % Select Files in _infoNums_ Folder

for i = 1:size(D1,2)
    load(D1{i});
    [path,name,ext] = fileparts(D1{i});
    newBranchLen    = cell2mat(branchLen);
    
    clear M;
    maxSize = 0;
    
    for j = 1:size(individualSizes,1)
        for k = 1:size(individualSizes(j,:),2)
            temp = individualSizes{j,k};
            if ~isempty(temp)
                if length(temp) > maxSize
                    maxSize = length(temp);
                end
            end
        end
        
        M = nan(maxSize+1,size(individualSizes(j,:),2));
        M(1,:) = newBranchLen;
        
        for k = 1:size(individualSizes(j,:),2)
            temp = individualSizes{j,k};
            if ~isempty(temp)
                M(2:length(temp)+1,k) = temp;
            end
        end
    end
   % figure, imagesc(M(2:end,:),[0 7]);
   % title(name);
   % drawnow;
    writematrix(M,fullfile(path,sprintf('%s_newFormat.txt',name)),'Delimiter','tab');
  
  Mcell{i} = M;
end
