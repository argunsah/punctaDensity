clear all
close all
clc

% Part3: Reads _info.mat Files and presents it for manual tracing

currentDir = pwd;
addpath(genpath(currentDir));

% Select the data to be analyzed
D1          = uipickfiles; % Select previously extracted Files (_info.mat Files)

% Change the save folder accordingly
saveFolder  = currentDir;
% or write the full folder name
% i.e. : saveFolder  = '/Users/Argunsah/Documents/punctaDensityResults';

for i = 1:size(D1,2)

   infoNum  = struct; 

   load(D1{i}); 
   smallTomatoBWF = puncta;
   
   imTemp = imadjust(adapthisteq(max(smallTomatoBWF,[],3)),[0 1]);

   figure, imshow(imcomplement(imTemp)); title(saveDataName,'Interpreter','none');
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

                    branchLen{i,xx}         = len3D;
                    numDendrites(i)         = xx;
                    shortestLinesAll{i,xx}  = meanSL;

                    b = xx;
            end         
       end
       catch ME
           ME
           delete(hh);
       end
   end
   
   saveas(gcf,fullfile(saveFolder,sprintf('%s_figure.png',saveDataName)),'png');
   
   infoNum.individualSizes  = individualSizes;
   infoNum.branchLen        = branchLen;
   infoNum.numDendrites     = numDendrites(i);
   infoNum.shortestLinesAll = shortestLinesAll;
   infoNum.fileName         = D1{i};

   save(fullfile(saveFolder,sprintf('%s_infoNums_%s.mat',saveDataName,date)), '-struct', 'infoNum');
   clear infoNum;
   clear individualSizes;
   clear branchLen;
   clear shortestLinesAll;
end