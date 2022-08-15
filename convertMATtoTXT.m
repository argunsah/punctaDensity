% Ali Ozgur Argunsah, Zurich, 2022

% Part4: Conversion on InfoNums Files to New Format for Analysis in R (txt - tab separated)

currentDir = pwd;
addpath(genpath(currentDir));

D1 = uipickfiles; % Select Files in _infoNums_ Folder

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
    writematrix(M,fullfile(path,sprintf('%s_newFormat.txt',name)),'Delimiter','tab');

    Mcell{i} = M;
end
