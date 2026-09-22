clear; clc; close all;
dirPath = uigetdir;
dapiList = dir([dirPath, '\', '*_cell.tiff']);
rfpList = dir([dirPath, '\', '*_RFP.tiff']);
gfpList = dir([dirPath, '\', '*_GFP.tiff']);
eRNAList = dir([dirPath, '\', '*_eRNA.tiff']);
allData = [];
for n = 1 : length(dapiList)
  fileToRead = [dirPath, '\', dapiList(n).name];
  rawData = importdata(fileToRead);

  cellNumber = size(unique(rawData), 1) - 1;
  fileData = zeros(cellNumber, 3); % RFP, GFP, eRNA
  if cellNumber >= 1
    cellTag = unique(rawData);
    cellTag = cellTag(2:end);

    rfpData = importdata([dirPath, '\', rfpList(n).name]);
    rfpNumber = size(unique(rfpData), 1) - 1;

    gfpData = importdata([dirPath, '\', gfpList(n).name]);
    gfpNumber = size(unique(gfpData), 1) - 1;

    eRNAData = importdata([dirPath, '\', eRNAList(n).name]);
    eRNANumber = size(unique(eRNAData), 1) - 1;
    
    for m = 1 : cellNumber
      [x, y] = find(rawData == cellTag(m));
      xCell = x(boundary(x,y)); yCell = y(boundary(x,y));

      % RFP
      if rfpNumber > 0
        [xRFP, yRFP] = find(rfpData ~= 0);
        in = inpolygon(xRFP, yRFP, xCell, yCell);
        if sum(in) > 0
          fileData(m, 1) = 1; % RFP positive
        end

%           in = inpolygon(xRFP, yRFP, x(boundary(x,y)), y(boundary(x,y)));
% 
%           plot(x(boundary(x,y)), y(boundary(x,y)),'LineWidth',2) % polygon
%           axis equal
% 
%           hold on
%           plot(xRFP(in), yRFP(in),'r+') % points inside
%           plot(xRFP(~in), yRFP(~in),'bo') % points outside
        
      end

      % GFP
      if gfpNumber > 0
        [xGFP, yGFP] = find(gfpData ~= 0);
        in = inpolygon(xGFP, yGFP, xCell, yCell);
        if sum(in) > 0
          fileData(m, 2) = 1; % GFP positive
        end
      end

      % eRNA
      if eRNANumber > 0
        for k = 1 : eRNANumber
          [xRNA, yRNA] = find(eRNAData ~= 0);
          in = inpolygon(xRNA, yRNA, xCell, yCell);
          if sum(in) > 0
            fileData(m, 3) = fileData(m, 3) + 1; % eRNA detected
            % calculate distance to GFP/RFP if they are positive?
%             if fileData(m, 1) == 1 || fileData(m, 2) == 1
%               1
%             end
          end
        end
      end
    end
    
  end

  allData = [allData; fileData];
end

boxData = strings(size(allData, 1), 1);
for n = 1 : size(allData, 1)
  if sum(allData(n, 1:2)) == 2
    boxData(n) = "GFP and RFP Positive";
  elseif allData(n, 1) == 1
    boxData(n) = "RFP Positive";
  elseif allData(n, 2) == 1
    boxData(n) = "GFP Positive";
  else
    boxData(n) = "None";
  end
end
boxplot(allData(:, 3), boxData, 'Symbol', '');
xLoc = zeros(size(allData, 1), 1);
hold on;
for n = 1 : size(allData, 1)
  if boxData(n) == "None"
    xLoc(n) = 1;
  elseif boxData(n) == "GFP Positive"
    xLoc(n) = 2;
  elseif boxData(n) == "RFP Positive"
    xLoc(n) = 3;
  elseif boxData(n) == "GFP and RFP Positive"
    xLoc(n) = 4;
  end
%   plot(xLoc(n), allData(n, 3), 'mo')
end
swarmchart(xLoc, allData(:,3));
ylabel("eRNA numbers")

% plotData = [allData, boxData];
% save([dirPath, '\', 'plotData.mat'], 'allData')