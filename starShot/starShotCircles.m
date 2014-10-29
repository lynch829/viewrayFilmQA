function [isocenterRadius, residualNorms] = starShotCircles(imageFileName, numBeams, dpi, invert)
%% [isocenterRadius, residualNorms] = starShotCircles(imageFileName, numBeams, dpi)
%
% Arguments:
%
% imageFileName: string containing path to image file to be analyzed
% numBeams: number of beams (not beam legs) on the star shot film
% dpi: Scanner resolution in dots per inch
% invert: Boolean.  Does the image need to be inverted?  If beams are dark and
% the background is bright, set invert to true.  Otherwise, set invert to
% false.
%
% Return Values:
%
% isocenterRadius: Radius of the minimum inscribed circle
% residualNorms: Norm of the residuals after fitting a centerline to each
% beam with polyfit.
%
%
% Dependencies:
% Image Processing Toolbox
% peakdet
% lineSegmentIntersect
% interpImg

%% Load scanned film image

filmImg = imread(imageFileName);

%% Invert, normalize, display and crop image

cropWindow = figure;

% Automatically enhance contrast for display only.  Pixel values used in the
% determination of beam centerlines are unchanged.

if (invert)
[filmNormR,filmRect] = imcrop(imadjust(1 - mat2gray(filmImg(:,:,1))));

else
[filmNormR,filmRect] = imcrop(imadjust(mat2gray(filmImg(:,:,1))));
end
    

filmImg = imcrop(filmImg,filmRect);

% % Check channel for highest response
% if (invert)
% filmNormR = 1 - mat2gray(filmImg(:,:,1));
% filmNormG = 1 - mat2gray(filmImg(:,:,2));
% filmNormB = 1 - mat2gray(filmImg(:,:,3));
% else
% filmNormR = mat2gray(filmImg(:,:,1));
% filmNormG = mat2gray(filmImg(:,:,2));
% filmNormB = mat2gray(filmImg(:,:,3));
% end
%     
% if max(filmNormR(:)) < max(filmNormG(:))   
%     filmNorm = filmNormG;
% else
%     filmNorm = filmNormR;
% end
% if max(filmNorm(:)) < max(filmNormB(:))
%     filmNorm = filmNormB;
% end

% Use red channel
filmNorm = filmNormR;
close(cropWindow) 

%% Display processed image

% Automatically enhance contrast for display only.  Pixel values used in the
% determination of beam centerlines are unchanged.
imDisp = figure;
imDispAxes=imagesc(imadjust(filmNorm)); colormap gray; axis image;
hold on;

% Calculate circles along which beam profiles are taken
title('Select approximate radiation isocenter.','FontSize',24);
% Get approximate isocenter
approxCenter = ginput(1);
title('');



% Prompt for start of sampling.
 title('Choose point.','FontSize',24)
 samplingStart = ginput(1);
 title('');
 

%% Grid search around approximate isocenter

isoSearchWidth = 10;
isoStepSize = 1;
isoX = [approxCenter(1) - (isoSearchWidth/2): isoStepSize: approxCenter(1) + (isoSearchWidth/2)];
isoY = [approxCenter(2) - (isoSearchWidth/2): isoStepSize: approxCenter(2) + (isoSearchWidth/2)];
[isoX,isoY] = meshgrid(isoX,isoY);
isoX = isoX(:);
isoY = isoY(:);



gridResults = struct([]);
bestResid = inf;
bestResidInd = 0;
isoSearch = waitbar(0,'Determining beam centerlines...');


for isoSearchInd = 1:length(isoX)
    
approxCenter = [isoX(isoSearchInd), isoY(isoSearchInd)];

% Determine limiting dimension of image and calculate maximum radius.  
edgeDistances = [fliplr(size(filmNorm)) - approxCenter; approxCenter];
maxRadius = .925 * min(edgeDistances(:));

% Set minimum radius
minRadius =  .80 * maxRadius;


% Calculate angle between the x axis and chosen point.  Start sampling from
% this angle.
 
 sampleRadius = norm(approxCenter - samplingStart);
 xOffset = samplingStart(1) - approxCenter(1);
 thetaStart = acosd(xOffset / sampleRadius);

 % Angle to start sampling pixel values
%thetaStart = 0;

% Get points along circles
theta = [thetaStart:0.1:(360 + thetaStart)];
radii = [minRadius : (maxRadius - minRadius) / 4 : maxRadius];
%radii = [ .9 * maxRadius, maxRadius];

% 3D Matrix to hold pixel coordinates of circles.  First column is X,
% second is Y. 3rd dimension is circle number in order of increading radius.

circles = zeros(length(theta),2,length(radii));

% Compute circle pixel coordinates
for ind = 1:length(radii)
	radius = radii(ind);
	[circles(:,1,ind), circles(:,2,ind)] = circPoints(theta,radius,approxCenter);
end	

% Prompt user for confirmation (?), if acceptable continue


% Calcultate the coordinates of the beam centers lying along each circle

allBeamCenters = [];
allBeamFWHM = [];
beamFits = zeros(numBeams,2);
residualNorms = zeros(1,numBeams);
clear beamCenters;

for ind = 1:length(radii)

    [beamCenters, beamSamples] = circProfiles(circles(:,1,ind),circles(:,2,ind),filmNorm,numBeams);
    allBeamCenters = cat(3,allBeamCenters,beamCenters);
    allBeamFWHM = cat(3,allBeamFWHM,beamSamples);
end


for ind = 1:numBeams
[beamFits(ind,:), S] = polyfit(allBeamCenters(ind,1,:), allBeamCenters(ind,2,:),1);
residualNorms(ind) = S.normr;
end


% Store results from this run
gridResults(isoSearchInd).center = approxCenter;
gridResults(isoSearchInd).beamFits = beamFits;
gridResults(isoSearchInd).residualNorms = residualNorms;
gridResults(isoSearchInd).allBeamCenters = allBeamCenters;
gridResults(isoSearchInd).allBeamFWHM = allBeamFWHM;
gridResults(isoSearchInd).circles = circles;

if max(residualNorms) < bestResid;
    bestResidInd = isoSearchInd;
    bestResid = max(residualNorms);
end

waitbar(isoSearchInd/length(isoX),isoSearch);

end
close(isoSearch);


%% Plot centerlines with best residuals
for ind = 1:length(radii)
plot(gridResults(bestResidInd).circles(:,1,ind),gridResults(bestResidInd).circles(:,2,ind),'r.');
end

% Plot beam centers
for ind = 1:numBeams 
plot(squeeze(gridResults(bestResidInd).allBeamCenters(ind,1,:)),squeeze(gridResults(bestResidInd).allBeamCenters(ind,2,:)),'o','MarkerSize',5,'MarkerFaceColor','b','Color','b');
plot(squeeze(gridResults(bestResidInd).allBeamFWHM(ind,1,:)),squeeze(gridResults(bestResidInd).allBeamFWHM(ind,2,:)),'o','MarkerSize',5,'MarkerFaceColor','g','Color','g');
plot(squeeze(gridResults(bestResidInd).allBeamFWHM(ind,3,:)),squeeze(gridResults(bestResidInd).allBeamFWHM(ind,4,:)),'o','MarkerSize',5,'MarkerFaceColor','g','Color','g');

end


% Plot beams
X = [1:size(filmNorm,2)];
for ind = 1:numBeams
	plot(X,polyval(gridResults(bestResidInd).beamFits(ind,:),X),'r','LineWidth',1);
end

%% Compute minimum intersecting circle

[radiationRadius, radiationCenter] = minCircle(gridResults(bestResidInd).beamFits,size(filmNorm));

% Plot circle
thetaDisp = [1:360];
plot(radiationCenter(1), radiationCenter(2),'bo','MarkerSize',5,'MarkerFaceColor','b')
plot(radiationCenter(1) + radiationRadius * cosd(thetaDisp),radiationCenter(2) + radiationRadius * sind(thetaDisp),'b')

% Convert radius in pixels to mm using input dpi
isocenterRadius = (radiationRadius / dpi) * 25.4;

% Round to three significant figures
isocenterRadius = round(isocenterRadius * 100) / 100;

% Get laser isocenter
laserCenter = laserIntersect(filmImg);
plot(laserCenter(1),laserCenter(2),'ko','MarkerSize',5,'MarkerFaceColor','k');

distToLaser = norm(laserCenter - radiationCenter);
mmDist = (distToLaser / dpi) * 25.4;
mmDist = round(mmDist * 100) / 100;
mmDistStr = num2str(mmDist);

% Display radius on plot
mmRadius = num2str(isocenterRadius);
xlabel(['R = ' mmRadius ' mm. ' mmDistStr ' mm from laser isocenter.' ],'FontSize',24)



    

end

