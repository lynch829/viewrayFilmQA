function [beamCenters, beamSamples] = circProfiles(circX,circY,filmNorm,numBeams)
%% circProfiles
% [beamCenters, beamPoints] = circProfiles(circX,circY,filmNorm,numBeams)
%
% Arguements --
% circX: X coordinates of circle lying within the image to be analyzed.
% circY: Y coordinates
% filmNorm: Image from which pixel values are sampled along the circle
% given by circX and circY
% numBeams: Number of beams in the star shot image.
%
% Return Values --
% beamCenters: (numBeams *2) x 2 matrix containing the X and Y coordinates 
% of points on the centeral axis of each beam.  1 point per beam leg.
% beamSamples: (numBeams * 4) x 2 matrix containing the X and Y coordinates 
% of the two equidistant points used in the determination of each beam
% center point.


%% Compute pixel intensity values along circle given by circX and circY

pixelVals = zeros(length(circY),1);

% Bilinear interpolation of neighboring pixels
for ind = 1:length(pixelVals)
    pixelVals(ind) = interpImg(filmNorm,[circY(ind),circX(ind)]);
end


%% Filter pixel values

% Low-pass filter by convolving with a Gaussian
filterSigma = 4;
filterSize = 15;
x = linspace( -(filterSize/2), (filterSize/2), filterSize);
gaussFilter = exp( -x .^2 / (2 * filterSigma ^2));
gaussFilter = gaussFilter / sum(gaussFilter);

pixelVals = conv(pixelVals,gaussFilter,'same');

% Normalize pixel values so that peak detection is independent of input
% image intensity

pixelVals = mat2gray(pixelVals);

%% Identify local maxima of the pixel values (points on beam centerline) 

% Set starting value for delta, the difference threshold used to determine
% what is and is not a local maxima.

delta = 0.075;
[maxTab, minTab] = peakdet(pixelVals,delta);

% Verify that correct number of peaks were detected.  If not, search for
% correct delta value.  Stop when correct number of peaks are detected.

stepSize = .001;

if length(maxTab) < (numBeams * 2)

    % Lower delta value
    while length(maxTab) < (numBeams * 2)
        delta = delta - stepSize;
        [maxTab, minTab] = peakdet(pixelVals,delta);
    end
    
elseif length(maxTab) > (numBeams * 2)
    
    % Raise delta value
     while length(maxTab) > (numBeams * 2)
        delta = delta + stepSize;
        [maxTab, minTab] = peakdet(pixelVals,delta);
     end
     
end
  
    


%% Determine beam centerline points

% For each peak, take two points on either side of it at equal intensity.
% Calculate the midpoint and use it as a point on the beam centerline.


beamWidth = mean(diff(maxTab(:,1)));
beamSpacing = round(beamWidth/2);

beamLegCenters = zeros(numBeams *2,2);
legPoints = zeros(numBeams*2,4);

% Percentage of peak value to sample points at
percPeakVal = 0.925;
    
debugIso = zeros(numBeams*2,2);   

for ind = 1: numBeams * 2;
    
    peakPos = maxTab(ind,1);
    peakVal = maxTab(ind,2);
    isointensityPoints = zeros(1,2);

    % First beamleg
    if ind == 1  
    startIndex = 1;
	stopIndex = peakPos + beamSpacing;
    
    % Last beam leg
	elseif ind == (numBeams * 2)
    startIndex = peakPos - beamSpacing;
	stopIndex = length(pixelVals);	
    
    % All others
    else  
    startIndex = peakPos - beamSpacing;
	stopIndex = peakPos + beamSpacing;
    end
    
    % Locate two points with percPealVal % of peak intensity on either side    
    [~,isointensityPoints(1)] = min(abs(pixelVals(startIndex:peakPos) - (percPeakVal * peakVal)));
	[~,isointensityPoints(2)] = min(abs(pixelVals(peakPos:stopIndex) - (percPeakVal * peakVal)));
    
    
    isointensityPoints(1) = isointensityPoints(1) + startIndex - 1;
    isointensityPoints(2) = isointensityPoints(2) + peakPos - 1;
	
    debugIso(ind,1) = isointensityPoints(1);
    debugIso(ind,2) = isointensityPoints(2);
    legPoints(ind,:) = [circX(isointensityPoints(1)), circY(isointensityPoints(1)), circX(isointensityPoints(2)), circY(isointensityPoints(2))];
    
    % Take the point along the circle midway between the two isointensity
    % points as a point on the beam centerline
    centerIndex = round(mean(isointensityPoints));
    beamLegCenters(ind,:) = [circX(centerIndex), circY(centerIndex)];

    
    
end


% Group points on corresponding beam legs
beamCenters = zeros(numBeams,2,2);
beamSamples = zeros(numBeams,4,2);

for ind = 1:numBeams
    
    beamCenters(ind,:,1) = beamLegCenters(ind,:);
    beamCenters(ind,:,2) = beamLegCenters(ind + numBeams,:);
    
    beamSamples(ind,:,1) = legPoints(ind,:);
    beamSamples(ind,:,2) = legPoints(ind + numBeams,:);
end
    
    
end