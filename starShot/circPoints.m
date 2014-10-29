function [circX, circY] = circPoints(theta, radius, approxCenter)
circX = approxCenter(1) + radius * cosd(theta);
circY = approxCenter(2) + radius * sind(theta);
end