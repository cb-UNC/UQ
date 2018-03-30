% Processes the data from the pump sweep

function PUMP_SWEEP_VIZ(sweepH, sweepVX, sweepVY, expH, expVX, expVY, X, Y)
% Size the arrays
N = size(sweepH);
sweepN = N(2)                                                    %#ok<NOPRT>

% Plot some random fields
vQSubplots(sweepVX, sweepVY, X, Y, sweepN);

% Take averages of the data
meanH = mean(sweepH, 2);
meanVX = mean(sweepVX,2);
meanVY = mean(sweepVY,2);

% Compare the means to the expectation
[errorX, errorY] = aveCompVs(expVX, expVY, meanVX, meanVY, X, Y); 
errorV = (errorX + errorY)/2                                     %#ok<NOPRT>
vArray = [expVX, expVY];
relErrorV = errorV/norm(vArray)                                  %#ok<NOPRT,NASGU>

errorH = norm(expH - meanH)                                      %#ok<NOPRT>
relErrorH = errorH/norm(expH)                                    %#ok<NOPRT,NASGU>

end

function vQSubplots(sweepVX, sweepVY, X, Y, sweepN)
vX = sweepVX(:,1);
vY = sweepVY(:,1);

figure
subplot(2,2,1)
quiver(X,Y,vX,vY);

vX = sweepVX(:,sweepN);
vY = sweepVY(:,sweepN);

subplot(2,2,2)
quiver(X,Y,vX,vY);

vX = sweepVX(:,2);
vY = sweepVY(:,2);

subplot(2,2,3)
quiver(X,Y,vX,vY);

vX = sweepVX(:,3);
vY = sweepVY(:,3);

subplot(2,2,4)
quiver(X,Y,vX,vY);

end

function [errorX, errorY] = aveCompVs(expVX, expVY, meanVX, meanVY, X, Y)
% Plot averages next to one another
figure
subplot(1,2,1)
quiver(X,Y,expVX,expVY)
title('Expectation Velocity')

subplot(1,2,2)
quiver(X,Y,meanVX,meanVY)
title('Mean of Random Velocities')

errorX = norm(expVX-meanVX);
errorY = norm(expVY-meanVY);

end

