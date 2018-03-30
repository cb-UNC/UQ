function PUMP_SWEEP_PT_STAT(sweepVX, sweepVY, X, Y)
% Initial particle position (PX,PY)
PX = 500;        % m
PY = 500;        % m

% Change in time measure
dt = 1e3;          % days

% Numerical checks
check = 1;
maxIt = 50000;
it = 1;

for j=1:4
    % Velocity fields used
    VX = sweepVX(:,j);
    VY = sweepVY(:,j);
    while (check==1)
        [PX,PY] = P_Update(VX,VY,PX,PY,X,Y,dt);
        
        if (j==1), PXR_1(it) = PX; PYR_1(it) = PY; end
        if (j==2), PXR_2(it) = PX; PYR_2(it) = PY; end
        if (j==3), PXR_3(it) = PX; PYR_3(it) = PY; end
        if (j==4), PXR_4(it) = PX; PYR_4(it) = PY; end

        it = it+1;
        if (PX < 0)
            check = 0;
        end
        if (it > maxIt)
            check = 0;
        end
    end
    disp(it)
    
    % Reset variables
    check = 1;
    it = 1;
    PX = 500;
    PY = 500;
end

figure
plot(PXR_1,PYR_1,PXR_2,PYR_2,PXR_3,PYR_3,PXR_4,PYR_4)



end

function [PX,PY] = P_Update(VX,VY,PX,PY,X,Y,dt)
N = size(VX,1);
CD = 1000;      % Closest Distance (CD)
for i=1:N
    DP = ((PX-X(i))^2 + (PY-Y(i))^2 )^(1/2);    % Distance to Particle (DP)
    if (DP < CD) 
        IC = i;
        CD = DP;
    end
end

PX = PX + VX(IC)*dt;
PY = PY + VY(IC)*dt;
end