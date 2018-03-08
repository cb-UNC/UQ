% This code used to solve the pump problem presented in the UQ notes. The
% permeability is set to be the expectation based on the lognormal
% statistics provided throughout.
%
% Rows go in the x-direction, columns in the y-direction.
% The model's indices are numbered from left to right, top to bottom such
% that, for node i in x and j in y, the center node is C*(j-1)+i. The
% indices which are used are labeled below
%
%       IC = C*(j-1)+i;     IC is the index for the center node
%       PY = IC + C;        PY is the index for +1 in the Y direction
%       PX = IC + 1;        PX is the index for +1 in the X direction
%       MY = IC - C;        MY is the index for -1 in the Y direction
%       MX = IC - 1;        MX is the index for -1 in the X direction
%
% This indexing scheme can be visualized as below
%                   
%                 (i-C)
%                   |
%                   |
%       (i-1) <---- i ----> (i+1)
%                   |
%                   |
%                 (i+C)
%
% 

% Darcy flow, steady-state, variable K, 2D
function DF_SS_VK_2D_PUMP_RAND
% System parameters
K_EXP = 0.664;                  % m/d
LENGTH = 1000;                  % m
HEIGHT = 1000;                  % m
WIDTH  = 50;                    % m  
Href = 10;                      % m
Hriv = 30;                      % m



% Model Parameters
R = 21;             % Nodes per Row
C = 21;             % Nodes per Column
N = R*C;            % Total number of nodes
dx = LENGTH/C;      % m
dy = HEIGHT/R;      % m  

% Compute these in advance for simplicity
dx2 = dx*dx;
dy2 = dy*dy;

% Compute the X- and Y-coordinates for each node
X = zeros(N,1);
Y = zeros(N,1);
hdx = dx/2;
hdy = dy/2;
for i=1:C
    for j=1:R
        IC = C*(j-1)+i;
        X(IC) = hdx + (i-1)*dx;
        Y(IC) = hdy + (j-1)*dy;
    end
end

% Pump Parameters
Q = 6000;                       % m3/d
T = K_EXP * WIDTH;              % m2/d
radE = 500;                     % m

% Pump/reference coordinates
i = 1; j = ceil(R/2); IC = C*(j-1)+i;
X_REF  = X(IC)-dx;   Y_REF  = Y(IC);
X_PUMP = X_REF-radE; Y_PUMP = Y_REF;

% Construct permeability field
K = zeros(N,1);
mu = -7.57; sigma = sqrt(0.8);
for i=1:C
    for j=1:R
        IC = C*(j-1)+i;
        K(IC) = lognrnd(mu,sigma);
    end
end


%%%%%%%%%%% Start the calculations to determine the head field %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Allocate for M and b
M = zeros(N,N);     % Model matrix of discretized equations
b = zeros(N,1);     % RHS of model

% Populate b
for j=1:R
    LI = C*(j-1)+1;
    RI = C*(j-1)+C;
    
    X_LI = X(LI)-dx; Y_LI = Y(LI);
    rad = ((X_PUMP-X_LI)^2 + (Y_PUMP-Y_LI)^2 )^(1/2);
    
    H_LI = (Q*log(rad/radE))/(2*pi*T) + Href;
    H_RI = Hriv;
    
    b(LI) = -(K(LI)/(dx2) )*H_LI;            % Left
    b(RI) = -(K(RI)/(dx2) )*H_RI;            % Right
end

% Populate M
for i=1:C
    for j=1:R
        % Prepare indices
        IC = C*(j-1)+i;
        PY = IC + C; PX = IC + 1;
        MY = IC - C; MX = IC - 1;
        
        % Calculate Effective Viscosities
        if (i>1), K_MX = K_Eff(K(IC),K(MX)); end     % Excludes Left 
        if (i<C), K_PX = K_Eff(K(IC),K(PX)); end     % Excludes Right 
        if (j>1), K_MY = K_Eff(K(IC),K(MY)); end     % Excludes Top 
        if (j<R), K_PY = K_Eff(K(IC),K(PY)); end     % Excludes Bottom

        % Populate M
        if (i>1), M(IC,MX) = K_MX/dx2; M(IC,IC) = M(IC,IC)-K_MX/dx2; end % Excludes Left 
        if (i<C), M(IC,PX) = K_PX/dx2; M(IC,IC) = M(IC,IC)-K_PX/dx2; end % Excludes Right 
        if (j>1), M(IC,MY) = K_MY/dy2; M(IC,IC) = M(IC,IC)-K_MY/dy2; end % Excludes Top
        if (j<R), M(IC,PY) = K_PY/dy2; M(IC,IC) = M(IC,IC)-K_PY/dy2; end % Excludes Bottom

        % Add Dirichlet left and right
        if (i==1), M(IC,IC) = M(IC,IC)-K(IC)/dx2; end % Front
        if (i==C), M(IC,IC) = M(IC,IC)-K(IC)/dx2; end % Back
    end
end

opts.SYM = true;
H = linsolve(M,b,opts);

% Visualize the data 
figure 
Visualize(H,R,C,dx,dy);

figure
VisualizeLatCent(H,R,C,Y);

%%%%%%%%%%%% End the calculation of the head field %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end

function K_Eff = K_Eff(K_1,K_2)
    K_Eff = 2/(1/K_1 + 1/K_2);
end


% Visualize Head Values
function Visualize(H,R,C,dx,dy)
% Create axes
X = zeros(C,1);
    for i=1:C , X(i) = i*dx; end
Y = zeros(R,1);
    for j=1:R , Y(j) = j*dy; end

% Create viewable H matrix
H_View = zeros(R,C);
for i=1:C
    for j=1:R
        IC = C*(j-1)+i;
        H_View(j,i) = H(IC);
    end
end

figure
pcolor(X,Y,H_View);
ylabel(colorbar,'Pressure Head (m)');

end

function H_Lat = VisualizeLatCent(H,R,C,Y)
H_Lat = zeros(R,1);
Y_Lat = zeros(R,1);

i = ceil(C/2); 
for j=1:R
    IC = C*(j-1)+i;
    H_Lat(j) = H(IC);
    Y_Lat(j) = Y(IC);
end

plot(Y_Lat,H_Lat);

end

