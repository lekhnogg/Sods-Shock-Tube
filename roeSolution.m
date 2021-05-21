%roeSolution.m
%Purpose: uses Godunov's scheme and Roe's algorithm to numerically   
%         simulate the time evolution of Sod's shock tube 
%Parameters: (passed through par structure)
%            densL       initial density in the left chamber
%            vxL         initial velocity in the left chamber
%            presL       initial pressure in the left chamber
%            densR       initial density in the right chamber
%            vxR         initial velocity in the right chamber
%            presR       initial pressure in the right chamber
%            x0          the diaphragm position
%            cellCoords  the x-positions of the cells in the 
%                        shock tube grid
%            gamma       the specific heat ration for the ideal gas
%            t           the time at which the function will terminate
%                        the simulation and return the new state of
%                        the system
%            maxCycles   the maximum number of times that the simulation 
%                        will apply Godunov's schem (prevents infinite
%                        loops if the timesteps become too small to
%                        reach tLim)
%            cfl         the Courant-Friedrichs-Lewy number for 
%                        determining timestep size
%Outputs: density  a vector containing the final result for the gas
%                  density in every cell (same size as cellCoords param)
%         velocity a vector containing the final result for the gas
%                  velocity in every cell (same size as cellCoords param) 
%         pressure a vector containing the final result for the gas
%                  pressure in every cell (same size as cellCoords param)
%
function [density,velocity,pressure,cellCoords] = roeSolution(par)
% Establish value for epsilon
eps = 0.5;
cellCoords = par.cellCoordstest;
% Create vectors to store the fluxes through each side of the cells
fluxL = zeros(3,par.numXCellstest+1);
fluxR = zeros(3,par.numXCellstest+1);
% Before starting the simulation, we need to create a grid/vector that
% matches the initial state of the shock tube system. 

% Use the list of x-positions to determine the size
% of the simulation grid
numXCellstest = size(par.cellCoordstest, 2);
% Only run the simulation if there are at least two cells
assert(numXCellstest > 1);
% Calculate the cell width (could also just use par.dx parameter)
dx = par.cellCoordstest(1,2) - par.cellCoordstest(1,1)

% Determine the index for the boundary between
% the left and right chambers
diaphragmIndex = ceil((par.x0-par.cellCoordstest(1,1))/dx)

% Create the density vector
leftDens = par.densL*ones([1 diaphragmIndex]);
rightDens = par.densR*ones([1 numXCellstest-diaphragmIndex]);
currDens = [leftDens rightDens];

% Create the velocity vector
leftV = par.vxL*ones([1 diaphragmIndex]);
rightV = par.vxR*ones([1 numXCellstest-diaphragmIndex]);
currV = [leftV rightV];


% Create the pressure vector
leftPres = par.presL*ones([1 diaphragmIndex]);
rightPres = par.presR*ones([1 numXCellstest-diaphragmIndex]);
currPres = [leftPres rightPres];

% Create a 3xnumXCellstest matrix that stores the conservative variables
% for each x-position (density, momentum, and energy, in that order)
U = zeros(3,numXCellstest);
U(1,:) = currDens;
U(2,:) = currV.*currDens;
U(3,:) = (1/2).*(currDens).*(currV.^2) + (currPres./(par.gamma-1));
% Start the simulation at t=0
t = 0.0;

% Keep track of how many times we have updated the simulation grid
cycle = 0;

% Run the simulation until we reach the desired time 
% or the maximum number of cycles
while t < par.t && cycle <= par.maxCycles
    
    % We first add a ghost cell to each side of the state vectors.
    % Since we use transparent boundary conditions, we simply have to 
    % copy the state variables in the left and rightmost cells.
    
    % Create density vector with ghost cells
    size(currDens);
    currDens = [currDens(1,1) currDens currDens(1,size(currDens,2))];
    % Copy the densities in every i-th cell 
    %(i.e., the left cell densities)
    dens_i = currDens(1,1:numXCellstest+1);
    % Copy the densities in every i+1-th cell 
    % (i.e., the right cell densities)
    dens_i1 = currDens(1,2:numXCellstest+2);
    
    % Create velocity vector with ghost cells
    currV = [currV(1,1) currV currV(1,size(currV,2))];
    % Copy the velocities in every i-th cell 
    % (i.e., the left cell velocities)
    vel_i = currV(1,1:numXCellstest+1);
    % Copy the velocities in every i+1-th cell 
    % (i.e., the right cell velocities)
    vel_i1 = currV(1,2:numXCellstest+2);
    
    % Create pressure vector with ghost cells
    currPres = [currPres(1,1) currPres currPres(1,size(currPres,2))];
    % Copy the pressure in every i-th cell 
    % (i.e., the left cell pressures)
    pres_i = currPres(1,1:numXCellstest+1);
    % Copy the pressure in every i+1-th cell 
    % (i.e., the right cell pressures)
    pres_i1 = currPres(1,2:numXCellstest+2);
    
    % Use density, velocity, and pressure vectors 
    % to determine the energy in every cell
    
    % Copy the energies in every i-th cell 
    % (i.e., the left cell energies)
    ene_i = (1/2).*dens_i.*(vel_i).^2 + pres_i./(par.gamma-1);
    
    % Copy the energies in every i+1-th cell 
    % (i.e., the right cell energies)
    ene_i1 =(1/2).*dens_i1.*(vel_i1).^2 + pres_i1./(par.gamma-1);
    
    % Use the density, pressure, and energy vectors 
    % to determine the enthalpy in every cell
    
    % Copy the enthalpies in every i-th cell 
    % (i.e., the left cell enthalpies)
    enth_i = (ene_i + pres_i)./dens_i;
    % Copy the enthalpies in every i+1-th cell 
    % (i.e., the right cell enthalpies)
    enth_i1 = (ene_i1 + pres_i1)./dens_i1;
    
    % Compute the roe-averaged velocity across every cell boundary
    roeV = (((dens_i).^(1/2)).*vel_i+((dens_i1).^(1/2)).*vel_i1)./((dens_i)+(dens_i1));
    % Compute the roe-averaged enthalpy across every cell boundary
    roeEnth = (((dens_i).^(1/2)).*enth_i+((dens_i1).^(1/2)).*enth_i1)./((dens_i)+(dens_i1));
    % Compute the roe-averaged sound speed across every cell boundary
    roeCS = ((par.gamma-1).*(roeEnth-((1/2)*((roeV).^2)))).^(1/2);
    
    % Compute the eigenvalues for every inter-cell interface
    valOne = abs(roeV-roeCS);
    valTwo = abs(roeV);
    valThree = abs(roeV+roeCS);
    
    % Apply entropy fix (with an epsilon of 0.5)
    valOne(abs(valOne)<eps) = eps +((valOne(abs(valOne)<eps)).^2)./(2*eps); 
    valTwo(abs(valTwo)<eps) = eps +((valTwo(abs(valTwo)<eps)).^2)./(2*eps);
    valThree(abs(valThree)<eps) = eps +((valThree(abs(valThree)<eps)).^2)./(2*eps);

    % Create the eigenvectors for every inter-cell interface
    eigOne = [ones(1,length(roeCS)); (roeV-roeCS); (roeEnth - roeV.*roeCS)];
    eigTwo = [ones(1,length(roeCS)); roeV ;((1/2)*(roeV.^(2)))];
    eigThree = [ones(1,length(roeCS)); (roeV + roeCS) ;(roeEnth + roeV.*roeCS)];
    
    U1 = (dens_i1-dens_i);
    U2 = (dens_i1.*vel_i1)-(dens_i.*vel_i);
    U3 = (ene_i1-ene_i);

    % Compute the wave amplitudes
    ampTwo = ((par.gamma-1)./(((roeCS).^2))).*((U1.*(roeEnth-(roeV.^2)))+(roeV.*U2)-U3);
    ampOne = (1./(2.*roeCS)).*((U1.*(roeV+roeCS) - U2 - roeCS.*ampTwo));
    ampThree = U1 - (ampOne + ampTwo);
    
    % Take the sum of the eigenvectors multiplied by their 
    % respective eigenvalues and amplitudes
    sumOne = abs(eigOne).*ampOne.*valOne;
    sumTwo = abs(eigTwo).*ampTwo.*valTwo;
    sumThree = abs(eigThree).*ampThree.*valThree;
    totSum = sumOne+sumTwo+sumThree;
    
    % Create the flux vector F_L for every left cell
    % (should be a 3x(numXCellstest+1) matrix)
    fluxL(1,:) = dens_i.*vel_i;
    fluxL(2,:) = (dens_i.*(vel_i.^2))+pres_i;
    fluxL(3,:) = vel_i.*(ene_i+pres_i);

    % Create the flux vector F_R for every right cell
    % (should be a 3x(numXCellstest+1) matrix)
    fluxR(1,:) = dens_i1.*vel_i1;
    fluxR(2,:) = (dens_i1.*(vel_i1.^2))+pres_i1;
    fluxR(3,:) = vel_i1.*(ene_i1+pres_i1);
    
    % Calculate the flux through every inter-cell interface
    
    % Get the inter-cell fluxes for the left boundaries    
    intFluxL = (1/2).*(fluxL+fluxR)-(1/2).*totSum;
    % Get the inter-cell fluxes for the right boundaries 
    intFluxR = (1/2).*(fluxL+fluxR)-(1/2).*totSum;
    size(intFluxR);
    
    % Determine the max wavespeed
    maxWav = max([max(valOne(:)) max(valTwo(:)) max(valThree(:))]);
    
    % Calculate the timestep size
    tStep = par.cfl*(par.dxtest/maxWav);
    
    % Use Godunov's formula to update the conservative variable
    U(1,:) = U(1,:) +(tStep./par.dxtest).*(intFluxL(1,1:par.numXCellstest)-intFluxR(1,2:par.numXCellstest+1));
    U(2,:) = U(2,:) +(tStep./par.dxtest).*(intFluxL(2,1:par.numXCellstest)-intFluxR(2,2:par.numXCellstest+1));
    U(3,:) = U(3,:) +(tStep./par.dxtest).*(intFluxL(3,1:par.numXCellstest)-intFluxR(3,2:par.numXCellstest+1));
    size(currPres);
    % Derive the new primitive variable values from 
    % the updated conservative variable matrix
    currDens = U(1,:);
    currV = U(2,:)./U(1,:);
    currPres = (par.gamma-1).*((U(3,:)-((1/2).*(U(1,:).*currV.^2))));
    
    % Update the time
    t = t+tStep;
    
    % Update the cycle number
    cycle = cycle +1;

end

%After the simulation terminates, output the final 
%density, velocity, and pressure vectors
density = currDens;
velocity = currV;
pressure = currPres;

end

