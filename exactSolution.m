%exactSolution.m
%Purpose: finds the exact solution for Sod's shock tube 
%Parameters: (passed through the par structure)
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
%            t           the time at which the function evaluates
%                        the exact solution
%Outputs: density  a vector containing the exact solution for the gas
%                  density in every cell (same size as par.cellCoords param)
%         velocity a vector containing the exact solution for the gas
%                  velocity in every cell (same size as par.cellCoords param) 
%         pressure a vector containing the exact solution for the gas
%                  pressure in every cell (same size as par.cellCoords param)
function [density,velocity,pressure,cellCoords] = exactSolution(par)
cellCoords = par.cellCoords;
%create empty vectors for the density, velocity,
%and pressure outputs
density = zeros([1 size(par.cellCoords,2)]);
velocity = zeros([1 size(par.cellCoords,2)]);
pressure = zeros([1 size(par.cellCoords,2)]);


%use the initial conditions to determine the density, 
%velocity, and pressure in each of the five regions
%described in the assignment document

%Region 1
r1d = par.densL;
r1v = par.vxL;
r1p = par.presL;

%density, velocity, and pressure equal to 
%the initial left chamber density, velocity, and pressure. 

%also calculate the sound speed in Region 1

%Skip Region 2 for now (depends on x-coordinates)
%Regions 3 and 4
A = 2/(par.densR*(par.gamma+1));
B = par.presR*(par.gamma-1)/(par.gamma+1);
tol = 0.01;
startVal = (par.presL+par.presR)/2;
pStar = pStarSolver(startVal,tol,par);

r3p = pStar;
r3d = r1d*(r3p/r1p)^(1/par.gamma);
vxStar = par.vxR + (pStar - par.presR)*(A/(pStar+B))^(1/2);
r3v = vxStar;

r4p = pStar;
r4d = par.densR*(((par.presR*(par.gamma-1)+pStar*(par.gamma+1))/(pStar*(par.gamma-1)+par.presR*(par.gamma+1))));
r4v = vxStar;

%use the Newton's method solver to find p-star
%(with a tolerance of 0.01 and an initial p-star guess of (p_L+p_R)/2)
%(Note: the value of p-star has to be somewhere between p_L and p_R,
%       so the average (p_L + p_R)/2 is a good place to start 
%       looking for a solution.)

%use p-star to find the Region 3 density and the Region 4 density

%use p-star to find vx-star, the velocity in Regions 3 and 4

%Region 5
r5d = par.densR;
r5v = par.vxR;
r5p = par.presR;

%Set density, velocity, and pressure equal to the
%initial right chamber density, velocity, and pressure. 

%Next we need to find the start and end of each region.
%Note: Before adding data to the density, velocity, and pressure vectors,
%      use assert statements to ensure that the start and end of the
%      region you're working on both fall within the domain of the 
%      shock tube system. Including these assert statements will help
%      you answer Question 1 of this assignment. 

%ensure that the shock tube grid has at least 2 cells
assert(size(par.cellCoords, 2) > 1);
%calculate the width of each cell (could also just use par.dx parameter)
dx = par.cellCoords(1,2) - par.cellCoords(1,1);

%the size of Region 1 depends on the speed of 
%the head of the rarefaction wave.
csOne = ((par.gamma*par.presL)/r1d)^(1/2);
vHead = r1v - csOne;
regOneR = par.x0 + par.t*vHead;

%Copy the Region 1 density, velocity, and pressure
%into the appropriate cells in the output vectors.
density(1:(floor(regOneR*length(density)))) = r1d;
velocity(1:(floor(regOneR*length(velocity)))) = r1v;
pressure(1:(floor(regOneR*length(pressure)))) = r1p;

%the size of Region 2 depends on the speed of 
%the tail of the rarefaction wave.
csThree = ((par.gamma*r3p)/r3d)^(1/2);
vTail = (r3v - csThree);
regTwoR = par.x0 + par.t*vTail;

%Use the x-positions provided in par.cellCoords to find the density, velocity,
%and pressure everywhere inside Region 2. (Note: Region 2 is the only
%region where the density, velocity, and pressure are not uniform.)
r2x = regOneR:par.dx:regTwoR;
r2d = zeros(1,length(r2x));
r2v = zeros(1,length(r2x));
r2p = zeros(1,length(r2x));

for i = 1:length(r2x)
    r2d(i) = r1d*((2/(par.gamma+1))+((par.gamma-1)/(csOne*(par.gamma+1)))*(r1v - (r2x(i) - par.x0)./par.t))^(2/(par.gamma-1));
    r2v(i) = (2/(par.gamma+1))*(csOne+(par.gamma-1)*r1v/2+(r2x(i)-par.x0)./par.t);
    r2p(i) = r1p*((2/(par.gamma+1))+((par.gamma-1)/(csOne*(par.gamma+1)))*(r1v - (r2x(i) - par.x0)./par.t))^((2*par.gamma)/(par.gamma-1));
end

density((floor(regOneR*length(density))):(floor(regTwoR*length(density)))) = r2d;
velocity((floor(regOneR*length(velocity))):(floor(regTwoR*length(velocity)))) = r2v;
pressure((floor(regOneR*length(pressure))):(floor(regTwoR*length(pressure)))) = r2p;

%the size of Region 3 depends on the speed of the right-travelling
%contact discontinuity.
vCont = vxStar;
regThreeR = par.x0 + par.t*vCont;

%Copy the Region 3 density, velocity, and pressure
%into the output vectors.
density((floor(regTwoR*length(density))):(floor(regThreeR*length(density)))) = r3d;
velocity((floor(regTwoR*length(velocity))):(floor(regThreeR*length(velocity)))) = r3v;
pressure((floor(regTwoR*length(pressure))):(floor(regThreeR*length(pressure)))) = r3p;

%the size of Region 4 depends on the speed of the shock wave.
csFive = ((par.gamma*r5p)/r5d)^(1/2);
vShock = r5v + csFive*(((((par.gamma+1)*r4p)/(2*par.gamma*r5p))+((par.gamma-1)/(2*par.gamma)))^(1/2));
regFourR = (par.x0 + par.t*vShock);

%Copy the Region 4 density, velocity, and pressure
%into the output vectors.

density((floor(regThreeR*length(density))):(floor(regFourR*length(density)))) = r4d;
velocity((floor(regThreeR*length(velocity))):(floor(regFourR*length(velocity)))) = r4v;
pressure((floor(regThreeR*length(pressure))):(floor(regFourR*length(pressure)))) = r4p;
%Region 5 occupies the rest of the shock tube grid

%Copy the Region 5 density, velocity, and pressure
%into the remaining cells.
density((floor(regFourR*length(density))):par.xMax*length(density)) = r5d;
velocity((floor(regFourR*length(velocity))):par.xMax*length(density)) = r5v;
pressure((floor(regFourR*length(pressure))):par.xMax*length(density)) = r5p;


end

