%shockTube.m
%main program for Assignment 2 shock tube simulation

function shockTube
clc
%Retrieve the shock tube parameters from the setup function.
par = setup();

%Use the functions exactSolution.m and roeSolution.m to find, respectively,
%the exact and simulated results for the state of the shock tube at the
%specified time.

[exactDensity, exactVelocity, exactPressure,xCoords] = exactSolution(par);
size(exactDensity)
[simDensity, simVelocity, simPressure,xCoordstest] = roeSolution(par);
size(simDensity)

%Use plots to compare the exact and simulated results.
%plot(xCoords,exactDensity);
%plot(xCoords,exactVelocity);
plot(xCoords,exactPressure);
hold on;
%plot(xCoordstest,simDensity);
%plot(xCoordstest,simVelocity);
plot(xCoordstest,simPressure);
hold off;
%title('Sod''s Shock Tube at t = 0.2 (Exact Solution)')
title('Sod''s Shock Tube at t = 0.2 (Exact vs. Roe)')
%title('Sod''s Shock Tube at t = 0.2 (Roe)')
xlabel('x')
%ylabel('Density')
%ylabel('Velocity')
ylabel('Pressure')
legend('Exact','Roe')
end

function par = setup()
    %Set the initial shock tube state here.
    par.densL = 1.0; %left chamber density
    par.vxL  = 0.75; %left chamber velocity
    par.presL = 1.0; %left chamber pressure

    par.densR = 0.125; %right chamber density
    par.vxR  = 0.0; %right chamber velocity
    par.presR = 0.1; %right chamber pressure

    %Set gas constant.
    par.gamma = 1.4; %specific heat ratio
    
    %It is also convenient to calcuate the initial
    %sound speed in the left chamber.
    par.csL = sqrt(par.gamma*par.presL/par.densL);

    %Set simulation parameters.
    par.cfl = 0.5; %Courant-Friedrichs-Lewy number

    par.maxCycles = 1000; %maximum number of cycles/iterations

    par.numXCells = 400; %number of cells in the simulation grid
    par.numXCellstest = 400; %test variable
    par.xMin = 0.0; %min x-position 
    par.xMax = 1.0; %max x-position
    par.dx = (par.xMax-par.xMin)/(par.numXCells-1); %cell width
    par.dxtest = (par.xMax-par.xMin)/(par.numXCellstest-1); %test variable
    %Create a vector that containes the x-positions of every 
    %cell in the simulation grid.
    par.cellCoords = par.xMin:par.dx:par.xMax;
    par.cellCoordstest = par.xMin:par.dxtest:par.xMax;

    %Set the diaphragm position.
    par.x0 = 0.5;

    %Choose the time at which you want to evaluate the 
    %shock tube solution. 
    par.t = 0.2;
end