%pStarSolver.m
%Purpose: uses Newton's method to iteratively compute the pressure in
%         Regions 3 and 4 of the exact shock tube solution 
%Parameters: startVal    initial guess for p-star
%            tol         the tolerance value for Newton's method
%                        (iterations will stop when the absolute value
%                         of the p-star function is less than or
%                         equal to tol) 
%            (The following parameters are passed 
%             through the par structure.)
%            densR       initial density in the right chamber
%            vxL         initial velocity in the left chamber
%            vxR         initial velocity in the right chamber
%            presL       initial pressure in the left chamber
%            presR       initial pressure in the right chamber
%            csL         initial sound speed in the left chamber
%            gamma       the specific heat ratio for the ideal gas
%            tol         the tolerance value for Newton's method
%                        (iterations will stop when the absolute value
%                         of the p-star function is less than or
%                         equal to tol) 
%Outputs: pStar the solution for p-star returned 
%               by Newton's method

function pStar = pStarSolver(startVal,tol,par);
f = startVal - (pStarFunc(startVal,par)/pStarDeriv(startVal,par));
fOld = 100;
while abs(fOld-f) >= tol
    fOld = f;
    f = f - pStarFunc(f,par)/pStarDeriv(f,par);
%use the provided functions (pStarFunc and pStarDeriv) to write 
%a Newton's method solver that will terminate when the value 
%returned by pStarFunc is sufficiently close to zero
%(i.e., when |pStarFunc| <= tol)
end
pStar = f;
end

%returns the value of the function that we need to minimize
function f = pStarFunc(pStar,par)
A = 2/(par.densR*(par.gamma+1));
B = par.presR*(par.gamma-1)/(par.gamma+1);
f = (pStar - par.presR)*sqrt(A/(pStar+B)) + (2*par.csL/(par.gamma-1))*((pStar/par.presL)^((par.gamma-1)/(2*par.gamma))-1) + par.vxR - par.vxL;
end

%returns the derivative of the function that we need to minimize
function fPrime = pStarDeriv(pStar,par)
A = 2/(par.densR*(par.gamma+1));
B = par.presR*(par.gamma-1)/(par.gamma+1);
fPrime = sqrt(A/(pStar+B)) - A*(pStar-par.presR)/(2*(pStar+B)^2*sqrt(A/(pStar+B))) + (par.csL/(par.gamma*par.presL))*(pStar/par.presL)^((par.gamma-1)/(2*par.gamma)-1);
end