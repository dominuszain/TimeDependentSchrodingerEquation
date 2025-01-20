// code of zain ul abideen
clear; clc;

// constants
hcut = 1
m = 1

// size of the grid matrices
lent = 100
lenx = 30

// see what works best, keeping in mind the convergence condition for the finite-difference method
delt = 0.01

// The Matrices.
R = zeros(lent, lenx)
I = zeros(lent, lenx)

//Gaussian initial wave function
function [y] = g(x)
    y = -10 * exp(-((x - 5) .^ 2) ./ (2 * 2 .^ 2));
endfunction

// tailor-made for this application
x = linspace(0, 10, lenx)
v = linspace(0, 10, lenx)

// defining important values
delx = x(2)-x(1)
V = g(v)

// initial conditions.
function [y]=f(x)
    y = exp(- (x - 4) .^ 2)
endfunction

R(1, :) = f(x)
I(1, :) = 0

// boundary condition. Might need to impliment some complicated form in the future to get more stable solutions
R(:, 1) = 0; R(:, lenx) = 0;
I(:, 1) = 0; I(:, lenx) = 0;


// defining a figure and giving it id
scf(0)

// solving the PDE.
for j = 1 : lent
    for i = 2 : lenx - 1
        I(j + 1, i) = I(j, i) - (V(i) .* R(j, i) .* delt ./ hcut) + (hcut .* delt ./ (2 .* m .* delx .* delx)) .* (R(j, i - 1) - 2 .* R(j, i) + R(j, i + 1))
        R(j + 1, i) = R(j, i) + (V(i) .* I(j, i) .* delt ./ hcut) - (hcut .* delt ./ (2 .* m .* delx .* delx)) .* (I(j, i - 1) - 2 .* I(j, i) + I(j, i + 1))
    end
    
    phy = sum(R(j, :) .* R(j, :) + I(j, :) .* I(j, :))
    
    R(j, :) = R(j, :) ./ sqrt(phy)
    I(j, :) = I(j, :) ./ sqrt(phy)
    
    arr(j, :) = R(j, :) .* R(j, :) + I(j, :) .* I(j, :)
    
    phy = sum(R(j, :) .* R(j, :) + I(j, :) .* I(j, :))
    
    // clearing the figure
    clf()
    
    // getting info about the current axes
    axes = gca()
    
    // setting the adjustable axes of the plot to constant values
    axes.data_bounds = [0, 0; 10, 0.3]

    plot(x, arr(j, :), "-b", "LineWidth", 2)
    plot(x, (g(x)+10) ./ 20, "-r", "LineWidth", 1)
    xtitle("Time Evolution of Probability Density", "x", "ψ*(x)ψ(x)")
    xgrid()
    
    // defines a naming scheme for plots upto 10,000
    if j < 10 then
        k = string("0000")+string(j)
    elseif j >= 10 & j < 100 then
        k = string("000")+string(j)
    elseif j>= 100 & j < 1000 then
        k = string("00")+string(j)
    elseif j>=1000 & j<10000 then
        k = string("0")+string(j)
    else k = string(j)
    end
    
    xs2jpg(0, k+".jpg")
    
     
end

// see the relevant paper for theoretical formulation.
