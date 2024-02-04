// code of zain ul abideen
clear; clc;

// constants
hcut = 1
m = 1

// size of the grid matrices.
lent = 20000
lenx = 30

// see what works best.
delt = 0.001

// The Matrices.
R = zeros(lent, lenx)
I = zeros(lent, lenx)

// defining a potential function.
function [y]=g(x)
    y = (x - 4.5) .^ 2
endfunction

// tailor-made for this application.
x = linspace(0, 10, lenx)
v = linspace(0, 10, lenx)

delx = x(2)
V = g(v)

// initial conditions.
function [y]=f(x)
    y = exp(- (x - 5) .^ 2)
endfunction

R(1, :) = f(x)
I(1, :) = 0

// boundary condition.
R(:, 1) = 0; R(:, lenx) = 0;
I(:, 1) = 0; I(:, lenx) = 0;

// solving the PDE.
for j = 2 : lent
    for i = 2 : lenx - 1
        I(j, i) = I(j - 1, i) - (V(i) .* R(j - 1, i) .* delt ./ hcut) + (hcut .* delt ./ (2 .* m .* delx .* delx)) .* (R(j - 1, i - 1) - 2 .* R(j - 1, i) + R(j - 1, i + 1))
        R(j, i) = R(j - 1, i) + (V(i) .* I(j - 1, i) .* delt ./ hcut) - (hcut .* delt ./ (2 .* m .* delx .* delx)) .* (I(j - 1, i - 1) - 2 .* I(j - 1, i) + I(j - 1, i + 1))
    end
    
    phy = sum(R(j, :) .* R(j, :) + I(j, :) .* I(j, :))
    
    R(j, :) = R(j, :) ./ sqrt(phy)
    I(j, :) = I(j, :) ./ sqrt(phy)
    
    arr(j, :) = R(j, :) .* R(j, :) + I(j, :) .* I(j, :)
    
    disp(phy)
    plot(x, arr(j,:))
    
end
//
//for k = 1 : 10 : lent
//    plot(x, arr(k, :), '-b')     
//end

// see the relevant paper for theoretical formulation.
