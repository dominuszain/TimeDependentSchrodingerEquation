// code of zain ul abideen
clear; clc;

// constants
hcut = 1
m = 1

// size of the grid matrices.
lent = 5000
lenx = 51

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
end

// plotting all iterations.
//for k = 1 : 20 : lent
//    plot(x, R(k, :), '-b')
//    plot(x, I(k, :), ':r')      
//end

for k = 1 : lent
    arr(k, :) = R(k, :) .* R(k, :) + I(k, :) .* I(k, :)    
end

//for k = 1 : 10 : lent
//    plot(x, arr(k, :), '-b')     
//end

x1 = zeros(1, lent)
x2 = zeros(1, lent)
x3 = zeros(1, lent)
x4 = zeros(1, lent)
x5 = zeros(1, lent)
x6 = zeros(1, lent)
x7 = zeros(1, lent)
x8 = zeros(1, lent)
x9 = zeros(1, lent)
x10 = zeros(1, lent)
x11 = zeros(1, lent)
x12 = zeros(1, lent)
x13 = zeros(1, lent)
x14 = zeros(1, lent)
x15 = zeros(1, lent)
x16 = zeros(1, lent)
x17 = zeros(1, lent)
x18 = zeros(1, lent)
x19 = zeros(1, lent)
x20 = zeros(1, lent)
x21 = zeros(1, lent)
x22 = zeros(1, lent)
x23 = zeros(1, lent)
x24 = zeros(1, lent)
x25 = zeros(1, lent)
x26 = zeros(1, lent)
x27 = zeros(1, lent)
x28 = zeros(1, lent)
x29 = zeros(1, lent)
x30 = zeros(1, lent)
x31 = zeros(1, lent)
x32 = zeros(1, lent)
x33 = zeros(1, lent)
x34 = zeros(1, lent)
x35 = zeros(1, lent)
x36 = zeros(1, lent)
x37 = zeros(1, lent)
x38 = zeros(1, lent)
x39 = zeros(1, lent)
x40 = zeros(1, lent)
x41 = zeros(1, lent)
x42 = zeros(1, lent)
x43 = zeros(1, lent)
x44 = zeros(1, lent)
x45 = zeros(1, lent)
x46 = zeros(1, lent)
x47 = zeros(1, lent)
x48 = zeros(1, lent)
x49 = zeros(1, lent)
x50 = zeros(1, lent)
x51 = zeros(1, lent)

x1(1, :) = x(1)
x2(1, :) = x(2)
x3(1, :) = x(3)
x4(1, :) = x(4)
x5(1, :) = x(5)
x6(1, :) = x(6)
x7(1, :) = x(7)
x8(1, :) = x(8)
x9(1, :) = x(9)
x10(1, :) = x(10)
x11(1, :) = x(11)
x12(1, :) = x(12)
x13(1, :) = x(13)
x14(1, :) = x(14)
x15(1, :) = x(15)
x16(1, :) = x(16)
x17(1, :) = x(17)
x18(1, :) = x(18)
x19(1, :) = x(19)
x20(1, :) = x(20)
x21(1, :) = x(21)
x22(1, :) = x(22)
x23(1, :) = x(23)
x24(1, :) = x(24)
x25(1, :) = x(25)
x26(1, :) = x(26)
x27(1, :) = x(27)
x28(1, :) = x(28)
x29(1, :) = x(29)
x30(1, :) = x(30)
x31(1, :) = x(31)
x32(1, :) = x(32)
x33(1, :) = x(33)
x34(1, :) = x(34)
x35(1, :) = x(35)
x36(1, :) = x(36)
x37(1, :) = x(37)
x38(1, :) = x(38)
x39(1, :) = x(39)
x40(1, :) = x(40)
x41(1, :) = x(41)
x42(1, :) = x(42)
x43(1, :) = x(43)
x44(1, :) = x(44)
x45(1, :) = x(45)
x46(1, :) = x(46)
x47(1, :) = x(47)
x48(1, :) = x(48)
x49(1, :) = x(49)
x50(1, :) = x(50)
x51(1, :) = x(51)

//comet([x1', x2', x3', x4', x5', x6', x7', x8', x9', x10', x11', x12', x13', x14', x15', x16', x17', x18', x19', x20', x21', x22', x23', x24', x25', x26', x27', x28', x29', x30', x31', x32', x33', x34', x35', x36', x37', x38', x39', x40', x41', x42', x43', x44', x45', x46', x47', x48', x49', x50', x51'], [arr(:, 1), arr(:, 2), arr(:, 3), arr(:, 4), arr(:, 5), arr(:, 6), arr(:, 7), arr(:, 8), arr(:, 9), arr(:, 10), arr(:, 11), arr(:, 12), arr(:, 13), arr(:, 14), arr(:, 15), arr(:, 16), arr(:, 17), arr(:, 18), arr(:, 19), arr(:, 20), arr(:, 21), arr(:, 22), arr(:, 23), arr(:, 24), arr(:, 25), arr(:, 26), arr(:, 27), arr(:, 28), arr(:, 29), arr(:, 30), arr(:, 31), arr(:, 32), arr(:, 33), arr(:, 34), arr(:, 35), arr(:, 36), arr(:, 37), arr(:, 38), arr(:, 39), arr(:, 40), arr(:, 41), arr(:, 42), arr(:, 43), arr(:, 44), arr(:, 45), arr(:, 46), arr(:, 47), arr(:, 48), arr(:, 49), arr(:, 50), arr(:, 51)])

comet([x16', x17', x18', x19', x20', x21', x22', x23', x24', x25', x26', x27', x28', x29', x30', x31', x32', x33', x34', x35', x36'], [arr(:, 16), arr(:, 17), arr(:, 18), arr(:, 19), arr(:, 20), arr(:, 21), arr(:, 22), arr(:, 23), arr(:, 24), arr(:, 25), arr(:, 26), arr(:, 27), arr(:, 28), arr(:, 29), arr(:, 30), arr(:, 31), arr(:, 32), arr(:, 33), arr(:, 34), arr(:, 35), arr(:, 36)], "colors", ["black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black"])

// see the relevant paper for theoretical formulation.
