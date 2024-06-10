function [ce50, ce50shift] = ce50calc(ce, bis, out)
% this function takes the user input Ce and corresponding BIS (single point!) and
% uses the analytical Sigmoid Emax equation to calculate Ce50 and Ce50shift
% if out == 1, we also plot the Sig Emax curve

% This is the old algorithm for Ce50 calc
% basebis = 93;
% shiftratio = 0.3; % this is the ratio of ce50shift/ce50 - a novel Propofol Dream concept :)
% gamma = 1.47;

% Updated Ce50calc algorithm 28/5/2024
basebis = 90;
shiftratio = 0;
gamma = 3.934368;


if nargin == 2
    out = 0;
end

ce50 = ce / ((basebis/bis-1)^(1/gamma)+shiftratio);
ce50shift = ce50 * shiftratio;

if out > 0
    cemat = [0:0.01:6];

    for i = 1:1:length(cemat)
        if cemat(i) - ce50shift <= 0
            bisplot(i) = basebis;
        elseif cemat(i) - ce50shift > ce50
            bisplot(i) = basebis*((ce50^1.47))/(ce50^1.47+(cemat(i)-ce50shift)^1.47);
        else
            bisplot(i) = basebis*((ce50^1.89))/(ce50^1.89+(cemat(i)-ce50shift)^1.89);
        end
    end
    plot(ce, bis, 'ro', cemat, bisplot, 'g-')
end