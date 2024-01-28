function [ce50new, bis] = agebis(age,ce)

cemat = 0:0.2:8;
basebis = 93;
bismat = zeros(1,length(cemat));
bismatmod = zeros(1,length(cemat));
bismatshift = zeros(1,length(cemat));

clf
%hold on

ceadjust = 1;
stimadjust = 0;

if nargin == 0
    clc
    disp ('Propofol Dreams v2.0')
    disp (' ')
    disp (' ')
    age = input ('Age (yr) ');
    ce = 4;
    ce50 = 3.08*exp(-0.00635*(age-35))*ceadjust;
    disp (['Ce50 = ' num2str(ce50)])
    disp (' ')
    disp (' ')
    bismod = input ('BIS value for adjustment ');
    cemod = input ('Ce value for adjustment ');
    %%calculate new PD curve by modifying Ce50
    if cemod > ce50
        ce50new = cemod/((basebis/bismod-1)^(1/1.47));
    else
        ce50new = cemod/((basebis/bismod-1)^(1/1.89));
    end

    %%calculate new PD curve by shifting
    if ce50new > ce50 %pt more tolerant - need to shift curve to the RIGHT
        if cemod > ce50
            ce50shift = ce50*((basebis/bismod-1)^(1/1.47));
        else
            ce50shift = ce50*((basebis/bismod-1)^(1/1.89));
        end
    else %pt more sensitive, shifting isn't appropriate
        ce50shift = 0;
    end


    for n = 1:1:length(cemat)
        if cemat(n) < stimadjust %stimadjust is an offset/shift factor, normally set it to 0
            bismat(n) = basebis;
        elseif cemat(n)-stimadjust > ce50
            bismat(n) = basebis*((ce50^1.47))/(ce50^1.47+(cemat(n)-stimadjust)^1.47);
        else
            bismat(n) = basebis*((ce50^1.89)/(ce50^1.89+(cemat(n)-stimadjust)^1.89));
        end
        if cemat(n) < stimadjust
            bismatmod(n) = basebis;
        elseif cemat(n)-stimadjust > ce50new
            bismatmod(n) = basebis*((ce50new^1.47))/(ce50new^1.47+(cemat(n)-stimadjust)^1.47);
        else
            bismatmod(n) = basebis*((ce50new^1.89)/(ce50new^1.89+(cemat(n)-stimadjust)^1.89));
        end
    end
    if ce50shift == 0
        plot(cemat,bismat, cemat, bismatmod, '--')
    else
        stimadjust = cemod - ce50shift;
        for n = 1:1:length(cemat)
            if cemat(n) < stimadjust %stimadjust is an offset/shift factor, normally set it to 0
                bismatshift(n) = basebis;
            elseif cemat(n)-stimadjust > ce50
                bismatshift(n) = basebis*((ce50^1.47))/(ce50^1.47+(cemat(n)-stimadjust)^1.47);
            else
                bismatshift(n) = basebis*((ce50^1.89)/(ce50^1.89+(cemat(n)-stimadjust)^1.89));
            end
        end
        plot(cemat,bismat, cemat, bismatmod, 'r--', cemat, bismatshift, 'k:')
    end
    xlabel('Effect site concentration')
    ylabel('Predicted BIS')
    bis = basebis*((ce50^1.47)/(ce50^1.47+ce^1.47));

else
    ce50 = 3.08*exp(-0.00635*(age-35));

    if ce >ce50
        bis = basebis*((ce50^1.47)/(ce50^1.47+ce^1.47));
    else
        bis = basebis*((ce50^1.89)/(ce50^1.89+ce^1.89));
    end
end




