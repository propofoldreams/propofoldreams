function [V, infn] = tci(CpT,bolus,Tmat,Vmat,kmat,maxinfnrate,ttpe)
%The function **tci** generates takes
% input = (CpT,bolus,Tmat,Vmat,kmat,maxinfnrate,+/-ttpe)
% output = [drug in each compartment ; infusion rate] at each given time step as
% per the time matrix, Tmat
%
%baseline parameters, setting everything up
V1 = Vmat(1);
V2 = Vmat(2);
V3 = Vmat(3);
k10 = kmat(1);
k12 = kmat(2);
k21 = kmat(3);
k13 = kmat(4);
k31 = kmat(5);
ke0 = kmat(6);
siz = size(Tmat);
infn = zeros(siz(1),1);

V = zeros(siz(1),6);
V(:,1) = Tmat(:,1);

if nargin == 6
    ttpe = 0; %we are in plasma target mode!
end

%set up bolus
bolustime = bolus/maxinfnrate*3600;

% run this loop first time to optimse the induction regime to compensate
% for peaking function algorithm
if ttpe == 0 %CpT - skip this and just use the input bolus
else %optimise bolus
    for t = 2:1:siz(1)
        tstep = V(t,1)-V(t-1,1);
        dV1 = (k21*V(t-1,3)+k31*V(t-1,4)-V(t-1,2)*(k10+k12+k13))*tstep/60; %delta V1 compartment from redistribution
        if V(t,1) <= bolustime + 1
            infn(t) = maxinfnrate;
        elseif V(t,1) <= bolustime + 2
            sugginfn = 3600 * (CpT*V1-dV1-V(t-1,2)) / tstep;
            if sugginfn < 0
                infn(t) = maxinfnrate * (bolustime + 2 - V(t,1));
            else
                infn(t) = maxinfnrate * (bolustime + 2 - V(t,1)) + sugginfn * (V(t,1) - bolustime - 1);
            end
        else
            sugginfn = 3600 * (CpT*V1-dV1-V(t-1,2)) / tstep;
            if sugginfn < 0
                infn(t) = 0; %delta V1 compartment from infusion
            else
                V = V(1:1:t-1,:);
                break;
            end
        end
        V(t,2) = V(t-1,2) + dV1 + infn(t)*tstep/3600; %iterative calc V1 drug
        V(t,3) = V(t-1,3) + (k12*V(t-1,2)-k21*V(t-1,3))*tstep/60; %iterative V2
        V(t,4) = V(t-1,4) + (k13*V(t-1,2)-k31*V(t-1,4))*tstep/60; %iterative V3
        V(t,5) = V(t,2)/V1;
        V(t,6) = V(t-1,6) + (V(t-1,5)-V(t-1,6))*ke0/60;
    end

    if V(end,6) ~= 0 %prevents error with super short bolustimes (e.g. w remi)
        %disp(V(end,6)) %debug lines
        bolus = bolus * CpT / V(end,6);
    end

    %reset everything
    %plot(V(:,1), V(:,6)) %debug lines
    V = zeros(siz(1),6);
    V(:,1) = Tmat(:,1);
    infn = zeros(siz(1),1);

    bolustime = bolus/maxinfnrate*3600;
end



%run this loop again for true induction
for t = 2:1:siz(1) %theoretically in the future, I can code this into accepting non-zero starting V states
    tstep = V(t,1)-V(t-1,1);
    dV1 = (k21*V(t-1,3)+k31*V(t-1,4)-V(t-1,2)*(k10+k12+k13))*tstep/60; %delta V1 compartment from redistribution
    if V(t,1) <= bolustime + 1
        infn(t) = maxinfnrate;
    elseif V(t,1) <= bolustime + 2
        sugginfn = 3600 * (CpT*V1-dV1-V(t-1,2)) / tstep;
        if sugginfn < 0
            infn(t) = maxinfnrate * (bolustime + 2 - V(t,1));
        else
            infn(t) = maxinfnrate * (bolustime + 2 - V(t,1)) + sugginfn * (V(t,1) - bolustime - 1);
        end
    else
        sugginfn = 3600 * (CpT*V1-dV1-V(t-1,2)) / tstep;
        if sugginfn < 0
            infn(t) = 0; %delta V1 compartment from infusion
        else
            infn(t) = sugginfn;
            if ttpe > 0 %i.e. we are in effect site targeting
                if V(t,6) > CpT %continue pausing infn until Ce reaches CeT
                    infn(t) = 0;
                end
            end
        end
    end
    V(t,2) = V(t-1,2) + dV1 + infn(t)*tstep/3600; %iterative calc V1 drug
    V(t,3) = V(t-1,3) + (k12*V(t-1,2)-k21*V(t-1,3))*tstep/60; %iterative V2
    V(t,4) = V(t-1,4) + (k13*V(t-1,2)-k31*V(t-1,4))*tstep/60; %iterative V3
    V(t,5) = V(t,2)/V1;
    V(t,6) = V(t-1,6) + (V(t-1,5)-V(t-1,6))*ke0/60;
end