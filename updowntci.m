function [V, infn] = updowntci(bolus,Vmat,kmat,maxinfnrate,ttpe)
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


%crossovertime = 0;
pausetime = 0;
CpT = 4; %assume CpT = 4 to start

Tmat = (1:1:14400)';
siz = size(Tmat);
infn = zeros(siz(1),1);
infn(1) = maxinfnrate;
V = zeros(siz(1),6);
V(:,1) = Tmat(:,1);

if nargin == 6
    ttpe = 0;
end

%set up bolus
bolustime = bolus/maxinfnrate*3600;

% run this loop first time to optimse the induction regime to compensate
% for peaking function algorithm
if ttpe == 0 %CpT - skip this and just use the input bolus
else %optimise bolus
    for t = 2:1:3600
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
        bolus = bolus * CpT / V(end,6);
    end

    %reset everything
    V = zeros(siz(1),6);
    V(:,1) = Tmat(:,1);
    infn = zeros(siz(1),1);

    bolustime = bolus/maxinfnrate*3600;
end



%run this loop again for true induction
for t = 2:1:3600
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

%next hour droptarget CeT = 2
CpT = 2;

for t = 3601:1:7200
    tstep = 1; %1 sec time step
    %tstep = V(t,1)-V(t-1,1); %this is for the generic algo
    dV1 = (k21*V(t-1,3)+k31*V(t-1,4)-V(t-1,2)*(k10+k12+k13))*tstep/60; %delta V1 compartment from redistribution
    sugginfn = 3600 * (CpT*V1-dV1-V(t-1,2)) / tstep;
    if sugginfn < 0
        infn(t) = 0; %delta V1 compartment from infusion
    else
        infn(t) = sugginfn;
        if ttpe > 0 %i.e. we are in effect site targeting
            if V(t-1,6) > CpT %continue pausing infn until Ce reaches CeT
                infn(t) = 0;
            end
        end
    end
    V(t,2) = V(t-1,2) + dV1 + infn(t)*tstep/3600; %iterative calc V1 drug
    V(t,3) = V(t-1,3) + (k12*V(t-1,2)-k21*V(t-1,3))*tstep/60; %iterative V2
    V(t,4) = V(t-1,4) + (k13*V(t-1,2)-k31*V(t-1,4))*tstep/60; %iterative V3
    V(t,5) = V(t,2)/V1;
    V(t,6) = V(t-1,6) + (V(t-1,5)-V(t-1,6))*ke0/60;
end

%next hour, go back up to CeT 3
CpT = 3;

%figure out the bolus to get from target of 2 to 3
currinfn = infn(7200);
bolus = bolus/4 * (1 + currinfn/maxinfnrate);
bolustime = bolus/maxinfnrate*3600;

% if crossovertime == 0
%     %bolustime = round(bolus/4/maxinfnrate*3600) + 0.75*pausetime;
%     %bolus = bolus/4 + bolustime*currinfn/3600;
%     %complex plasma targeting w compensation algo
%     bolus = bolus/4; %simple bolus regime
%     bolustime = round(bolus/maxinfnrate*3600);
% else
%     bolustime = 0.5*round(bolus/4/maxinfnrate*3600) + 0.5*ttpe; %effect site targeting
%     bolus = bolus/4 + bolustime*currinfn/3600;
%     bolustime = round(bolus/maxinfnrate*3600);
% end


for t = 7201:1:10800
    tstep = 1; %1 sec time step
    %tstep = V(t,1)-V(t-1,1); %this is for the generic algo
    dV1 = (k21*V(t-1,3)+k31*V(t-1,4)-V(t-1,2)*(k10+k12+k13))*tstep/60; %delta V1 compartment from redistributio
    if V(t,1) <= 7200 + bolustime
        infn(t) = maxinfnrate;
    elseif V(t,1) <= 7201 + bolustime
        sugginfn = 3600 * (CpT*V1-dV1-V(t-1,2)) / tstep;
        if sugginfn < 0
            infn(t) = maxinfnrate * (7201 + bolustime - V(t,1));
        else
            infn(t) = maxinfnrate * (7201 + bolustime - V(t,1)) + sugginfn * (V(t,1) - bolustime - 7202);
        end
    else
        sugginfn = 3600 * (CpT*V1-dV1-V(t-1,2)) / tstep;
        if sugginfn < 0
            infn(t) = 0; %delta V1 compartment from infusion
        else
            infn(t) = sugginfn;
        end
    end
    V(t,2) = V(t-1,2) + dV1 + infn(t)*tstep/3600; %iterative calc V1 drug
    V(t,3) = V(t-1,3) + (k12*V(t-1,2)-k21*V(t-1,3))*tstep/60; %iterative V2
    V(t,4) = V(t-1,4) + (k13*V(t-1,2)-k31*V(t-1,4))*tstep/60; %iterative V3
    V(t,5) = V(t,2)/V1;
    V(t,6) = V(t-1,6) + (V(t-1,5)-V(t-1,6))*ke0/60;
end

%finally drop target CpT = 2
CpT = 2;

for t = 10801:1:14400
    tstep = 1; %1 sec time step
    %tstep = V(t,1)-V(t-1,1); %this is for the generic algo
    dV1 = (k21*V(t-1,3)+k31*V(t-1,4)-V(t-1,2)*(k10+k12+k13))*tstep/60; %delta V1 compartment from redistribution
    sugginfn = 3600 * (CpT*V1-dV1-V(t-1,2)) / tstep;
    if sugginfn < 0
        infn(t) = 0; %delta V1 compartment from infusion
    else
        infn(t) = sugginfn;
        if ttpe > 0 %i.e. we are in effect site targeting
            if V(t-1,6) > CpT %continue pausing infn until Ce reaches CeT
                infn(t) = 0;
            end
        end
    end
    V(t,2) = V(t-1,2) + dV1 + infn(t)*tstep/3600; %iterative calc V1 drug
    V(t,3) = V(t-1,3) + (k12*V(t-1,2)-k21*V(t-1,3))*tstep/60; %iterative V2
    V(t,4) = V(t-1,4) + (k13*V(t-1,2)-k31*V(t-1,4))*tstep/60; %iterative V3
    V(t,5) = V(t,2)/V1;
    V(t,6) = V(t-1,6) + (V(t-1,5)-V(t-1,6))*ke0/60;
end