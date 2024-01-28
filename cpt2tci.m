function [V, infn] = cpt2tci(cptprofile,Vmat,kmat,maxinfnrate,Vbaseline)
%% Converts CpT profile into V matrix and infusion profile
V1m = Vmat(1);
V2m = Vmat(2);
V3m = Vmat(3);
k10m = kmat(1);
k12m = kmat(2);
k21m = kmat(3);
k13m = kmat(4);
k31m = kmat(5);
siz = size(cptprofile,1);

V = zeros(siz(1),5);
infn = zeros(siz(1),1);
V(:,1) = cptprofile(:,1);
tstep =  (cptprofile(2,1) - cptprofile(1,1));

if nargin == 4
    Vbaseline = V(1,:);
else
    V(1,:) = Vbaseline;
end


for t = 2:1:siz
    %start the infusion
    dV1 = (k21m*V(t-1,3)+k31m*V(t-1,4)-V(t-1,2)*(k10m+k12m+k13m))*tstep/60; %delta V1 compartment from redistribution
    sugginfn = 3600 * (cptprofile(t,2)*V1m-dV1-V(t-1,2)) / tstep;
    if sugginfn > maxinfnrate
        infn(t) = maxinfnrate;
    elseif sugginfn < 0
        infn(t) = 0;
    else
        infn(t) = sugginfn;
    end
    V(t,2) = V(t-1,2) + dV1 + infn(t)*tstep/3600; %iterative calc V1 drug
    V(t,3) = V(t-1,3) + (k12m*V(t-1,2)-k21m*V(t-1,3))*tstep/60; %iterative V2
    V(t,4) = V(t-1,4) + (k13m*V(t-1,2)-k31m*V(t-1,4))*tstep/60; %iterative V3
    V(t,5) = V(t,2)/V1m;
end