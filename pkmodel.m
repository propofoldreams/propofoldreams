function V = pkmodel(infn,Tmat,Vmat,kmat,Vbaseline)
%this function takes an infusion regime and outputs the corresponding
%plasma and effect site concentrations for a given input model

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
V = zeros(siz(1),6);
V(:,1) = Tmat(:,1);

if nargin == 4
    Vbaseline = V(1,:);
else
    if length(Vbaseline) == 5
        V(1,:) = [Vbaseline 0];
    else
        V(1, :) = Vbaseline;
    end
end

for t = 2:1:siz(1)
    tstep = V(t,1)-V(t-1,1);
    dV1 = (k21*V(t-1,3)+k31*V(t-1,4)-V(t-1,2)*(k10+k12+k13))*tstep/60; %delta V1 compartment from redistribution    
    V(t,2) = V(t-1,2) + dV1 + infn(t)*tstep/3600; %iterative calc V1 drug
    V(t,3) = V(t-1,3) + (k12*V(t-1,2)-k21*V(t-1,3))*tstep/60; %iterative V2
    V(t,4) = V(t-1,4) + (k13*V(t-1,2)-k31*V(t-1,4))*tstep/60; %iterative V3
    V(t,5) = V(t,2)/V1;
    V(t,6) = V(t-1,6) + (V(t-1,5)-V(t-1,6))*ke0/60;
end