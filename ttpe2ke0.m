function ke0 = ttpe2ke0(ttpe, Vmat, kmat)
% ttpe2ke0 algorithm
% this algorithm calculates the ke0 from the TTPE

V1 = Vmat(1);
V2 = Vmat(2);
V3 = Vmat(3);
k10 = kmat(1);
k12 = kmat(2);
k21 = kmat(3);
k13 = kmat(4);
k31 = kmat(5);
ke0 = kmat(6);

peakV = zeros(400,6);
peakV(:,1) = 1:1:400;
peakV(1,2) = 10; %put 10 units of drug into peakV at the start
for d = 2:1:400
    dpeakV1 = (k21*peakV(d-1,3)+k31*peakV(d-1,4)-peakV(d-1,2)*(k10+k12+k13))/60; %delta peakV1 compartment from redistribution
    peakV(d,2) = peakV(d-1,2) + dpeakV1;
    peakV(d,3) = peakV(d-1,3) + (k12*peakV(d-1,2)-k21*peakV(d-1,3))/60; %iterative V2
    peakV(d,4) = peakV(d-1,4) + (k13*peakV(d-1,2)-k31*peakV(d-1,4))/60; %iterative V3
    peakV(d,5) = peakV(d,2)/V1; %plasma conc
    peakV(d,6) = peakV(d-1,6) + (peakV(d-1,5)-peakV(d-1,6))*ke0/60; %equilibration into effect site
end
[peakconc, ttpeloc] = max(peakV(:,6));
ttpe = peakV(ttpeloc, 1);

%% work in progress