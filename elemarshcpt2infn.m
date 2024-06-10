function [V, infn] = elemarshcpt2infn(cptprofile,ptdata)
%% Converts EleMarsh CpT profile into V matrix and infusion profile
%
% [V, infn] = elemarshcpt2tci(cptprofile,ptdata)
% Syntax cptprofile is the time-CpT matrix (where time is in seconds and CpT is when CpT has been changing)
% ptdata can be either the ABW or a vector [age weight height sex]
%
% if the full patient data is used then the 6th column of V would be the
% Eleveld Ce

maxinfnrate = 12000;
tstep =  1; %default time step = 1 second

%% Extract pt data
if length(ptdata) == 1
    EleMarshABW = ptdata;
else
    age = ptdata(1);
    weight = ptdata(2);
    height = ptdata(3);
    sex = ptdata(4);
    bmi = weight/(height/100)^2;
    % EleMarsh ABW guess using PDreams linear regression algo (updated 29/12/23)
    if sex == 0     %female
        EleMarshABW = round(17.5+0.9912*weight-0.001305*weight^2+1.528e-6*weight^3+1.006e-4*height^2-3.690e-4*weight*bmi-0.2682*age+1.560e-3*age^2-0.003543*age*weight+2.322e-6*age*weight^2+0.001080*age*bmi-0.07786*bmi);
    else     %male
        EleMarshABW = round(16.07+0.9376*weight-0.001383*weight^2+1.684e-6*weight^3+1.292e-4*height^2-3.801e-4*weight*bmi-0.2617*age+1.614e-3*age^2-0.003841*age*weight+3.927e-6*age*weight^2+0.001340*age*bmi-0.09995*bmi);
    end
end
disp (['ABW = ' num2str(EleMarshABW) ' kg'])

k10m = 0.119;
k12m = 0.112;
k21m = 0.055;
k13m = 0.042;
k31m = 0.0033;
V1m = 0.228 * EleMarshABW;
V2m = 0.463 * EleMarshABW;
V3m = 2.893 * EleMarshABW;

%% Let's convert CpT profile into second by second
siz = size(cptprofile,1);
for i = 1:1:siz-1
    currenttimestep = cptprofile(i,1);
    nexttimestep = cptprofile(i+1,1);
    cptprofileformat(currenttimestep:1:nexttimestep-1,1) = (currenttimestep:1:nexttimestep-1)';
    cptprofileformat(currenttimestep:1:nexttimestep-1,2) = cptprofile(i,2);
end
cptprofileformat(cptprofile(siz,1):1:(cptprofile(siz,1)),1) = (cptprofile(siz,1):1:(cptprofile(siz,1)))';
cptprofileformat(cptprofile(siz,1):1:(cptprofile(siz,1)),2) = cptprofile(end,2);
cptprofile = cptprofileformat; %let's start using the re-formated CpT profile

%% Let's get into the calculations!
siz = size(cptprofile,1);
V = zeros(siz,6);
infn = zeros(siz,1);
V(:,1) = 1:1:siz;



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
infn = [(1:1:length(infn))' infn];

if length(ptdata) == 4
    % Do Eleveld Ce calculations
    V1 = 6.28*(weight/(weight + 33.6))/(0.675675675676);
    V2 = 25.5 * (weight/70)*exp(-0.0156*(age-35));
    V3 = 273*exp(-0.0138*age)*(sex*((0.88+(0.12)/(1+(age/13.4)^(-12.7)))*(9270*weight/(6680+216*weight/(height/100)^2)))+(1-sex)*((1.11+(-0.11)/(1+(age/7.1)^(-1.1)))*(9270*weight/(8780+244*weight/(height/100)^2))))/54.4752059601377;
    Cl1 = ((sex*1.79+(1-sex)*2.1)*((weight/70)^0.75)*((age*52.143+40)^9.06)/((age*52.143+40)^9.06+42.3^9.06))*exp(-0.00286*age);
    Cl2 = 1.75*(((25.5*(weight/70)*exp(-0.0156*(age-35)))/25.5)^0.75)*(1+1.3*(1-(age*52.143+40)/((age*52.143+40)+68.3)));
    Cl3 = 1.11*(((sex*((0.88+(0.12)/(1+(age/13.4)^(-12.7)))*(9270*weight/(6680+216*weight/(height/100)^2)))+(1-sex)*((1.11+(-0.11)/(1+(age/7.1)^(-1.1)))*(9270*weight/(8780+244*weight/(height/100)^2))))*exp(-0.0138*age)/54.4752059601377)^0.75)*((age*52.143+40)/((age*52.143+40)+68.3)/0.964695544);
    k10 = Cl1 / V1;
    k12 = Cl2/V1;
    k21 = Cl2/V2;
    k13 = Cl3/V1;
    k31 = Cl3/V3;
    ke0 = 0.146*(weight/70)^(-0.25);
    %store the V and k's into a vector
    VmatE = [V1 ; V2 ; V3];
    kmatE = [k10 ; k12 ; k21 ; k13 ; k31 ; ke0];
    Veleveld = pkmodel(infn(:,2), infn(:,1),VmatE, kmatE);
    V(:,6) = Veleveld(:,6);
end