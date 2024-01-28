function Velemarsh = elemarshCeLive(ptdata)
% calculates the Eleveld Ce profile based on EleMarsh inputs IN REAL TIME
% inputs:
% ptdata is a vector [age weight height sex]

% infndata is a 2col matrix where 1st column is time, 2nd column is the user CpT input
% in this function, we assume that the pt entered infndata IN REAL TIME

% this is the test dataset that pt entered
infndata = [1	5.1
505 5.6
1000 3];


% convert this into sec-by-sec dataset
cptprofile = zeros (infndata(end,1), 2);
cptprofile(:,1) = (1:1:size(cptprofile, 1))'; %create time steps in 1 sec intervals

for i = 1:1:length(infndata)-1
    numsteps = infndata(i+1,1) - infndata(i,1);
    cptprofile(infndata(i,1):1:infndata(i,1)+numsteps-1,2) = ones(numsteps,1)*infndata(i,2);
end
cptprofile(infndata(end,1), 2) = infndata(end,2);

maxinfnrate = 12000;

%% Extract pt data
age = ptdata(1);
weight = ptdata(2);
height = ptdata(3);
sex = ptdata(4);
bmi = weight/(height/100)^2;

%% Eleveld parameters
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
Vmat = [V1 ; V2 ; V3];
kmat = [k10 ; k12 ; k21 ; k13 ; k31 ; ke0];


%% EleMarsh constant parameters
k10m = 0.119;
k12m = 0.112;
k21m = 0.055;
k13m = 0.042;
k31m = 0.0033;

%% EleMarsh ABW guess using PDreams linear regression algo (updated 29/12/23)
if sex == 0     %female
    EleMarshABW = round(17.5+0.9912*weight-0.001305*weight^2+1.528e-6*weight^3+1.006e-4*height^2-3.690e-4*weight*bmi-0.2682*age+1.560e-3*age^2-0.003543*age*weight+2.322e-6*age*weight^2+0.001080*age*bmi-0.07786*bmi);
else     %male
    EleMarshABW = round(16.07+0.9376*weight-0.001383*weight^2+1.684e-6*weight^3+1.292e-4*height^2-3.801e-4*weight*bmi-0.2617*age+1.614e-3*age^2-0.003841*age*weight+3.927e-6*age*weight^2+0.001340*age*bmi-0.09995*bmi);
end

%check with user if they are happy with the PDreams App ABW (i.e. lrego)
chk = input(['Did you use ' num2str(EleMarshABW) 'kg as your EleMarsh adjusted body weight? (1 = Yes, 0 = No) ']);
if chk == 0
    EleMarshABW = input('Enter the EleMarsh adjusted body weight you used: ');
end
V1m = 0.228 * EleMarshABW;
V2m = 0.463 * EleMarshABW;
V3m = 2.893 * EleMarshABW;

kmatMarsh = [k10m ; k12m ; k21m ; k13m ; k31m ; 0];
VmatMarsh = [V1m ; V2m ; V3m];

%% Set up time matrix

tstep = 1; %1 second time steps
TT = 0; % the TT variable is the TIME KEEPING variable in this algo

%% Let's start looping "in real time" through each time step and do the following calculations steps
% (1) Read the user Cp input (from matrix cptprofile)
% (2) Use Cp into EleMarsh TCI model to calc infusion rate
% (3) Use infusion rate into Eleveld TCI model to calc instantaneous Ce
% (4) Report the Eleveld Ce

%% FIRST loop is for BEFORE user FIRST START the infusion, i.e. PUMP PAUSE

for i = 1:1:size(cptprofile,1)
    TT = TT + 1;
    disp (['Time step = ' num2str(TT)])
    Velemarsh(i, :) = zeros(1,5);
    Velemarsh(i, 1) = TT;
    Veleveld(i, :) = zeros(1,6);
    Veleveld(i, 1) = TT;
    infn(i, 1) = 0;
    disp (['CeT = ' num2str(Veleveld(TT,6))])
    disp (' ')
    if cptprofile(i,2) > 0 % detect a rise in CpT - i.e. user is starting induction
        break % break out of the pause loop and go into NEXT loop, which is for induction
    end
end

%% SECOND loop is when user starts the infusion

for t = TT+1:1:size(cptprofile,1)
    TT = TT + 1; % step in time
    Velemarsh(t,1) = TT;
    Veleveld(t,1) = TT;
    disp (['Time step = ' num2str(TT)])

    if TT == 2842
        disp(sum(infn)/3600)        
        Velemarsh(:,3) = Velemarsh(:,3) * 895 / (sum(infn)/3600);
        Velemarsh(:,4) = Velemarsh(:,3) * 895 / (sum(infn)/3600);
        infn = infn * 895 / (sum(infn)/3600);
        disp(sum(infn)/3600)
        pause
    end
    
    dV1 = (k21m*Velemarsh(t-1,3)+k31m*Velemarsh(t-1,4)-Velemarsh(t-1,2)*(k10m+k12m+k13m))*tstep/60; %delta V1 compartment from redistribution
    sugginfn = 3600 * (cptprofile(t,2)*V1m-dV1-Velemarsh(t-1,2)) / tstep;
    if sugginfn > maxinfnrate
        infn(t) = maxinfnrate;
    elseif sugginfn < 0
        infn(t) = 0;
    else
        infn(t) = sugginfn;
    end
    
    Velemarsh(t,2) = Velemarsh(t-1,2) + dV1 + infn(t)*tstep/3600; %iterative calc V1 drug
    Velemarsh(t,3) = Velemarsh(t-1,3) + (k12m*Velemarsh(t-1,2)-k21m*Velemarsh(t-1,3))*tstep/60; %iterative V2
    Velemarsh(t,4) = Velemarsh(t-1,4) + (k13m*Velemarsh(t-1,2)-k31m*Velemarsh(t-1,4))*tstep/60; %iterative V3
    Velemarsh(t,5) = Velemarsh(t,2)/V1m;
    % now calc eleveld tci for the Ce
    dV1 = (k21*Veleveld(t-1,3)+k31*Veleveld(t-1,4)-Veleveld(t-1,2)*(k10+k12+k13))*tstep/60; %delta V1 compartment from redistribution
    Veleveld(t,2) = Veleveld(t-1,2) + dV1 + infn(t)*tstep/3600; %iterative calc V1 drug
    Veleveld(t,3) = Veleveld(t-1,3) + (k12*Veleveld(t-1,2)-k21*Veleveld(t-1,3))*tstep/60; %iterative V2
    Veleveld(t,4) = Veleveld(t-1,4) + (k13*Veleveld(t-1,2)-k31*Veleveld(t-1,4))*tstep/60; %iterative V3
    Veleveld(t,5) = Veleveld(t,2)/V1;
    Veleveld(t,6) = Veleveld(t-1,6) + (Veleveld(t-1,5)-Veleveld(t-1,6))*ke0/60;
    disp (['CeT = ' num2str(Veleveld(TT,6))])
    disp (' ')
end

plot(Veleveld(:,1), Veleveld(:,6), 'r-', Velemarsh(:,1), Velemarsh(:,5), 'g.-', infndata(:,1), infndata(:,2), 'bo')