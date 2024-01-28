function V = elemarshCe(infndata,ptdata)
% calculates the Eleveld Ce profile based on EleMarsh inputs
% CpT-time profile
% inputs:
% infndata is a 2col matrix where 1st column is time, 2nd column is the user CpT input
% ptdata is a vector [age weight height sex]

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
tmax = max(infndata(:,1));
tsteps = size(infndata,1);
Tmat = (1:1:tmax)';
V = zeros(tmax,6);
V(:,1) = Tmat(:,1);

%% Let's do calculation for first time step based on user Cp input

% N.B. if the first time entered is NOT 0 then assume pump was NOT running
% before this first time step i.e. leave everything equal to 0

if infndata(1,1) > 0
    infn = zeros(infndata(1,1),1);
else
    infn = [];
end

% Start calculation at the first time step
% extract user input for the first time step
timeinput = infndata(1,1);
cptinput = infndata(1,2);
durFirsttimestep = infndata(2,1) - infndata(1,1); %duration of first time step

% calculate the Marsh bolus using the empirical algo derived from Agila pumps
bolus = (0.2197*cptinput^2+15.693*cptinput)/70*weight;

%generate EleMarsh infusion regime
[Velemarsh, infncalc] = tci(cptinput,bolus,(timeinput:1:timeinput+durFirsttimestep)',VmatMarsh,kmatMarsh,maxinfnrate);

infncalc(1) = []; %get rid of the first zero (which is a "byproduct" of the tci function)

infn = [infn ; infncalc];

%% Loop over each user entry 
% Convert CpT change into EleMarsh infn regime, then convert infn regime
% into Eleveld Ce

Vstate = Velemarsh(end,:); %this is the final state of the EleMarsh model at the END of the first time step

for t = 2:1:length(infndata(:,1))-1
    timeinput = infndata(t,1);
    cptinput = infndata(t,2);
    durtimestep = infndata(t+1,1) - infndata(t,1); %duration of next time step
    Tmat = (timeinput:1:timeinput+durtimestep)';

    % generate EleMarsh infusion regime for the next time step
    [Velemarsh, infncalc] = titrate(cptinput,Vstate,Tmat,VmatMarsh,kmatMarsh,maxinfnrate,0);

    infncalc(1) = []; %get rid of the first zero of infn regime, which is a "byproduct" of the titrate fn
    infn = [infn ; infncalc];

    % store the final state of EleMarsh model
    Vstate = Velemarsh(end,:);

end

% plug EleMarsh infusion regime into Eleveld V matrix
V = pkmodel(infn,(1:1:infndata(end,1))',Vmat,kmat);

% get rid of the useless columns
V(:,2:1:5) = [];

% debug plot - comment out if required
plot(V(:,1), V(:,2), 'r-', infndata(:,1), infndata(:,2), 'bo')
