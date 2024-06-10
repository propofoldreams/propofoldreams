function wakeupestimate (infnrate, cp, bis, ptdata)
%% An algorithm to estimate the Wake Up Ce and time from following inputs
% (1) infnrate = current infusion rate (in mL/hr)
% (2) cp = current EleMarsh Cp
% (3) bis = [Ce BIS]
% (4) ptdata = [age weight height sex]

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
%VmatE = [V1 ; V2 ; V3];
%kmatE = [k10 ; k12 ; k21 ; k13 ; k31 ; ke0];

%% EleMarsh parameters
% EleMarsh ABW guess using PDreams linear regression algo (updated 29/12/23)
if sex == 0     %female
    ABW = round(17.5+0.9912*weight-0.001305*weight^2+1.528e-6*weight^3+1.006e-4*height^2-3.690e-4*weight*bmi-0.2682*age+1.560e-3*age^2-0.003543*age*weight+2.322e-6*age*weight^2+0.001080*age*bmi-0.07786*bmi);
else     %male
    ABW = round(16.07+0.9376*weight-0.001383*weight^2+1.684e-6*weight^3+1.292e-4*height^2-3.801e-4*weight*bmi-0.2617*age+1.614e-3*age^2-0.003841*age*weight+3.927e-6*age*weight^2+0.001340*age*bmi-0.09995*bmi);
end
V1m = 0.228 * ABW;
V2m = 0.463 * ABW;
V3m = 2.893 * ABW;
k10m = 0.119;
k12m = 0.112;
k21m = 0.055;
k13m = 0.042;
k31m = 0.0033;
Cl1m = k10 * V1;
Cl2m = k21 * V2;
Cl3m = k31 * V3;
ke0m = 0;

%store the V and k's into a vector
%VmatM = [V1m ; V2m ; V3m];
%kmatM = [k10m ; k12m ; k21m ; k13m ; k31m ; ke0m];

% Convert infnrate into mg/second
infnrate_mlperhr = infnrate; %backup
infnrate = infnrate * 10 / 60;

%% PDreams WUT Estimation Algorithm

% PK COMPONENT
% the infn rate is what is needed to replenish the V1 compartment
% this is comprised of A1*(k10 + k12 + k13) - (A2*k21 + A3*k31)
% assumptions:
% (1) k31 is very small, let's assume A2 = A3 (generally A3 >> A2)
% (2) elemarsh has reached equilibrium so EleMarsh Cp = Eleveld Ce
% infnrate = Ce*V1*(k10 + k12 + k13) - (A2*k21 + A3*k31)
% thus: A2 = [Ce*V1* (k10 + k12 + k13) - infnrate]/(k21 + k31) [approx]
% plug A2 into 3 compartment model to work out decrement profile for A1

% PD COMPONENT
% as per the algorithm in ce50calc

% put above together, we get wake up times!
% back calculate to work out the EleMarsh Cp decrement threshold, et voila

%% Pharmacodynamics
gamma = 3.934368;
basebis = 93;
ce50 = 3.08*exp(-0.00635*(age-35));
[ce50fit, ce50shiftfit] = ce50calc(bis(1), bis(2));
cebis50 = (93/50-1)^(1/gamma)*ce50fit+ce50shiftfit;
cebis60 = (93/60-1)^(1/gamma)*ce50fit+ce50shiftfit;
cebis70 = (93/70-1)^(1/gamma)*ce50fit+ce50shiftfit;
cebis80 = (93/80-1)^(1/gamma)*ce50fit+ce50shiftfit;

%% Pharmacokinetics
A1 = cp*V1; % drug amount in V1 compartment
A2 = (A1*(k10+k12+k13) - infnrate)/(k21+k31); %drug amount in V2 compartment

Vwake = [1 A1 A2 A2 cp cp];

t = 1;
while Vwake(t,6) > cebis60
    t = t + 1;
    dV1 = (k21*Vwake(t-1,3)+k31*Vwake(t-1,4)-Vwake(t-1,2)*(k10+k12+k13))/60; %delta V1 compartment from redistribution
    Vwake(t,1) = Vwake(t-1,1) + 1;
    Vwake(t,2) = Vwake(t-1,2) + dV1; %iterative calc V1 drug
    Vwake(t,3) = Vwake(t-1,3) + (k12*Vwake(t-1,2)-k21*Vwake(t-1,3))/60; %iterative V2
    Vwake(t,4) = Vwake(t-1,4) + (k13*Vwake(t-1,2)-k31*Vwake(t-1,4))/60; %iterative V3
    Vwake(t,5) = Vwake(t,2)/V1;
    Vwake(t,6) = Vwake(t-1,6) + (Vwake(t-1,5)-Vwake(t-1,6))*ke0/60;
end

decteleveldaug(1) = t; %wake up time = how long after pump stopped before patient wakes up

while Vwake(t,6) > cebis70
    t = t + 1;
    dV1 = (k21*Vwake(t-1,3)+k31*Vwake(t-1,4)-Vwake(t-1,2)*(k10+k12+k13))/60; %delta V1 compartment from redistribution
    Vwake(t,1) = Vwake(t-1,1) + 1;
    Vwake(t,2) = Vwake(t-1,2) + dV1; %iterative calc V1 drug
    Vwake(t,3) = Vwake(t-1,3) + (k12*Vwake(t-1,2)-k21*Vwake(t-1,3))/60; %iterative V2
    Vwake(t,4) = Vwake(t-1,4) + (k13*Vwake(t-1,2)-k31*Vwake(t-1,4))/60; %iterative V3
    Vwake(t,5) = Vwake(t,2)/V1;
    Vwake(t,6) = Vwake(t-1,6) + (Vwake(t-1,5)-Vwake(t-1,6))*ke0/60;
end

decteleveldaug(2) = t; %wake up time = how long after pump stopped before patient wakes up

%% Convert Eleveld Ce to EleMarsh Cp

% Now that we worked out ETA BIS60 and ETA BIS70, we can work backwards and
% deduce the corresponding "wake up" EleMarsh Cp
% This is based on the assumption that at equilibrium:
% (1) EleMarsh infnrate = Eleveld infnrate
% (2) EleMarsh Cp = Eleveld Ce

A1m = cp*V1m; % drug amount in V1 compartment
A2m = (A1m*(k10m+k12m+k13m) - infnrate)/(k21m+k31m); %drug amount in V2 compartment
Vwakem = [1 A1m A2m A2m cp];

for t = 2:1:decteleveldaug(1)
    dV1m = (k21m*Vwakem(t-1,3)+k31m*Vwakem(t-1,4)-Vwakem(t-1,2)*(k10m+k12m+k13m))/60; %delta V1 compartment from redistribution
    Vwakem(t,1) = Vwakem(t-1,1) + 1;
    Vwakem(t,2) = Vwakem(t-1,2) + dV1m; %iterative calc V1 drug
    Vwakem(t,3) = Vwakem(t-1,3) + (k12m*Vwakem(t-1,2)-k21m*Vwakem(t-1,3))/60; %iterative V2
    Vwakem(t,4) = Vwakem(t-1,4) + (k13m*Vwakem(t-1,2)-k31m*Vwakem(t-1,4))/60; %iterative V3
    Vwakem(t,5) = Vwakem(t,2)/V1m;
end
cebis60m = Vwakem(end,5);

for t = decteleveldaug(1)+1:1:decteleveldaug(2)
    dV1m = (k21m*Vwakem(t-1,3)+k31m*Vwakem(t-1,4)-Vwakem(t-1,2)*(k10m+k12m+k13m))/60; %delta V1 compartment from redistribution
    Vwakem(t,1) = Vwakem(t-1,1) + 1;
    Vwakem(t,2) = Vwakem(t-1,2) + dV1m; %iterative calc V1 drug
    Vwakem(t,3) = Vwakem(t-1,3) + (k12m*Vwakem(t-1,2)-k21m*Vwakem(t-1,3))/60; %iterative V2
    Vwakem(t,4) = Vwakem(t-1,4) + (k13m*Vwakem(t-1,2)-k31m*Vwakem(t-1,4))/60; %iterative V3
    Vwakem(t,5) = Vwakem(t,2)/V1m;
end
cebis70m = Vwakem(end,5);

%% Let's do some duration calculations
% calculate ETA 20 mL and 50 mL

Vwakemdur = [1 A1m A2m A2m cp 0];

t = 1;
while Vwakemdur(end,6) < 20 % calculate ETA 20 mL
    t = t + 1;
    dV1m = (k21m*Vwakemdur(t-1,3)+k31m*Vwakemdur(t-1,4)-Vwakemdur(t-1,2)*(k10m+k12m+k13m))/60; %delta V1 compartment from redistribution
    sugginfn = 3600 * (cp*V1m-dV1m-Vwakemdur(t-1,2));
    Vwakemdur(t,1) = Vwakemdur(t-1,1) + 1;
    Vwakemdur(t,2) = Vwakemdur(t-1,2) + dV1m + sugginfn/3600; %iterative calc V1 drug
    Vwakemdur(t,3) = Vwakemdur(t-1,3) + (k12m*Vwakemdur(t-1,2)-k21m*Vwakemdur(t-1,3))/60; %iterative V2
    Vwakemdur(t,4) = Vwakemdur(t-1,4) + (k13m*Vwakemdur(t-1,2)-k31m*Vwakemdur(t-1,4))/60; %iterative V3
    Vwakemdur(t,5) = Vwakemdur(t,2)/V1m;
    Vwakemdur(t,6) = Vwakemdur(t-1,6) + sugginfn/3600/10;
end
eta20mL = t;

while Vwakemdur(end,6) < 50 % calculate ETA 50 mL
    t = t + 1;
    dV1m = (k21m*Vwakemdur(t-1,3)+k31m*Vwakemdur(t-1,4)-Vwakemdur(t-1,2)*(k10m+k12m+k13m))/60; %delta V1 compartment from redistribution
    sugginfn = 3600 * (cp*V1m-dV1m-Vwakemdur(t-1,2));
    Vwakemdur(t,1) = Vwakemdur(t-1,1) + 1;
    Vwakemdur(t,2) = Vwakemdur(t-1,2) + dV1m + sugginfn/3600; %iterative calc V1 drug
    Vwakemdur(t,3) = Vwakemdur(t-1,3) + (k12m*Vwakemdur(t-1,2)-k21m*Vwakemdur(t-1,3))/60; %iterative V2
    Vwakemdur(t,4) = Vwakemdur(t-1,4) + (k13m*Vwakemdur(t-1,2)-k31m*Vwakemdur(t-1,4))/60; %iterative V3
    Vwakemdur(t,5) = Vwakemdur(t,2)/V1m;
    Vwakemdur(t,6) = Vwakemdur(t-1,6) + sugginfn/3600/10;
end
eta50mL = t;

%% DISPLAY DATA
clc
disp ('-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=')
disp (' ')
disp ('Propofol Dreams EleMarsh Calculator v1.1')
disp (' ')
disp ('-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=')
disp (' ')
if sex == 1
    disp (['A' num2str(age) ' W' num2str(weight) ' H' num2str(height) ' Male'])
else
    disp (['A' num2str(age) ' W' num2str(weight) ' H' num2str(height) ' Female'])
end
disp (['EleMarsh ABW: ' num2str(ABW) ' kg'])
disp (['Current infusion rate: ' num2str(infnrate_mlperhr) ' mL/hr'])
disp (' ')
disp (['Baseline Eleveld Ce50: ' num2str(round(ce50*100)/100)])
disp (['Adjusted Eleveld Ce50: ' num2str(round(ce50fit*100)/100) ' | Shift: ' num2str(round(ce50shiftfit*100)/100)])
disp (' ')
disp ('-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=')
disp (' ')
disp ('Duration Mode')
disp (['ETA 20 mL = ' num2str(floor(eta20mL/60)) ' min ' num2str(mod(eta20mL,60)) ' sec'])
disp (['ETA 50 mL = ' num2str(floor(eta50mL/60)) ' min ' num2str(mod(eta50mL,60)) ' sec'])
disp (' ')
disp ('-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=')
disp (' ')
disp ('Decrement Times')
disp (['Decrement time to eBIS60 = ' num2str(floor(decteleveldaug(1)/60)) ' min ' num2str(mod(decteleveldaug(1),60)) ' sec'])
disp (['Eleveld Ce = ' num2str(round(cebis60,2)) '  |  EleMARSH Cp = ' num2str(round(cebis60m,2))])
disp (' ')
disp (['Decrement time to eBIS70 = ' num2str(floor(decteleveldaug(2)/60)) ' min ' num2str(mod(decteleveldaug(2),60)) ' sec'])
disp (['Eleveld Ce = ' num2str(round(cebis70,2)) '  |  EleMARSH Cp = ' num2str(round(cebis70m,2))])
disp (' ')
disp ('-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=')

%% Let's do some sexy plots!!

plot (Vwakem(:,1)/60, Vwakem(:,5), 'r-', Vwake(:,1)/60, Vwake(:,6), 'b--')
xlim([0 ceil(Vwake(end,1)/60)])
ylim ([floor(Vwakem(end,5)) ceil(Vwakem(1,5))])
ylabel ('Concentration (mcg/mL)')
xlabel ('Time (min)')
yregion(floor(min(Vwakem(:,5))), cebis70, FaceColor="red")
yregion(cebis70, cebis60, FaceColor="yellow")
yregion(cebis60, ceil(Vwakem(1,5)), FaceColor="green")
legend ('Cp', 'Ce', 'Location', 'northeast')