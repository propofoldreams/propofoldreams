function elemarshcalc(cptprofile,ptdata,timebis)

% This is the user interface that unifies multiple elemarsh functions
% it basically serves as a daily driver
% Takes the inputs:
% (1) time-cpt profile
% (2) ptdata - [age weight height sex]
% (3) timebis - [time(seconds) bis]
% Analyses the above and produces the Eleveld Ce as well as augmented
% pharmacodynamic parameters (sigmoid Emax)

%% Extract pt data
age = ptdata(1);
weight = ptdata(2);
height = ptdata(3);
sex = ptdata(4);
bmi = weight/(height/100)^2;
lbm = (1.07*weight-148*(weight/height)^2)*(1-sex)+sex*(1.1*weight-128*(weight/height)^2); % James Equation for LBM

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
VmatE = [V1 ; V2 ; V3];
kmatE = [k10 ; k12 ; k21 ; k13 ; k31 ; ke0];


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
VmatM = [V1m ; V2m ; V3m];
kmatM = [k10m ; k12m ; k21m ; k13m ; k31m ; ke0m];

%% First convert the cptprofile into infusion regime and do Eleveld Ce calculations
[Velemarsh, infn] = elemarshcpt2infn(cptprofile,ptdata);
endtime = Velemarsh(end,1);

%% Let's do some Eleveld eBIS calculations based on pt
basebis = 93;
ce50 = 3.08*exp(-0.00635*(age-35));

if Velemarsh(end,6) > ce50
    endbis = basebis * ce50^1.47/(ce50^1.47 + Velemarsh(end,6)^1.47);
else
    endbis = basebis * ce50^1.89/(ce50^1.89 + Velemarsh(end,6)^1.89);
end

%% This is the PROPOFOL DREAMS dose-response augmentation algorithm
if nargin == 3 % only run this section if peeg vector is entered by user
    % take the timebis vector and augment the ce50 fitting; the ce50fit and
    % ce50fitshift values using the ce50calc function

    [ce50fit, ce50shiftfit] = ce50calc(Velemarsh(timebis(1),6), timebis(2));
else
    ce50fit = ce50;
    ce50shiftfit = 0;
end

%% WAKE UP TIME prediction algorithm
% pump stop time at bisread(k,1)
% we iterate until wake up Ce (i.e. eBIS = 60) is reached

% set up the baseline parameters
cebis30 = (93/30-1)^(1/1.89)*ce50fit+ce50shiftfit;
cebis40 = (93/40-1)^(1/1.89)*ce50fit+ce50shiftfit;
cebis50 = (93/50-1)^(1/1.47)*ce50fit+ce50shiftfit;
cebis60 = (93/60-1)^(1/1.47)*ce50fit+ce50shiftfit;
cebis70 = (93/70-1)^(1/1.47)*ce50fit+ce50shiftfit;
cebis80 = (93/80-1)^(1/1.47)*ce50fit+ce50shiftfit;
if Velemarsh(end,6) > ce50fit
    endbis = basebis * ce50fit^1.47/(ce50fit^1.47 + (Velemarsh(end,6)-ce50shiftfit)^1.47);
else
    endbis = basebis * ce50fit^1.89/(ce50fit^1.89 + (Velemarsh(end,6)-ce50shiftfit)^1.89);
end
Vwake = Velemarsh(end,:);

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

%% PLOT DATA
plot(Velemarsh(:,1)/60, Velemarsh(:,6), 'r-', Velemarsh(:,1)/60, Velemarsh(:,5), 'b--')
%yline(cebis70, 'r-.', {'eBIS70'})
xlabel ('Time (min)')
ylabel ('Concentration (mcg/mL)')
yregion(0,cebis70,FaceColor="r")
yregion(cebis70,cebis60,FaceColor="y")
yregion(cebis60,cebis40,FaceColor="g")
yregion(cebis40, cebis30, FaceColor="y")
legend ('Eleveld Ce', 'EleMarsh CpT', 'Location', 'northeast')


%% DISPLAY DATA
clc
disp ('-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=')
disp (' ')
disp ('Propofol Dreams EleMarsh Calculator v1.0')
disp (' ')
disp ('-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=')
disp (' ')
disp (['A' num2str(age) ' W' num2str(weight) ' H' num2str(height) ' Sex' num2str(sex)])
disp (['EleMarsh ABW' num2str(ABW)])
disp (' ')
disp (['Input end time: ' num2str(endtime) ' secs'])
disp (' ')
disp (['Baseline Eleveld Ce50: ' num2str(round(ce50*100)/100)])
if nargin == 3
    disp (['Augmented Eleveld Ce50: ' num2str(round(ce50fit*100)/100)])
    disp (['Augmented Eleveld Ce50 SHIFT: ' num2str(round(ce50shiftfit*100)/100)])
end
disp (' ')
disp (['Eleveld Ce = ' num2str(round(Velemarsh(end,6)*100)/100) ' [augBIS' num2str(round(endbis)) ']'])
disp (' ')
disp ('-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=')
disp (' ')
disp ('Decrement Times')
if nargin == 3
    disp (' ')
    disp (['Decrement time to eBIS60 = ' num2str(floor(decteleveldaug(1)/60)) ' min ' num2str(mod(decteleveldaug(1),60)) ' sec @ Ce = ' num2str(cebis60)])
    disp (['Decrement time to eBIS70 = ' num2str(floor(decteleveldaug(2)/60)) ' min ' num2str(mod(decteleveldaug(2),60)) ' sec @ Ce = ' num2str(cebis70)])
else
    disp (['Decrement time to eBIS60 = ' num2str(floor(decteleveldaug(1)/60)) ' min ' num2str(mod(decteleveldaug(1),60)) ' sec @ ' num2str(cebis60)])
    disp (['Decrement time to eBIS70 = ' num2str(floor(decteleveldaug(2)/60)) ' min ' num2str(mod(decteleveldaug(2),60)) ' sec @ ' num2str(cebis70)])
end
disp (' ')
disp ('-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=')
