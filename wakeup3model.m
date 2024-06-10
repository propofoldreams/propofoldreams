function wakeup3model (infn, ptdata, timebis)

%% Function that compares the decrement time for a given infusion regime
% infn is a 2 column matrix where:
% col 1 = time
% col 2 = infn rate
%
% timebis is a 2 row 2 column matrix where:
% col 1 = time | col 2 = SE or BIS
% row 1 = maintenance
% row 2 = when pEEG speed up

% set baseline parameters
stoptime = 0;

%% Extract pt data
age = ptdata(1);
weight = ptdata(2);
height = ptdata(3);
sex = ptdata(4);
bmi = weight/(height/100)^2;
lbm = (1.07*weight-148*(weight/height)^2)*(1-sex)+sex*(1.1*weight-128*(weight/height)^2); % James Equation for LBM

% Extract timebis matrix data
maintaintime = timebis(1,1);
maintainbis = timebis(1,2);
peegtime = timebis(2,1);
peegbis = timebis(2,2);

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
% Eleveld PD
basebis = 93;
ce50 = 3.08*exp(-0.00635*(age-35));
%store the V and k's into a vector
VmatE = [V1 ; V2 ; V3];
kmatE = [k10 ; k12 ; k21 ; k13 ; k31 ; ke0];

%% Schnider parameters
V1s = 4.27;
V2s = 18.9-0.391*(age-53);
V3s = 238;
k10s = (0.443+0.0107*(weight-77)-0.0159*(lbm-59)+0.0062*(height-177));
k12s = 0.0035;
k21s = (1.29-0.024*(age-53))/(18.9-0.391*(age-53));
k13s = 0.196;
k31s = 0.0035;
ke0s = 0.456;
%store the V and k's into a vector
VmatS = [V1s ; V2s ; V3s];
kmatS = [k10s ; k12s ; k21s ; k13s ; k31s ; ke0s];


%% EleMarsh (ABW) parameters
% EleMarsh ABW guess using PDreams linear regression algo (updated 29/12/23)
if sex == 0     %female
    EleMarshABW = round(17.5+0.9912*weight-0.001305*weight^2+1.528e-6*weight^3+1.006e-4*height^2-3.690e-4*weight*bmi-0.2682*age+1.560e-3*age^2-0.003543*age*weight+2.322e-6*age*weight^2+0.001080*age*bmi-0.07786*bmi);
else     %male
    EleMarshABW = round(16.07+0.9376*weight-0.001383*weight^2+1.684e-6*weight^3+1.292e-4*height^2-3.801e-4*weight*bmi-0.2617*age+1.614e-3*age^2-0.003841*age*weight+3.927e-6*age*weight^2+0.001340*age*bmi-0.09995*bmi);
end

V1m = 0.228 * EleMarshABW;
V2m = 0.463 * EleMarshABW;
V3m = 2.893 * EleMarshABW;
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

%% Infusion stop time

% Parse the infn regime to find when did pump stop

for i = 1:1:size(infn,1)
    if infn(end-i+1,2) == 0
        stoptime = stoptime + 1;
    else
        break
    end
end
stoptime = size(infn,1) - stoptime + 1;
wakeuptime = infn(end,1);
time2wakeup = wakeuptime - stoptime;


%% Let's calculate the relevant Ce/Cp at various time points for the 3 models
Veleveld = pkmodel(infn(:,2), infn(:,1),VmatE, kmatE);
Vschnider = pkmodel(infn(:,2), infn(:,1),VmatS, kmatS);
Vmarsh = pkmodel(infn(:,2), infn(:,1),VmatM, kmatM);

% Eleveld
eleveldmaintaince = Veleveld(maintaintime,6);
eleveldpeegce = Veleveld(peegtime,6);
eleveldwakece = Veleveld(end,6);

% EleMarsh
marshmaintaince = Vmarsh(maintaintime,5);
marshpeegce = Vmarsh(peegtime,5);
marshwakece = Vmarsh(end,5);


% Schnider
schnidermaintaince = Vschnider(maintaintime,6);
schniderpeegce = Vschnider(peegtime,6);
schniderwakece = Vschnider(end,6);


%% This is the PROPOFOL DREAMS dose-response augmentation algorithm
if nargin == 3 % only run this section if peeg vector is entered by user
    % take the timebis vector and augment the ce50 fitting; the ce50fit and
    % ce50fitshift values using the ce50calc function

    [ce50fit, ce50shiftfit] = ce50calc(eleveldmaintaince, maintainbis);

    % WAKE UP TIME prediction algorithm
    % pump stop time at bisread(k,1)
    % we iterate until wake up Ce (i.e. eBIS = 60) is reached

    % set up the baseline parameters
    % gamma = 1.47; % (this is from Eleveld population study)
    gamma = 3.93436825118642;
    basebis = 90;

    wakeupce60 = (basebis/60-1)^(1/gamma)*ce50fit+ce50shiftfit;
    wakeupce70 = (basebis/70-1)^(1/gamma)*ce50fit+ce50shiftfit;
    wakeupce80 = (basebis/80-1)^(1/gamma)*ce50fit+ce50shiftfit;
    eleveldbispeegaug = round(basebis*(ce50fit^gamma)/(ce50fit^gamma + (eleveldpeegce - ce50shiftfit)^gamma));
    eleveldbiswakeaug = round(basebis*(ce50fit^gamma)/(ce50fit^gamma + (eleveldwakece - ce50shiftfit)^gamma));
    Vwake = Veleveld(stoptime,:);

    t = 1;
    while Vwake(t,6) > wakeupce60
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

    while Vwake(t,6) > wakeupce70
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

    while Vwake(t,6) > wakeupce80
        t = t + 1;
        dV1 = (k21*Vwake(t-1,3)+k31*Vwake(t-1,4)-Vwake(t-1,2)*(k10+k12+k13))/60; %delta V1 compartment from redistribution
        Vwake(t,1) = Vwake(t-1,1) + 1;
        Vwake(t,2) = Vwake(t-1,2) + dV1; %iterative calc V1 drug
        Vwake(t,3) = Vwake(t-1,3) + (k12*Vwake(t-1,2)-k21*Vwake(t-1,3))/60; %iterative V2
        Vwake(t,4) = Vwake(t-1,4) + (k13*Vwake(t-1,2)-k31*Vwake(t-1,4))/60; %iterative V3
        Vwake(t,5) = Vwake(t,2)/V1;
        Vwake(t,6) = Vwake(t-1,6) + (Vwake(t-1,5)-Vwake(t-1,6))*ke0/60;
    end

    decteleveldaug(3) = t; %wake up time = how long after pump stopped before patient wakes up
else
    decteleveldaug = [0 0 0];
end


plot(Veleveld(:,1), Veleveld(:,6), 'r-', Vschnider(:,1), Vschnider(:,6), 'g-', Vmarsh(:,1), Vmarsh(:,5), 'b--')
xline(stoptime, '-', {'infusion stop'})
yline(wakeupce80, 'r-.', {'eBIS80'})
legend ('Eleveld', 'Schnider', 'EleMarsh', 'Location', 'north')

clc
disp ('-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=')
disp (' ')
disp (['A' num2str(age) ' W' num2str(weight) ' H' num2str(height) ' Sex' num2str(sex)])
disp (['EleMarsh ABW' num2str(EleMarshABW)])
disp (' ')
disp (['Infusion stop time: ' num2str(stoptime) ' secs'])
disp (['Patient wake up time: ' num2str(wakeuptime) ' secs'])
disp (['Wake up time: ' num2str(floor(time2wakeup/60)) ' min ' num2str(mod(time2wakeup,60)) ' sec'])
disp (' ')
disp (['Baseline Eleveld Ce50: ' num2str(round(ce50*100)/100)])
disp (['Augmented Eleveld Ce50: ' num2str(ce50fit)])
disp (' ')
disp ('-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=')
disp (' ')
disp ('Eleveld Model')
disp (['Maintenance Ce = ' num2str(round(eleveldmaintaince*1000)/1000) ' [BIS = ' num2str(maintainbis) ']'])
disp (['pEEG Speedup Ce = ' num2str(round(eleveldpeegce*1000)/1000) ' [eBIS = ' num2str(eleveldbispeegaug) ']'])
disp (['Wakeup Ce = ' num2str(round(eleveldwakece*1000)/1000) ' [eBIS = ' num2str(eleveldbiswakeaug) ']'])
disp (' ')
disp ('EleMARSH Model')
disp (['Maintenance Cp = ' num2str(round(marshmaintaince*1000)/1000)])
disp (['pEEG Speedup Cp = ' num2str(round(marshpeegce*1000)/1000)])
disp (['Wakeup Cp = ' num2str(round(marshwakece*1000)/1000)])
disp (' ')
disp ('Schnider Model')
disp (['Maintenance Ce = ' num2str(round(schnidermaintaince*1000)/1000)])
disp (['pEEG Speedup Ce = ' num2str(round(schniderpeegce*1000)/1000)])
disp (['Wakeup Ce = ' num2str(round(schniderwakece*1000)/1000)])
disp (' ')
disp ('-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=')
disp (' ')
disp ('Some wake up predictors')
disp (' ')
disp (['Actual time to wake = ' num2str(time2wakeup) ' sec'])
disp (' ')
disp (['Augmented Eleveld decrement time to eBIS60 = ' num2str(decteleveldaug(1)) ' @ ' num2str(wakeupce60)])
disp (['Augmented Eleveld decrement time to eBIS70 = ' num2str(decteleveldaug(2)) ' @ ' num2str(wakeupce70)])
disp (['Augmented Eleveld decrement time to eBIS80 = ' num2str(decteleveldaug(3)) ' @ ' num2str(wakeupce80)])
disp (' ')
disp ('-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=')
