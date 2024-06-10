function Veleveld = elemarshinfn2plot (infn, ptdata)

%% Function that compares the decrement time for a given infusion regime
% infn is a 2 column matrix where:
% col 1 = time
% col 2 = infn rate
%
% if you want to use eBIS augmentation then enter peeg vector too

% set baseline parameters
stoptime = 0;

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

%% EleMarsh (ABW) parameters

% EleMarsh ABW guess using PDreams linear regression algo (updated 29/12/23)
if sex == 0     %female
    EleMarshABW = round(17.5+0.9912*weight-0.001305*weight^2+1.528e-6*weight^3+1.006e-4*height^2-3.690e-4*weight*bmi-0.2682*age+1.560e-3*age^2-0.003543*age*weight+2.322e-6*age*weight^2+0.001080*age*bmi-0.07786*bmi);
else     %male
    EleMarshABW = round(16.07+0.9376*weight-0.001383*weight^2+1.684e-6*weight^3+1.292e-4*height^2-3.801e-4*weight*bmi-0.2617*age+1.614e-3*age^2-0.003841*age*weight+3.927e-6*age*weight^2+0.001340*age*bmi-0.09995*bmi);
end
disp (EleMarshABW)

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

Veleveld = pkmodel(infn(:,2), infn(:,1),VmatE, kmatE);
Vmarsh = pkmodel(infn(:,2), infn(:,1),VmatM, kmatM);

plot(Veleveld(:,1)/60, Veleveld(:,6), 'r-', Vmarsh(:,1)/60, Vmarsh(:,5), 'b--')

xlabel ('Time (min)')
yyaxis left
ylabel ('Eleveld Ce (mcg/mL)', 'Color', 'r')
yyaxis right
ylabel ('EleMarsh CpT (mcg/mL)', 'Color', 'b')
ax = gca;
ax.YAxis(1).Color = 'r';
ax.YAxis(2).Color = 'b';
%xline(stoptime, '-', {'infusion stop'})
yyaxis left
ylabel ('Eleveld Ce (mcg/mL)', 'Color', 'r')
yyaxis right
ylabel ('EleMarsh CpT (mcg/mL)', 'Color', 'b')
legend ('Eleveld', 'EleMarsh', 'Location', 'north')
