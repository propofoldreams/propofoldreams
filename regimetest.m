%%this script tests performance using an up/down CeT regime
%Induce with CeT at 4 for 1 hr
% then reduce to 2 for 1 hr
% then increase to 3 for 1 hr
% then reduce to 2 for 1 hr
function [DevMatrix, Vgold, Vmarsh] = regimetest (ptdata,ABW,Abolus,maxinfnrate)

%baseline parameters, setting everything up
age = ptdata(1);
weight = ptdata(2);
height = ptdata(3);
sex = ptdata(4);

%sets up both Marsh and Eleveld parameters

%marsh constant parameters
k10m = 0.119;
k12m = 0.112;
k21m = 0.055;
k13m = 0.042;
k31m = 0.0033;
V1m = 0.228 * ABW;
V2m = 0.463 * ABW;
V3m = 2.893 * ABW;


%eleveld parameters
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

%store the V and k's into a vector for potential easy transfer if needed
Vmat = [V1 ; V2 ; V3];
kmat = [k10 ; k12 ; k21 ; k13 ; k31 ; ke0];
kmatMarsh = [k10m ; k12m ; k21m ; k13m ; k31m ; 0];
VmatMarsh = [V1m ; V2m ; V3m];

%set up initial conditions for everything else
Tmat = (1:1:14400)';

%% Let's build the gold standard Eleveld first

[peakconc, ttpe] = peaking(Vmat,kmat);
bolus = 4/peakconc*10; %bolus assume CeT = 4

[Vgold, ~] = updowntci(bolus,Vmat,kmat,maxinfnrate,ttpe); %Eleveld TCI algo

%% Next let's build the Marsh mimic

[~, infnMarsh] = updowntci(Abolus,VmatMarsh,kmatMarsh,maxinfnrate,0); %Marsh up down algo

Vmarsh = pkmodel(infnMarsh,Tmat,Vmat,kmat); %plug the marsh infn regime into Eleveld

%% Calc deviations

Cediff = Vmarsh(:,6) - Vgold(:,6);
PE = Cediff(5:end) ./ Vgold((5:1:end),6);
mdpe = round(1000*median(PE))/10; %MDPE calc
mdape = round(1000*median(abs(PE)))/10; %MDAPE calc
maxape = round(1000*max(abs(PE)))/10; %maxAPE calc
%record the result in the DevMatrix in order
DevMatrix(1) = mdpe;
DevMatrix(2) = mdape;
DevMatrix(3) = maxape;

%disp (DevMatrix) %suppress displays

%plot the best fit %suppress plots
%plot (Vgold(:,1), Vgold(:,6), 'r', Vmarsh(:,1), Vmarsh(:,6),'b--')
plot (Vmarsh(:,1)/60, Vmarsh(:,6),'color', [.7 .7 .7], 'linewidth', 0.2)
xlabel('Time (min)')
ylabel('Effect site concentration (mcg.mL^{-1})')
