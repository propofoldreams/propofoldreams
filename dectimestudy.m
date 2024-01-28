function dectimestudy
%% CODE for the TCI decrement time study
% Basical we want to look at how various covariates affect the decrement
% time, including:
% Patient factors
% Infusion duration
% Ce50/Ce50shift

maxinfnrate = 12000;

%Young male BMI 23
ptdata = [50 120 173 1];

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


%% Set up time matrix
tmax = 6*3600+1800; %6.5 hr infusion
Tmat = (1:1:tmax)';
V = zeros(tmax,6);
V(:,1) = Tmat(:,1);

%% do CeBIS calculations
basebis = 93;
ce50 = 3.08*exp(-0.00635*(age-35));
cebis30 = (93/30-1)^(1/1.89)*ce50;
cebis40 = (93/40-1)^(1/1.89)*ce50;
cebis50 = (93/50-1)^(1/1.47)*ce50;
cebis60 = (93/60-1)^(1/1.47)*ce50;
cebis70 = (93/70-1)^(1/1.47)*ce50;
cebis80 = (93/80-1)^(1/1.47)*ce50;
disp (['BIS 30  = Ce ' num2str(round(cebis30*10)/10)])
disp (['BIS 40  = Ce ' num2str(round(cebis40*10)/10)])
disp (['BIS 50  = Ce ' num2str(round(cebis50*10)/10)])
disp (['BIS 60  = Ce ' num2str(round(cebis60*10)/10)])
disp (['BIS 70  = Ce ' num2str(round(cebis70*10)/10)])
disp (['BIS 80  = Ce ' num2str(round(cebis80*10)/10)])

%% Build the Eleveld model
maintcet = cebis30;
disp (' ')
disp (['Maintenance @ Ce = ' num2str(round(maintcet*10)/10)])
disp (' ')
[peakconc, ttpe] = peaking(Vmat,kmat);
bolus = maintcet/peakconc*10;

[Vgold, ~] = tci(maintcet,bolus,Tmat,Vmat,kmat,maxinfnrate,ttpe); %Eleveld TCI algo

disp ('Building Eleveld Model ... done!')

%% Calculate various decrements times

for i = 1:1:3
    if i == 1
        Vwake(1,:) = Vgold(1800, :); % Decrement after 30 minute infusion
    elseif i == 2
        Vwake(1,:) = Vgold(7200, :); % Decrement after 2 hr infusion
    else
        Vwake(1,:) = Vgold(21600, :); % Decrement after 6 hr infusion
    end

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

    dect60(i) = t; %this is decrement time to eBIS 60

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

    dect70(i) = t; %this is decrement time to eBIS 70

    while Vwake(t,6) > cebis80
        t = t + 1;
        dV1 = (k21*Vwake(t-1,3)+k31*Vwake(t-1,4)-Vwake(t-1,2)*(k10+k12+k13))/60; %delta V1 compartment from redistribution
        Vwake(t,1) = Vwake(t-1,1) + 1;
        Vwake(t,2) = Vwake(t-1,2) + dV1; %iterative calc V1 drug
        Vwake(t,3) = Vwake(t-1,3) + (k12*Vwake(t-1,2)-k21*Vwake(t-1,3))/60; %iterative V2
        Vwake(t,4) = Vwake(t-1,4) + (k13*Vwake(t-1,2)-k31*Vwake(t-1,4))/60; %iterative V3
        Vwake(t,5) = Vwake(t,2)/V1;
        Vwake(t,6) = Vwake(t-1,6) + (Vwake(t-1,5)-Vwake(t-1,6))*ke0/60;
    end

    dect80(i) = t; %this is decrement time to eBIS 80
end

disp ('Decrement time to 60')
disp (dect60/60)
disp (' ')
disp ('Decrement time to 70')
disp (dect70/60)
disp (' ')
disp ('Decrement time to 80')
disp (dect80/60)