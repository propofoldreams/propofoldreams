%this function gives the infusion regime for a given patient and model
confirm = 0;
maxpumprate = 1200; %mL/hr
tmax = 3600*4; %how long simulation runs in seconds

%set up variables
conc = zeros (tmax,1);
Tmat = (1:1:tmax)';
infn = zeros(tmax,1);
V = zeros(tmax,6);
V(:,1) = Tmat(:,1);

%% GETTING PATIENT DEETS
while confirm == 0
    clc
    disp ('Propofol Dreams v2.0')
    disp (' ')
    disp (' ')
    age = input ('Age (yr) ');
    weight = input ('Weight (kg) ');
    height = input ('Height (cm) ');
    sex = input ('Sex (0 = XX, 1 = XY) ');
    disp (' ')
    disp ('Model selection')
    disp ('0 = PROPOFOL Marsh')
    disp ('1 = PROPOFOL Eleveld')
    disp ('2 = REMIFENTANIL Eleveld')
    disp ('3 = DEXMED Hannivoort')    
    model = input ('Please select a model: ');
    if model == 0
        target = input ('Plasma target (mcg/mL) ');
    elseif model == 1
        target = input ('Effect site target (mcg/mL) ');
    elseif model == 2
        target = input ('Effect site target (ng/mL) ');
    elseif model == 3
        disp (' ')
        disp ('Light sedation = 0.2 - 0.3 ng/mL')
        disp ('Deep sedation > 1.9 ng/mL')
        disp ('Maximum 3 ng/mL')
        disp (' ')
        target = input ('Plasma target (ng/mL) ');
    end
    %confirm pt inputs
    clc
    disp (' ')
    disp ('Are the following inputs correct?')
    disp (' ')
    disp (['Age: ', num2str(age), ' yr'])
    disp (['Weight: ', num2str(weight), ' kg'])
    disp (['Height: ', num2str(height), ' cm'])
    if sex == 0
        disp ('Sex: Female')
    else
        disp ('Sex: Male')
    end
    if model == 0
        disp ('Model: Marsh')
        disp (['CpT: ' num2str(target) ' mcg/mL'])
    elseif model == 1
        disp ('Model: PROPOFOL Eleveld')
        disp (['CeT: ' num2str(target) ' mcg/mL'])
    elseif model == 2
        disp ('Model: REMIFENTANIL Eleveld')
        disp ('Concentration = 50 mcg/mL')
        disp (['CeT: ' num2str(target) ' ng/mL'])
    elseif model == 3
        disp ('Model: DEXMED Hannivoort')
        disp ('Concentration = 4 mcg/mL')
        disp (['CpT: ' num2str(target) ' ng/mL'])
    end
    disp (' ')
    confirm = input ('1 = Yes, 0 = No ');
    clc
end

%% SETTING UP MODEL PARAMETERS
%disp ('Model Parameters')
%disp (' ')

if model == 0
    drugconc = 10; %mg/mL
    maxhandrate = 2400*drugconc; %mg/hr (assuming double of pump)
    maxinfnrate = maxpumprate * drugconc; %mg/hr
    infn(1) = maxhandrate;

    V1 = 0.228 * weight;
    V2 = 0.463 * weight;
    V3 = 2.893 * weight;
    k10 = 0.119;
    k12 = 0.112;
    k21 = 0.055;
    k13 = 0.042;
    k31 = 0.0033;
    Cl1 = k10 * V1;
    Cl2 = k21 * V2;
    Cl3 = k31 * V3;
    ke0 = 0;
    %marsh parameters
elseif model == 1
    drugconc = 10; %mg/mL
    maxhandrate = 2400*drugconc; %mg/hr (assuming double of pump)
    maxinfnrate = maxpumprate * drugconc; %mg/hr
    infn(1) = maxhandrate;

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
elseif model == 2 %Eleveld Remifentanil
    drugconc = 50; %mcg/mL
    maxhandrate = 2400*drugconc; %mcg/hr (assuming double of pump)
    maxinfnrate = maxpumprate*drugconc; %mcg/hr
    infn(1) = maxhandrate;

    %reference individual = 35 yo, male, 70 kg, 170 cm
    ffm = sex*(0.88+(1-0.88)/(1+(age/13.4)^(-12.7)))*(9270*weight/(6680+216*(weight/(height/100)^2))) + (1-sex)*(1.11+(1-1.11)/(1+(age/7.1)^(-1.1)))*(9270*weight/(8780+244*(weight/(height/100)^2)));
    siz = ffm / 54.4752;
    ksex = sex + (1-sex)*(1+(0.470*age^6/(age^6 + 12^6))*(1-age^6/(age^6 + 45^6)));

    V1 = 5.81 * siz * exp(-0.00554*(age-35));
    V2 = 8.82 * siz * exp(-0.00327*(age-35))*ksex;
    V3 = 5.03 * siz * exp(-0.0315*(age-35))*exp(-0.0360*(weight-70));
    Cl1 = 2.58 * siz^0.75 * (weight^2/(weight^2+2.88^2)/(70^2/(70^2+2.88^2))) * ksex * exp(-0.00327*(age-35));
    Cl2 = 1.72 * (V2 / 8.82)^0.75 * exp(-0.00554*(age-35)) * ksex;
    Cl3 = 0.124 * (V3 / 5.03)^0.75 * exp(-0.00554*(age-35));

    k10 = Cl1 / V1;
    k12 = Cl2/V1;
    k21 = Cl2/V2;
    k13 = Cl3/V1;
    k31 = Cl3/V3;
    ke0 = 1.09*exp(-0.0289*(age-35));
elseif model == 3 %Hannivoort DEXMED
    drugconc = 4; %mcg/mL
    maxinfnrate = min(maxpumprate*drugconc, 6*weight*3600); %either maxpump rate but limit to 6 mcg/kg/hr
    maxhandrate = maxinfnrate; %same as maxinfnrate

    infn(1) = maxhandrate;

    V1 = 1.78 * weight/70;
    V2 = 30.3 * weight/70;
    V3 = 52.0 * weight/70;
    Cl1 = 0.686 * (weight/70)^0.75;
    Cl2 = 2.98 * (V2/30.3)^0.75;
    Cl3 = 0.602 * (V3/52.0)^0.75;

    k10 = Cl1/V1;
    k12 = Cl2/V1;
    k21 = Cl2/V2;
    k13 = Cl3/V1;
    k31 = Cl3/V3;
    ke0 = 0;
end
CpT = target;

%store the V and k's into a vector
Vmat = [V1 ; V2 ; V3];
kmat = [k10 ; k12 ; k21 ; k13 ; k31 ; ke0];

if model == 0
    disp ('Performing PROPOFOL Marsh CpT calculations')
    disp (' ')
    bolus = CpT*V1;
    ttpe = 0;
    [Vgold, Infngold] = tci(CpT,bolus,(1:1:900)',Vmat,kmat,maxinfnrate); %build the gold standard
elseif model == 1
    disp ('Performing PROPOFOL Eleveld CeT calculations')
    disp (' ')
    [peakconc, ttpe] = peaking(Vmat,kmat);
    bolus = CpT/peakconc*10;
    [Vgold, Infngold] = tci(CpT,bolus,(1:1:900)',Vmat,kmat,maxinfnrate,ttpe); %build the gold standard
elseif model == 2
    disp ('Performing REMIFENTANIL Eleveld CeT calculations')
    disp (' ')
    [peakconc, ttpe] = peaking(Vmat,kmat);
    bolus = CpT/peakconc*10;
    [Vgold, Infngold] = tci(CpT,bolus,(1:1:900)',Vmat,kmat,maxinfnrate,ttpe); %build the gold standard
elseif model == 3
    disp ('Performing DEXMED Hannivoort CpT calculations')
    disp (' ')
    bolus = CpT*V1;
    ttpe = 0;
    [Vgold, Infngold] = tci(CpT,bolus,(1:1:900)',Vmat,kmat,maxinfnrate); %build the gold standard
end

pumpbolustime = bolus/maxinfnrate*3600; %induction via pump bolusing

%% OPTIMISATION ALGO for First 15 minutes!
% let's sort out what happens in the first 15 min
%basically we need to work out a bolus + infusion regime that eliminates
%having to pause the pump during Ce targeting (i.e. makes it easier for users)

%to achieve this, let's test various bolus and infusion strategies by
%BRUTE FORCE!
%of course, we can also use a regression guess strategy to save computation time

VGuess = zeros(900,6);
VGuess(:,1) = (1:1:900)';

maxEbolus = ceil(bolus/drugconc)*drugconc; %nearest 10 mg/mcg bolus
minEbolus = round(bolus*0.85/drugconc)*drugconc;

%get rid of pauses
Infngold(Infngold == 0) = [];
Infngold(Infngold==maxinfnrate) = [];
maxInfn = ceil(sum(Infngold)/length(Infngold)/drugconc)*drugconc;
minInfn = round(Infngold(end)/drugconc)*drugconc;

%calc size of guess matrix
gsize = ((maxEbolus - minEbolus)/drugconc+1)*((maxInfn - minInfn)/drugconc+1);
GuessMatrix = zeros(gsize,3);
guessindex = 0;

%start guess algorithm
for BolusGuess = minEbolus:drugconc:maxEbolus
    bolustime = BolusGuess/maxhandrate*3600; %induction via hand bolusing
    for InfnGuess = minInfn:drugconc:maxInfn
        guessindex = guessindex + 1;
        % TCI-ing!
        for t = 2:1:900
            dV1 = (k21*V(t-1,3)+k31*V(t-1,4)-V(t-1,2)*(k10+k12+k13))/60; %delta V1 compartment from redistribution
            if VGuess(t,1) <= bolustime + 1
                infn(t) = maxhandrate;
            elseif VGuess(t,1) <= bolustime + 2
                infn(t) = maxhandrate * (bolustime + 2 - V(t,1)) + InfnGuess * (V(t,1) - bolustime - 1);
            else
                infn(t) = InfnGuess;
            end
            VGuess(t,2) = VGuess(t-1,2) + dV1 + infn(t)/3600; %iterative calc V1 drug
            VGuess(t,3) = VGuess(t-1,3) + (k12*VGuess(t-1,2)-k21*VGuess(t-1,3))/60; %iterative V2
            VGuess(t,4) = VGuess(t-1,4) + (k13*VGuess(t-1,2)-k31*VGuess(t-1,4))/60; %iterative V3
            VGuess(t,5) = VGuess(t,2)/V1;
            VGuess(t,6) = VGuess(t-1,6) + (VGuess(t-1,5)-VGuess(t-1,6))*ke0/60;

            %diff
            Cediff = VGuess(:,6) - Vgold(:,6);
            PE = Cediff(5:end) ./ Vgold((5:1:end),6);
            sumsqr = sqrt(sum(Cediff.^2)/length(Cediff)); %sumsqr calc
        end

        %record the result in the GuessMatrix in order
        GuessMatrix(guessindex,1) = BolusGuess;
        GuessMatrix(guessindex,2) = InfnGuess;
        GuessMatrix(guessindex,3) = sumsqr;
    end
end

%disp(GuessMatrix) %show the brute force results - for debug purposes only
%pause
minCe = min(GuessMatrix(:,3));
minloc = find (GuessMatrix(:,3)==minCe);

bolus = GuessMatrix (minloc,1);

if size(bolus,1) > 1
    bolus = bolus(end); %if there's more than 1 result, use the largest bolus
    minloc = minloc(end);
end

% infn algorithm for first "900" seconds (technically 899)
bolustime = bolus/maxhandrate*3600;
infn(1:1:floor(bolustime)) = maxhandrate*ones(floor(bolustime),1);
infn(ceil(bolustime)) = maxhandrate * (bolustime-floor(bolustime)) + GuessMatrix(minloc,2)*(1-(bolustime-floor(bolustime)));
infn(ceil(bolustime)+1:1:900) = GuessMatrix(minloc,2)*ones(900-ceil(bolustime),1);


%%
%pausetime = 0;
for t = 2:1:900
    dV1 = (k21*V(t-1,3)+k31*V(t-1,4)-V(t-1,2)*(k10+k12+k13))/60; %delta V1 compartment from redistribution
    V(t,2) = V(t-1,2) + dV1 + infn(t)/3600; %iterative calc V1 drug
    V(t,3) = V(t-1,3) + (k12*V(t-1,2)-k21*V(t-1,3))/60; %iterative V2
    V(t,4) = V(t-1,4) + (k13*V(t-1,2)-k31*V(t-1,4))/60; %iterative V3
    V(t,5) = V(t,2)/V1;
    V(t,6) = V(t-1,6) + (V(t-1,5)-V(t-1,6))*ke0/60;
end

%   if ttpe > 0 %effect site target for PROPOFOL
%   bolus = bolus * (-0.0008*age+0.8838); % let's scale the bolus with an age factor
%   keep this commented, use this linear approx for approx functions later

%% Averaging algo
% Keep cycling the process of
% (1) find the average infusion rate in 15 min chunks then
% (2) apply the avg infusion rate to the 15 min chunk

starttime = 900;
while starttime+899 < tmax
    for t = starttime:1:starttime+900
        dV1 = (k21*V(t-1,3)+k31*V(t-1,4)-V(t-1,2)*(k10+k12+k13))/60;
        infn(t) = 3600 * (CpT*V1-dV1-V(t-1,2)); %this is the suggested infusion for next time step

        V(t,2) = V(t-1,2) + dV1 + infn(t)/3600; %iterative calc V1 drug
        V(t,3) = V(t-1,3) + (k12*V(t-1,2)-k21*V(t-1,3))/60; %iterative V2
        V(t,4) = V(t-1,4) + (k13*V(t-1,2)-k31*V(t-1,4))/60; %iterative V3
        V(t,5) = V(t,2)/V1; %calc Cp
        V(t,6) = V(t-1,6) + (V(t-1,5)-V(t-1,6))*ke0/60; %calc Ce
    end

    avginfn = drugconc*round(sum(infn(starttime:1:starttime+899))/900/drugconc);
    if avginfn < 5 %if less than 5, then we want accuracy go up to 1 dp
        avginfn = drugconc/10*round(sum(infn(starttime:1:starttime+899))/900/drugconc*10);
    else
        avginfn = drugconc*round(sum(infn(starttime:1:starttime+899))/900/drugconc);
    end

    infn(starttime:1:starttime+900) = avginfn * ones(901,1);
    starttime = starttime + 900;
end

%% Let's start displaying the data
clc
disp ('Propofol Dreams v2.0')
disp (' ')
disp (['Age: ', num2str(age), ' yr'])
disp (['Weight: ', num2str(weight), ' kg'])
disp (['Height: ', num2str(height), ' cm'])
if sex == 0
    disp ('Sex: Female')
else
    disp ('Sex: Male')
end
disp (' ')
if model == 0
    disp ('Model: PROPOFOL Marsh')
    disp (['CpT: ' num2str(target) ' mcg/mL'])
elseif model == 1
    disp ('Model: PROPOFOL Eleveld')
    disp (['CeT: ' num2str(target) ' mcg/mL'])
elseif model == 2
    disp ('Model: REMIFENTANIL Eleveld')
    disp ('Concentration = 50 mcg/mL')
    disp (['CeT: ' num2str(target) ' ng/mL'])
elseif model == 3
    disp ('Model: DEXMED Hannivoort')
    disp ('Concentration = 4 mcg/mL')
    disp (['CpT: ' num2str(target) ' ng/mL'])
end
disp (' ')
disp ('-=-=-=-=-=-=-=Infusion Regime=-=-=-=-=-=-=-')
disp (' ')
disp (['Bolus ' num2str(round(bolus)) ' mcg'])
disp (' ')
disp ('Time (min) ----------- Infusion Rate (mL/hr)')
disp (['0          -----------   ' num2str(round(infn(899)/drugconc))])


oldinf = 0;
for t = 900:900:tmax-900
    if infn(t)/drugconc == oldinf % this is a simple way to clump infusion rates
    else
        oldinf = infn(t)/drugconc;
        if t/60<100
            disp ([num2str(t/60) '         -----------   ' num2str(infn(t)/drugconc)])
        elseif t/60<1000
            disp ([num2str(t/60) '        -----------   ' num2str(infn(t)/drugconc)])
        else
            disp ([num2str(t/60) '     -----------   ' num2str(infn(t)/drugconc)])
        end

    end
end
disp ([num2str(tmax/60) '        -----------   ' num2str(infn(tmax)/drugconc)])
disp (' ')
disp ('-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-')
if ttpe > 0
    plot (V(:,1), V(:,6));
else
    plot (V(:,1), V(:,5));
end