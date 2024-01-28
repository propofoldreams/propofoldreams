function [ABW, Abol] = augmarsh (ptdata)
%this function calculates the EleMarsh ABW and bolus for given pt data
%input assuming CeT = 4

%some presets here
maxpumprate = 1200; %mL/hr
propofolconc = 10; %mg/mL
maxinfnrate = maxpumprate * propofolconc; %mg/hr
CpT = 4; %target
tmax = 60*60*3; %how long simulation runs in seconds

%set up variables
Tmat = [(1:1:299) (300:5:3600) (3610:10:tmax)]';

%extract pt data
age = ptdata(1);
weight = ptdata(2);
height = ptdata(3);
sex = ptdata(4);
bmi = weight/(height/100)^2;

%marsh constant parameters
k10m = 0.119;
k12m = 0.112;
k21m = 0.055;
k13m = 0.042;
k31m = 0.0033;

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

%store the V and k's into a vector
Vmat = [V1 ; V2 ; V3];
kmat = [k10 ; k12 ; k21 ; k13 ; k31 ; ke0];
kmatMarsh = [k10m ; k12m ; k21m ; k13m ; k31m ; 0];

%%Build the gold standard Eleveld model
[peakconc, ttpe] = peaking(Vmat,kmat);
bolus = CpT/peakconc*10;
[Vgold, ~] = tci(CpT,bolus,Tmat,Vmat,kmat,maxinfnrate,ttpe); %Eleveld TCI algo


%Marsh ABW guess using PDreams linear regression algo (updated 29/12/23)
if sex == 0     %female
    PDABWguess = round(17.5+0.9912*weight-0.001305*weight^2+1.528e-6*weight^3+1.006e-4*height^2-3.690e-4*weight*bmi-0.2682*age+1.560e-3*age^2-0.003543*age*weight+2.322e-6*age*weight^2+0.001080*age*bmi-0.07786*bmi);
else     %male
    PDABWguess = round(16.07+0.9376*weight-0.001383*weight^2+1.684e-6*weight^3+1.292e-4*height^2-3.801e-4*weight*bmi-0.2617*age+1.614e-3*age^2-0.003841*age*weight+3.927e-6*age*weight^2+0.001340*age*bmi-0.09995*bmi);
end

%Marsh Abol guess using PDreams linear regression algo (updated 29/12/23)
if sex == 0     %female
    PDABolguess = round(38.01+3.096*weight-2.187e-3*weight^2+4.676e-6*weight^3-0.09256*height+4.097e-4*height^2-7.581e-4*weight*bmi-0.6293*age+0.005249*age^2-0.01542*age*weight-2.887e-5*age*weight^2+3.178e-7*age^2*weight^2+0.002109*age*bmi-0.1769*bmi);
else    %male
    PDABolguess = round(42.34+2.947*weight-1.996e-3*weight^2+4.323e-6*weight^3-0.09979*height+4.667e-4*height^2-8.554e-4*weight*bmi-0.6839*age+0.005336*age^2-0.01454*age*weight-2.864e-5*age*weight^2+2.875e-7*age^2*weight^2+0.002405*age*bmi-0.2078*bmi);
end

%% STARTING OPTIMISATION ALGORITHM

minEbolus = round(PDABolguess * 0.75);
maxEbolus = round(PDABolguess * 1.25);

%calc size of guess matrix
gsize = floor((maxEbolus - minEbolus)/2)*4;
GuessMatrix = zeros(gsize,4);
guessindex = 0;

for WGuess = PDABWguess-round(PDABWguess*0.1,0):1:PDABWguess+round(PDABWguess*0.1,0)
    %marsh parameters
    V1m = 0.228 * WGuess;
    V2m = 0.463 * WGuess;
    V3m = 2.893 * WGuess;
    VmatMarsh = [V1m ; V2m ; V3m];
    for BGuess = minEbolus:2:maxEbolus
        guessindex = guessindex + 1;    
        bolusMarsh = BGuess;

        %run the TCI algo
        [~, infnMarsh] = tci(CpT,bolusMarsh,Tmat,VmatMarsh,kmatMarsh,maxinfnrate); %generate marsh infn regime
        Vguess = pkmodel(infnMarsh,Tmat,Vmat,kmat); %plug the marsh infn regime into Eleveld
        Cediff = Vguess(:,6) - Vgold(:,6);
        PE = Cediff(5:end) ./ Vgold((5:1:end),6);
        sumsqr = sqrt(sum(Cediff.^2)/length(Cediff)); %sumsqr calc
        maxape = 100*max(abs(PE)); %maxAPE calc
        %record the result in the GuessMatrix in order
        GuessMatrix(guessindex,1) = WGuess;
        GuessMatrix(guessindex,2) = BGuess;
        GuessMatrix(guessindex,3) = sumsqr;
        GuessMatrix(guessindex,4) = maxape;
    end
end
minCe = min(GuessMatrix(:,3));
minloc = find (GuessMatrix(:,3)==minCe);
if ~isempty(minloc)
    %disp ('Multiple minima found')
    GuessMatrixNEW = GuessMatrix(minloc, :);
    [~,minmaxapeloc] = min(GuessMatrixNEW(:,4)); %if multiple min sqr error then pick the one with lowest maxAPE
    ABW = GuessMatrixNEW (minmaxapeloc,1);
    Abol = GuessMatrixNEW (minmaxapeloc,2);
else
    ABW = GuessMatrix (minloc,1);
    Abol = GuessMatrix (minloc,2);
end

%re-calculate the best bolus regime
%V1m = 0.228 * ABW;
%V2m = 0.463 * ABW;
%V3m = 2.893 * ABW;
%VmatMarsh = [V1m ; V2m ; V3m];
%[~, infnMarsh] = tci(CpT,Abol,Tmat,VmatMarsh,kmatMarsh,maxinfnrate);
%Vbest = pkmodel(infnMarsh,Tmat,Vmat,kmat);



