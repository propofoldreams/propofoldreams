%some presets here
clear
confirm = 0;
maxpumprate = 1200; %mL/hr
propofolconc = 10; %mg/mL
maxinfnrate = maxpumprate * propofolconc; %mg/hr
CpT = 4; %target
tmax = 60*60*3; %how long simulation runs in seconds

%set up variables
Tmat = [(1:1:299) (300:5:3600) (3610:10:tmax)]';

while confirm == 0
    clc
    disp ('Propofol Dreams v2.0')
    disp ('EleMarsh Model Builder')
    disp (' ')
    disp (' ')
    age = input ('Age (yr) ');
    weight = input ('Weight (kg) ');
    height = input ('Height (cm) ');
    sex = input ('Sex (0 = XX, 1 = XY) ');
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
    disp (' ')
    confirm = input ('1 = Yes, 0 = No ');
    clc
end
disp ('Calculating model parameters')
disp (' ')

%setting up the models
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

disp ('Eleveld Volumes')
disp (['V1 = ' num2str(Vmat(1)) '  V2 = ' num2str(Vmat(2)) '  V3 = ' num2str(Vmat(3))])
disp (' ')
disp ('Eleveld Rate Constants')
disp (kmat')
disp (' ')

%%Build the gold standard Eleveld model
[peakconc, ttpe] = peaking(Vmat,kmat);
bolus = CpT/peakconc*10;

[Vgold, EleveldInfn] = tci(CpT,bolus,Tmat,Vmat,kmat,maxinfnrate,ttpe); %Eleveld TCI algo

disp ('Building Eleveld Model ... done!')
disp (' ')

%Marsh ABW guess using Propofol Dreams linear regression algo
if sex == 0
    %female
    PDABWguess = round(15.24+1.033*weight-0.001552*weight^2+2.119e-6*weight^3+8.909e-5*height^2-4.423e-4*weight*bmi-0.1928*age+9.729e-4*age^2-0.003927*age*weight+1.779e-6*age*weight^2+0.001165*age*bmi-0.08306*bmi);
else
    %male
    PDABWguess = round(15.03+0.9526*weight-0.001513*weight^2+1.991e-6*weight^3+1.144e-4*height^2-4.308e-4*weight*bmi-0.2029*age+1.047e-3*age^2-0.003866*age*weight+3.305e-6*age*weight^2+0.001263*age*bmi-0.09866*bmi);
end

%Marsh Abol guess using Propofol Dreams linear regression algo
if sex == 0
    %female
    PDABolguess = round(41.31+3.063*weight-2.312e-3*weight^2+6.172e-6*weight^3-0.1026*height+4.375e-4*height^2-5.997e-4*weight*bmi-0.5831*age+0.004267*age^2-0.01399*age*weight-3.716e-5*age*weight^2+3.345e-7*age^2*weight^2+0.001912*age*bmi-0.1885*bmi);

else
    %male
    PDABolguess = round(47.92+2.983*weight-2.339e-3*weight^2+6.439e-6*weight^3-0.1693*height+6.393e-4*height^2-5.025e-4*weight*bmi-0.5454*age+0.003780*age^2-0.01376*age*weight-4.149e-5*age*weight^2+3.661e-7*age^2*weight^2+0.002259*age*bmi-0.2682*bmi);
end

%% here are the linear regression algo for 2% propofol
%
%
%
%
%if sex == 0
%    %female
%    PDABWguess = round(17.92+0.9685*weight-0.001057*weight^2+1.21e-6*weight^3+9.414e-5*height^2-4.366e-4*weight*bmi-0.2552*age+1.281e-3*age^2-0.003136*age*weight-6.911e-7*age*weight^2+0.001090*age*bmi-0.07710*bmi);
%else
%    %male
%    PDABWguess = round(16.42+0.9375*weight-0.001365*weight^2+1.838e-6*weight^3+1.256e-4*height^2-4.43e-4*weight*bmi-0.2637*age+1.557e-3*age^2-0.003772*age*weight+3.321e-6*age*weight^2+0.001355*age*bmi-0.09703*bmi);
%end

%Marsh Abol guess using Propofol Dreams linear regression algo
%if sex == 0
%    %female
%    PDABolguess = round(58.18+3.358*weight-3.991e-3*weight^2+8.628e-6*weight^3-0.3385*height+1.058e-3*height^2-0.5662*age+0.005773*age^2-0.01759*age*weight-2.245e-5*age*weight^2+3.48e-7*age^2*weight^2+0.002268*age*bmi-0.3347*bmi);
%
%else
%    %male
%    PDABolguess = round(36.37+3.024*weight-2.734e-3*weight^2+6.121e-6*weight^3+3.25e-4*height^2-1.287e-3*weight*bmi-1.023*age+0.01396*age^2-4.776e-5*age^3-0.01692*age*weight-1.859e-5*age*weight^2+2.843e-7*age^2*weight^2+0.002622*age*bmi-0.077*bmi);
%end

%% BIS Estimation
%EstBIS = agebis(age, CpT);
%disp (['Estimated BIS = ' num2str(round(EstBIS,0)) ' @ CpT = ' num2str(CpT)])
%disp (' ')

%% Eleveld output
disp (['Eleveld bolus = ' num2str(round(bolus)) ' mg'])
disp (['PDGuess Marsh weight = ' num2str(PDABWguess) ' kg'])
disp (['PDGuess Marsh bolus = ' num2str(PDABolguess) ' mg'])
disp (' ')

disp ('Starting optimisation algo ...')
tic
minEbolus = round(PDABolguess * 0.8);
maxEbolus = round(PDABolguess * 1.2);

%calc size of guess matrix
gsize = floor((maxEbolus - minEbolus)/2)*9;
disp (['Guessing ' num2str(gsize) ' combinations ...'])
GuessMatrix = zeros(gsize,6);
guessindex = 0;
for WGuess = PDABWguess-10:1:PDABWguess+10
    %marsh parameters
    V1m = 0.228 * WGuess;
    V2m = 0.463 * WGuess;
    V3m = 2.893 * WGuess;
    VmatMarsh = [V1m ; V2m ; V3m];
    for BGuess = minEbolus:2:maxEbolus
        %disp (['Guessing ... ' num2str(WGuess) ' kg and ' num2str(BGuess) ' mg bolus'])
        guessindex = guessindex + 1;    
        bolusMarsh = BGuess;

        %run the TCI algo
        [~, infnMarsh] = tci(CpT,bolusMarsh,Tmat,VmatMarsh,kmatMarsh,maxinfnrate); %generate marsh infn regime
        Vguess = pkmodel(infnMarsh,Tmat,Vmat,kmat); %plug the marsh infn regime into Eleveld
        Cediff = Vguess(:,6) - Vgold(:,6);
        PE = Cediff(5:end) ./ Vgold((5:1:end),6);
        sumsqr = sqrt(sum(Cediff.^2)/length(Cediff)); %sumsqr calc
        mdpe = 100*median(PE); %MDPE calc
        mdape = 100*median(abs(PE)); %MDAPE calc
        maxape = 100*max(abs(PE)); %maxAPE calc
        %record the result in the GuessMatrix in order
        GuessMatrix(guessindex,1) = WGuess;
        GuessMatrix(guessindex,2) = BGuess;
        GuessMatrix(guessindex,3) = sumsqr;
        GuessMatrix(guessindex,4) = mdpe;
        GuessMatrix(guessindex,5) = mdape;
        GuessMatrix(guessindex,6) = maxape;
    end
end

%now let's find the best combo

minCe = min(GuessMatrix(:,3));
minloc = find (GuessMatrix(:,3)==minCe);
if ~isempty(minloc)
    disp (' ')
    disp ('Multiple minima found')
    GuessMatrixNEW = GuessMatrix(minloc, :);
    [~,minmaxapeloc] = min(GuessMatrixNEW(:,6));
    bestWt = GuessMatrixNEW (minmaxapeloc,1);
    bestBol = GuessMatrixNEW (minmaxapeloc,2);
    minloc = minloc(minmaxapeloc);
else
    bestWt = GuessMatrix (minloc,1);
    bestBol = GuessMatrix (minloc,2);
end
toc
disp (' ')
disp ('All done!')

%re-calculate the best bolus regime
CpT = 4;

[Vgold, ~] = tci(CpT,bolus,Tmat,Vmat,kmat,maxinfnrate,ttpe);
bestBol = round(bestBol * 4 / 4);
V1m = 0.228 * bestWt;
V2m = 0.463 * bestWt;
V3m = 2.893 * bestWt;
VmatMarsh = [V1m ; V2m ; V3m];
[Vmarsh, infnMarsh] = tci(CpT,bestBol,Tmat,VmatMarsh,kmatMarsh,maxinfnrate);
Vbest = pkmodel(infnMarsh,Tmat,Vmat,kmat);
disp (' ')
disp ('-------------------------------')
    disp (' ')
    disp (['Age: ', num2str(age), ' yr'])
    disp (['Weight: ', num2str(weight), ' kg'])
    disp (['Height: ', num2str(height), ' cm'])
    disp (['BMI: ' num2str(round(bmi*10,0)/10)])
    if sex == 0
        disp ('Sex: Female')
    else
        disp ('Sex: Male')
    end
    disp (['CeT: ' num2str(round(CpT*10,0)/10)])
    disp (' ')
    disp ('-------------------------------')
    disp (' ')
    disp (['Optimal Marsh ABW = ' num2str(bestWt) ' kg'])
disp (['Optimal Bolus = ' num2str(bestBol) ' mg'])
disp (['Set Induction CpT at ' num2str(round(max(Vmarsh(:,5))*10)/10) ' mcg/mL'])
disp (' ')
disp (['MDPE = ' num2str(round(GuessMatrix(minloc,4))) '%'])
disp (['MDAPE = ' num2str(round(GuessMatrix(minloc,5))) '%'])
disp (['MaxAPE = ' num2str(round(GuessMatrix(minloc,6))) '%'])
    disp (' ')
    disp ('-------------------------------')

%plot the best fit

plot (Vgold(:,1), Vgold(:,6), 'r', Vbest(:,1), Vbest(:,6),'g--')
xlabel('time (s)')
ylabel('Ce (mcg/mL)')



