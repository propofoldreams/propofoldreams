function V = eleveldconc(infndata,ptdata,peeg)
% calculates the Eleveld Ce profile based on cumulative propofol
% volume-time profile
maxinfnrate = 12000;


%extract pt data
age = ptdata(1);
weight = ptdata(2);
height = ptdata(3);
sex = ptdata(4);
bmi = weight/(height/100)^2;

if nargin == 2
    peeg = 0;
end

%infndata matrix cleaning
if max(infndata(:,2)) < 15 %2nd column is Cp and 3rd column is cumvol, so let's swap 2nd and 3rd columns around
    infndata(:,4) = infndata(:,2);
    infndata(:,2) = [];
end

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

%set up time matrix
tmax = max(infndata(:,1));
tsteps = size(infndata,1);
Tmat = (1:1:tmax)';
V = zeros(tmax,7);
V(:,1) = Tmat(:,1);
infn = zeros(tmax,1);

%calc infn regime
for i = 1:1:tsteps
    if i == 1 %first time step
        % we need to first work out the bolus and infusion regime for
        % induction here! first build the gold standard
        cumvol = infndata(1,2); %how much total propofol infused?
        t = infndata(1,1); %how much time elapsed?
        [peakconc, ttpe] = peaking(Vmat,kmat);
        bolus = 3/peakconc*10;
        [~, Infngold] = tci(3,bolus,(1:1:t)',Vmat,kmat,maxinfnrate,ttpe); %calculate stereotypical infn regime with CeT = 3 corresponding to first time step
        infn(1:1:t) = Infngold/sum(Infngold/3600)*(cumvol); %scale the stereotypical infn regime to match that the cumvol of first timestep
    else
        cumvol = infndata(i,2) - infndata(i-1,2);
        t = infndata(i,1) - infndata(i-1,1) - 1;
        infn((infndata(i-1,1)+1):1:(infndata(i,1))) = cumvol/t*3600;
    end
end

for t = 2:1:tmax
    dV1 = (k21*V(t-1,3)+k31*V(t-1,4)-V(t-1,2)*(k10+k12+k13))/60; %delta V1 compartment from redistribution
    V(t,2) = V(t-1,2) + dV1 + infn(t)/3600; %iterative calc V1 drug
    V(t,3) = V(t-1,3) + (k12*V(t-1,2)-k21*V(t-1,3))/60; %iterative V2
    V(t,4) = V(t-1,4) + (k13*V(t-1,2)-k31*V(t-1,4))/60; %iterative V3
    V(t,5) = V(t,2)/V1;
    V(t,6) = V(t-1,6) + (V(t-1,5)-V(t-1,6))*ke0/60;
end

%% do some eBIS calculations
basebis = 93;
ce50 = 3.08*exp(-0.00635*(age-35));
for t = 2:1:tmax
    if V(t,6) > ce50
        V(t,7) = basebis*((ce50^1.47))/(ce50^1.47+V(t,6)^1.47);
    else
        V(t,7) = basebis*((ce50^1.89))/(ce50^1.89+V(t,6)^1.89);
    end
end

%% eBIS adjustment algorithm
% basic concept: user inputs BIS-Ce pairs for fitting
% Find the best combination of Ce50 and Ce50shift that define a sigmoid dose-response curve
% that best fit to the user input data using least squares

% Import all the bis data into a matrix named bisread whereby ...
k = size(infndata,1);
bisread = zeros(k,3);

%first column is time (in seconds)
bisread(1:1:k, 1) = infndata(1:1:k, 1); % 1st column is time

%second column is Ce
for i = 1:1:k
    bisread(i, 2) = V(infndata(i,1),6); % 2nd column is corresponding Eleveld Ce
end

%third column is BIS
bisread(1:1:k, 3) = peeg(1:1:k);

%% Let's do some DATA CLEANING
bisread(bisread(:,1) < 600, :) = []; % discard all the bis readings in the first 10 minutes (non-steady state)

%let's now check if the infn has ended
if sum(infndata (:, 2) == infndata(end,2)) > 1
    %yes, infn has ended and let's remove all the duplicates except one
    bisread(end-sum(infndata (:, 2) == infndata(end,2))+2:1:end, :) = [];
end

totweight = sum(bisread(:,1)); % this is for time-based scaling - i.e. give more weight to more recent BIS readings

% how many readings are there left in the bisread matrix?
k = size(bisread, 1);

% if there are no readings left after data cleaning, then just should use the ce50
% prediction from Eleveld
if k == 0
    ce50fit = ce50;
    ce50shiftfit = ce50*0.3;
else

    %% Bayesian BIS fitting search algorithm
    % aim is to find the best Ce50 and horizontal shift (Ce50shift) that match the dataset

    cemat = round((ce50*0.5)*100,0)/100:0.05:round((ce50*2)*100,0)/100;
    %ceshift = round((ce50*0.1)*100,0)/100:0.05:round((ce50)*100,0)/100;
    %cesqrdiff = zeros(length(cemat)*length(ceshift),3);
    p = 0;

    for n = 1:1:length(cemat)
        ce50test = cemat(n);
        ceshift = round((ce50test*0.3)*100,0)/100:0.05:round((ce50test)*100,0)/100; % shift by at least 30% of ce50
        for m = 1:1:length(ceshift)
            p = p + 1;
            ce50shift = ceshift(m);
            sqrdiffmat = zeros(1,k);
            for i = 1:1:k
                if bisread(i,2) - ce50shift < 0
                    sqrdiffmat(i) = basebis;
                elseif bisread(i,2) - ce50shift > ce50test
                    sqrdiffmat(i) = basebis*((ce50test^1.47))/(ce50test^1.47+(bisread(i,2)-ce50shift)^1.47);
                else
                    sqrdiffmat(i) = basebis*((ce50test^1.89))/(ce50test^1.89+(bisread(i,2)-ce50shift)^1.89);
                end
                sqrdiffmat(i) = (sqrdiffmat(i) - bisread(i,3))^2 *(bisread(i,1)/totweight); %weighted sqr difference
            end
            cesqrdiff(p,1) = ce50test;
            cesqrdiff(p,2) = ce50shift;
            cesqrdiff(p,3) = sum(sqrdiffmat);
        end
    end

    %find the ce50 with the smallest sqr difference
    minloc = cesqrdiff(:,3)==min(cesqrdiff(:,3));
    ce50fit = cesqrdiff(minloc,1);
    if length(ce50fit) > 1
        ce50fit = mean(ce50fit);
    end
    ce50shiftfit = cesqrdiff(minloc,2);
    if length(ce50shiftfit) > 1
        ce50shiftfit = mean(ce50shiftfit);
    end
end

%% debugging purposes
%  here I specify a MANUAL ce50fit and shift
%ce50fit = 2;
%ce50shiftfit = 0.6;

%% Generate the fitted predicted BIS curve for plotting
%eegfit is a vector containing the predicted BIS at each given time step
eegfit = zeros(tmax,1);
eegfit(1) = 99;
for t = 2:1:tmax
    if V(t,6)-ce50shiftfit < 0
        eegfit(t) = basebis;
    elseif V(t,6)-ce50shiftfit > ce50fit
        eegfit(t) = basebis*((ce50fit^1.47))/(ce50fit^1.47+(V(t,6)-ce50shiftfit)^1.47);
    else
        eegfit(t) = basebis*((ce50fit^1.89))/(ce50fit^1.89+(V(t,6)-ce50shiftfit)^1.89);
    end
end

%% Generate the fitted sigmoid dose-response curve for debugging purposes ONLY
%note the following code is for DEBUGGING purposes only
doseresp = zeros(401,2);
doseresp (:, 1) = 0:0.02:8;
for t = 1:1:401
    if doseresp(t,1)-ce50shiftfit < 0
        doseresp(t,2) = basebis;
    elseif doseresp(t,1)-ce50shiftfit > ce50fit
        doseresp(t,2) = basebis*((ce50fit^1.47))/(ce50fit^1.47+(doseresp(t,1)-ce50shiftfit)^1.47);
    else
        doseresp(t,2) = basebis*((ce50fit^1.89))/(ce50fit^1.89+(doseresp(t,1)-ce50shiftfit)^1.89);
    end
end
%Note the above code is for DEBUGGING purposes only

%% How good is the BIS fit?
% Calculate MDAPE for BIS fitting
predbis = zeros(size(bisread,1),1);
for i = 1:1:size(bisread,1)
    if bisread(i,2) - ce50shiftfit < 0
        predbis(i) = basebis;
    elseif bisread(i,2) - ce50shiftfit > ce50fit
        predbis(i) = basebis*((ce50fit^1.47))/(ce50fit^1.47+(bisread(i,2)-ce50shiftfit)^1.47);
    else
        predbis(i) = basebis*((ce50fit^1.89))/(ce50fit^1.89+(bisread(i,2)-ce50shiftfit)^1.89);
    end
end
mdape = median(abs((predbis - bisread(:,3))) ./ bisread(:,3))*100;

%% WAKE UP TIME prediction algorithm
% pump stop time at bisread(k,1)
% we iterate until wake up Ce (i.e. eBIS = 60) is reached

% set up the baseline parameters
wakeupce60 = (93/60-1)^(1/1.47)*ce50fit+ce50shiftfit;
wakeupce80 = (93/80-1)^(1/1.47)*ce50fit+ce50shiftfit;
t = bisread(k, 1);
Vwake(1, :) = V(t,:);

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

wakeuptime60 = t; %wake up time = how long after pump stopped before patient wakes up

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

wakeuptime80 = t; %wake up time = how long after pump stopped before patient wakes up

%% Clean up V matrix (eliminate the useless columns, i.e. of peripheral compartments and Cp, and get ready to present results
V(:,2:1:4) = [];


%% Display results
clc
disp ('Propofol Dreams - BAYESIAN Calculator')
disp ('-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=')
disp (' ')
disp (['Eleveld Baseline Ce50 = ' num2str(ce50)])
disp (' ')
if k > 1
    disp (['Best Ce50 fit = ' num2str(ce50fit)])
    disp (['Best Ce50 shift = ' num2str(ce50shiftfit)])
    disp (['MDAPE = ' num2str(round(mdape*10,0)/10) '%'])
    disp (' ')
    disp (' ')
    disp (['eBIS = 60   @   Ce = ' num2str((93/60-1)^(1/1.47)*ce50fit+ce50shiftfit) '   in ' num2str(floor(wakeuptime60/60)) ' min ' num2str(mod(wakeuptime60,60)) ' sec @ Time = ' num2str(bisread(k,1)+wakeuptime60)])
    disp (['eBIS = 80   @   Ce = ' num2str((93/80-1)^(1/1.47)*ce50fit+ce50shiftfit) '   in ' num2str(floor(wakeuptime80/60)) ' min ' num2str(mod(wakeuptime80,60)) ' sec @ Time = ' num2str(bisread(k,1)+wakeuptime80)])
    disp (['Predicted wake up in ' num2str(floor(floor((wakeuptime60*4+wakeuptime80)/5)/60)) ' min ' num2str(mod(floor((wakeuptime60*4+wakeuptime80)/5),60)) ' sec @ time = ' num2str(bisread(k,1)+floor((wakeuptime60*4+wakeuptime80)/5))]) %% note: bisread(k ,1) should be the time that PUMP STOPPED - otherwise wake up time not accurate
    if peeg(end) > 90
        disp (' ')
        disp (['Clinical wake up   @   Ce = ' num2str(V(end,3))])
        disp (['Clinical wake up time = ' num2str(V(end,1))])
    end
end
disp (' ')

%% Plotting
figure
%subplot (2,1,1); % this plot is for debug purposes only
%plot (doseresp(:,1), doseresp(:,2))
%subplot (2,1,2);

yyaxis left
plot (V(:,1)/60, V(:,3))
ylabel ('Ce')
xlabel ('Time (min)')

if k > 1
    yyaxis right
    %plot (V(:,1), V(:,4), 'r-', V(:,1), eegfit, 'g--', infndata(:,1), peeg, 'ko')
    plot (V(:,1)/60, eegfit, 'r--', infndata(:,1)/60, peeg, 'ko')
else
    yyaxis right
    plot (V(:,1)/60, V(:,4), 'r-', infndata(:,1)/60, peeg, 'ko')
end
ylabel ('Predicted BIS')
xline ((floor((wakeuptime60*4+wakeuptime80)/5)+bisread(k,1))/60, 'r.', {'wake up'})

%test case (pt 22)
% actcp = [526	5.6
% 833	4.3
% 1160	3.9
% 2270	3.4