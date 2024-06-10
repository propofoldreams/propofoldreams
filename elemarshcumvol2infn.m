function [V, infn] = elemarshcumvol2infn(infndata, ptdata, out)
% This algorithm back calculates the infusion regime that would have resulted in the
% infndata vector which contains:
% column 1 = time
% column 2 = Cp
% column 3 = cumulative volume


% this is the test dataset that pt entered
% ptdata = [52 150 187 1];
% infndata = [526	5.6	363
%     833	4.3	445
%     1160	3.9	528
%     2270	3.4	750
%     2974	3.3	879
%     4396	3.3	1140
%     5011	3.3	1250
%     7147	4	1680
%     8583	4	2000
%     10507	4	2420
%     10817	3.5	2470
%     10955	2.7	2470
%     11074	2.4	2470
%     11149	2.2	2470
%     11192	2.1	2470
%     11221	2.1	2470];

siz = size(infndata,1);
maxinfnrate = 12000;

if nargin == 2
    out = 0;
end

%% Column detection
% column 1 is always time
% if column 2 contains cumulative dose and col 3 contains Cp then do a flip

if sum(infndata(:,2)>20) > 1
    infndata = [infndata(:, 1) infndata(:, 3) infndata(:, 2)];
end

%% Extract pt data
age = ptdata(1);
weight = ptdata(2);
height = ptdata(3);
sex = ptdata(4);
bmi = weight/(height/100)^2;


%% EleMarsh constant parameters
k10m = 0.119;
k12m = 0.112;
k21m = 0.055;
k13m = 0.042;
k31m = 0.0033;

%% EleMarsh ABW guess using PDreams linear regression algo (updated 29/12/23)
if sex == 0     %female
    EleMarshABW = round(17.5+0.9912*weight-0.001305*weight^2+1.528e-6*weight^3+1.006e-4*height^2-3.690e-4*weight*bmi-0.2682*age+1.560e-3*age^2-0.003543*age*weight+2.322e-6*age*weight^2+0.001080*age*bmi-0.07786*bmi);
    PDABolguess = round(41.31+3.063*weight-2.312e-3*weight^2+6.172e-6*weight^3-0.1026*height+4.375e-4*height^2-5.997e-4*weight*bmi-0.5831*age+0.004267*age^2-0.01399*age*weight-3.716e-5*age*weight^2+3.345e-7*age^2*weight^2+0.001912*age*bmi-0.1885*bmi);
else     %male
    EleMarshABW = round(16.07+0.9376*weight-0.001383*weight^2+1.684e-6*weight^3+1.292e-4*height^2-3.801e-4*weight*bmi-0.2617*age+1.614e-3*age^2-0.003841*age*weight+3.927e-6*age*weight^2+0.001340*age*bmi-0.09995*bmi);
    PDABolguess = round(47.92+2.983*weight-2.339e-3*weight^2+6.439e-6*weight^3-0.1693*height+6.393e-4*height^2-5.025e-4*weight*bmi-0.5454*age+0.003780*age^2-0.01376*age*weight-4.149e-5*age*weight^2+3.661e-7*age^2*weight^2+0.002259*age*bmi-0.2682*bmi);
end

V1m = 0.228 * EleMarshABW;
V2m = 0.463 * EleMarshABW;
V3m = 2.893 * EleMarshABW;

kmatMarsh = [k10m ; k12m ; k21m ; k13m ; k31m ; 0];
VmatMarsh = [V1m ; V2m ; V3m];

%% Set up time matrix

% let's work it out for first time step
% we know the cum dose = cumvol (i.e. 363 mg)
% we know the time of interval is 526 sec
% we know the interval ends with Cp = 5.6

% work out 2 parameters such that first t seconds of time interval is
% maintained at CpT of x and then remaining of time interval (526 - t) maintained at
% CpT of 5.6 such that the total infused volume is 363 mg
Velemarsh = zeros(infndata(1,1), 5);
Velemarsh(:,1) = (1:1:infndata(1,1))';

indx = 0; % set up the indx keeping vector
infnstop = []; % set up infusion stop finding vector

for cpguess = infndata(1,2):0.2:1.7*infndata(1,2)
    % bolus = PDABolguess / 4 * cpguess;
    % bolustime = round(PDABolguess / maxinfnrate * 3600,0);

    for guesst = 1:1:infndata(1,1)-1
        indx = indx + 1;
        %construct the cpt profile and see what's the cumvol at end
        cptprofile = zeros (infndata(1,1), 2);
        cptprofile(:,1) = (1:1:size(cptprofile, 1))'; %create time steps in 1 sec intervals
        cptprofile(1:1:guesst,2) = cpguess;
        cptprofile(guesst+1:1:end,2) = infndata(1,2);

        %% run the calc
        [Velemarsh, infn] = cpt2tci (cptprofile,VmatMarsh, kmatMarsh, maxinfnrate);

        % calculate how well is the fit
        cvol = sum(infn)/3600;
        guessmat(indx,1) = guesst;
        guessmat(indx,2) = cpguess;
        guessmat(indx,3) = (cvol - infndata(1,3))^2;
    end
end

% we found the best fit for 1st time step
loc = find(guessmat(:,3)==min(guessmat(:,3)));
bestfit = guessmat(loc,:);
%disp (bestfit)

% let's re-produce this!
bolus = PDABolguess / 4 * bestfit(2);
bolustime = round(PDABolguess / maxinfnrate * 3600,0);
cptprofile = zeros (infndata(1,1), 2);
cptprofile(:,1) = (1:1:size(cptprofile, 1))'; %create time steps in 1 sec intervals
cptprofile(1:1:bolustime,2) = 99;
cptprofile(bolustime+1:1:bestfit(1),2) = bestfit(2);
cptprofile(bestfit(1)+1:1:end,2) = infndata(1,2);

[V, infn] = cpt2tci (cptprofile,VmatMarsh, kmatMarsh, maxinfnrate);


%% Loop over remaining time steps
for i = 2:1:siz
    Vstate = V(end,:);
    ttarget = infndata(i,1) - infndata(i-1,1)+1;
    voltarget = infndata(i,3)-infndata(i-1,3);
    if voltarget == 0
        infnstop = [infnstop ; i-1];
        cpttarget = 0;
    else
        cpttarget = Vstate(5);
    end


    %construct the cpt profile and see what's the cumvol at end
    cptprofile = zeros (ttarget, 2);
    cptprofile(:,1) = (infndata(i-1,1):1:infndata(i,1))'; %create time steps in 1 sec intervals
    cptprofile(:,2) = cpttarget;

    %% run the calc
   
    [Velemarsh, infnelemarsh] = cpt2tci (cptprofile,VmatMarsh, kmatMarsh, maxinfnrate,Vstate);

    if voltarget == 0
        infnelemarsh = infnelemarsh * 0;
    else
        infnelemarsh = infnelemarsh / sum(infnelemarsh) * voltarget * 3600;
    end

    % generate the V matrix and integrate it into the main matrix
    Velemarsh = pkmodel(infnelemarsh,Velemarsh(:,1), VmatMarsh, kmatMarsh,Vstate);
    Velemarsh(1,:) = [];
    Velemarsh(:,6) = [];
    infnelemarsh(1) = [];
    V = [V ; Velemarsh];
    infn = [infn ; infnelemarsh];
end

%% Figure out when the pump was stopped
% let's parse through the infusion stop vector to find where the infusion
% ended

if length(infnstop) == 0
    disp ('Infusion has not stopped')
else
infnstoptimeindx = infnstop(end); % this is the baseline
for t = 0:1:length(infnstop)-2
    if infnstop(end-t) - infnstop(end-t-1) >  1
        break
    else
        infnstoptimeindx = infnstop(end-t-1);
    end
end

infnstoptimemax = infndata(infnstoptimeindx,1);
infnstoptimemin = infndata(infnstoptimeindx-1,1);

indx = 0;
    Vstate = V(infnstoptimemin,:);
    cpttarget = Vstate(5);
    ttarget = infndata(end,1) - infnstoptimemin + 1;
    voltarget = infndata(infnstoptimeindx,3) - infndata(infnstoptimeindx-1,3);

for t = infnstoptimemin:1:infnstoptimemax
    indx = indx + 1;
    % set up what would happen if infusion were to terminate NOW
    stoptime = t - infnstoptimemin + 1;

    %construct the cpt profile and see what's the cumvol at end
    cptprofile = zeros (ttarget, 2);
    cptprofile(:,1) = (infnstoptimemin:1:infndata(end,1))'; %create time steps in 1 sec intervals
    cptprofile(1:1:stoptime,2) = cpttarget;
    cptprofile(stoptime:1:end,2) = 0;

    %% run the calc
    [Velemarsh, infnelemarsh] = cpt2tci (cptprofile,VmatMarsh, kmatMarsh, maxinfnrate, Vstate);

    infnelemarsh = infnelemarsh / sum(infnelemarsh) * voltarget * 3600;

    % generate the V matrix and integrate it into the main matrix
    Velemarsh = pkmodel(infnelemarsh,Velemarsh(:,1), VmatMarsh, kmatMarsh,Vstate);
    Velemarsh(1,:) = [];
    Velemarsh(:,6) = [];
    infnelemarsh(1) = [];

    % let's calculate the sqr errors
    sqrdiff(indx, 1) = t;
    sqrdiff(indx, 2) = 0;
    for k = infnstoptimeindx:1:size(infndata,1)
        ttest = infndata(k,1);
        sqrdiff(indx, 2) = sqrdiff(indx, 2) + (Velemarsh(find(Velemarsh(:,1)==ttest),5) - infndata(k,2))^2;
    end
end

loc = find(sqrdiff(:,2)==min(sqrdiff(:,2))); %this is the optimum stop time
disp (['Time infn stopped = ' num2str(sqrdiff(loc,1))])
%stoptime = sqrdiff(loc,1); %this is the time stamp where infn stops
stoptime = sqrdiff(loc,1) - infnstoptimemin + 1; % convert time stamp into duration

% reproduce the matrix
    Vstate = V(infnstoptimemin,:);
    cpttarget = Vstate(5);
    ttarget = infndata(end,1) - infnstoptimemin + 1;
    voltarget = infndata(infnstoptimeindx,3) - infndata(infnstoptimeindx-1,3);
    cptprofile = zeros (ttarget, 2);
    cptprofile(:,1) = (infnstoptimemin:1:infndata(end,1))'; %create time steps in 1 sec intervals
    cptprofile(1:1:stoptime,2) = cpttarget;
    cptprofile(stoptime+1:1:end,2) = 0;

    [Velemarsh, infnelemarsh] = cpt2tci (cptprofile,VmatMarsh, kmatMarsh, maxinfnrate, Vstate);

    infnelemarsh = infnelemarsh / sum(infnelemarsh) * voltarget * 3600;
    Velemarsh = pkmodel(infnelemarsh,Velemarsh(:,1), VmatMarsh, kmatMarsh,Vstate);
    Velemarsh(1,:) = [];
    Velemarsh(:,6) = [];
    infnelemarsh(1) = [];

    V(Velemarsh(1,1):1:end,:) = [];
    V = [V ; Velemarsh];
    infn(Velemarsh(1,1):1:end) = [];
    infn = [infn ; infnelemarsh];
end

if out == 1
% test plot
plot (V(:,1), V(:,5), 'b-', infndata(:,1), infndata(:,2), 'ko')
writematrix([(1:1:length(infn))' infn], 'infusion.csv')
end