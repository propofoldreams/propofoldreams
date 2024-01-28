%Wake up prediction function

confirm = 0;
maxpumprate = 1200; %mL/hr

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
    target = input ('Maintenance effect site target (ng/mL) ');
    stoptime = input ('Time when pump was stoppped (sec) ');
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
    disp ('Model: PROPOFOL Eleveld')
    disp (['CeT: ' num2str(target) ' mcg/mL'])
    disp (['Pump stopped at time ' num2str(stoptime) ' seconds'])
    disp (' ')
    confirm = input ('1 = Yes, 0 = No ');
    clc
end

%% SETTING UP variables
tmax = stoptime; %how long simulation runs in seconds
conc = zeros (tmax,1);
Tmat = (1:1:tmax)';
infn = zeros(tmax,1);

%% SETTING UP MODEL PARAMETERS
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
CpT = target;

%store the V and k's into a vector
Vmat = [V1 ; V2 ; V3];
kmat = [k10 ; k12 ; k21 ; k13 ; k31 ; ke0];

%% Build gold standard model
disp ('Performing PROPOFOL Eleveld CeT calculations')
disp (' ')
[peakconc, ttpe] = peaking(Vmat,kmat);
bolus = CpT/peakconc*10;
[Vgold, Infngold] = tci(CpT,bolus,Tmat,Vmat,kmat,maxinfnrate,ttpe); %build the gold standard

%% Decrement calculations
Vstate = Vgold(end,:);
Tdecmat = (1:1:3600)';
[Vdec, ~] = titrate(0,Vstate,Tdecmat,Vmat,kmat,maxinfnrate,1);
Vdec(1,:) = [];
V = [Vgold ; Vdec];

%% BIS calculations
basebis = 93;
ceadjust = 1; %for future use
stimadjust = 0; %for future use
bismat = zeros(length(V(:,6)),1);
bismatmod = zeros(length(V(:,6)),1);


bismod = input ('BIS value for adjustment ');
cemod = input ('Ce value for adjustment ');

% modify Ce50 based on input for given pt
if cemod > ce50
    ce50new = cemod/((basebis/bismod-1)^(1/1.47));
else
    ce50new = cemod/((basebis/bismod-1)^(1/1.89));
end

ce50 = 3.08*exp(-0.00635*(age-35))*ceadjust;
disp (['Ce50 = ' num2str(ce50)])
disp (['Modified Ce50 = ' num2str(ce50new)])
disp (' ')


for n = 1:1:length(V(:,6))
    if V(n,6) < stimadjust %stimadjust is an offset/shift factor, normally set it to 0
        bismat(n) = basebis;
    elseif V(n,6)-stimadjust > ce50
        bismat(n) = basebis*((ce50^1.47))/(ce50^1.47+(V(n,6)-stimadjust)^1.47);
    else
        bismat(n) = basebis*((ce50^1.89)/(ce50^1.89+(V(n,6)-stimadjust)^1.89));
    end
    if V(n,6) < stimadjust %stimadjust is an offset/shift factor, normally set it to 0
        bismatmod(n) = basebis;
    elseif V(n,6)-stimadjust > ce50
        bismatmod(n) = basebis*((ce50new^1.47))/(ce50new^1.47+(V(n,6)-stimadjust)^1.47);
    else
        bismatmod(n) = basebis*((ce50new^1.89)/(ce50new^1.89+(V(n,6)-stimadjust)^1.89));
    end
end

%% PLOT results
figure
yyaxis left
plot (V(:,1), V(:,6))
yyaxis right
plot (V(:,1), bismat, V(:,1), bismatmod)
