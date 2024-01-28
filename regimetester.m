function V = regimetester(infn)
%%for a given infusion regime, specify model and pt deets, outputs the Cp
%regime
confirm = 0;
tmax = size(infn, 1); %how long simulation runs in seconds

%set up variables
Tmat = (1:1:tmax)';

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
    disp ('1 = PROPOFOL Marsh')
    disp ('2 = PROPOFOL Eleveld')
    disp ('3 = PROPOFOL Schnider')
    model = input ('Please select a model: ');

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
    if model == 1
        disp ('Model: PROPOFOL Marsh')
    elseif model == 2
        disp ('Model: PROPOFOL Eleveld')
    elseif model == 3
        disp ('Model: PROPOFOL Schnider')
    end
    disp (' ')
    confirm = input ('1 = Yes, 0 = No ');
    clc
end

%% SETTING UP MODEL PARAMETERS
%disp ('Model Parameters')
%disp (' ')

if model == 1
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
elseif model == 2
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
elseif model == 3 %PROPOFOL Schnider
    lbm = sex*(1.07*weight-148*(weight/height)^2)+(1-sex)*(1.1*weight-128*(weight/height)^2);
    V1 = 4.27;
    V2 = 18.9-0.391*(age-53);
    V3 = 238;
    k10 = 0.443+0.0107*(weight-77)-0.0159*(lbm-59)+0.0062*(height-177);
    k12 = 0.302-0.0056*(age-53);
    k21 = (1.29-0.024*(age-53))/(18.9-0.391*(age-53));
    k13 = 0.196;
    k31 = 0.0035;
    ke0 = 0.456;
end

%store the V and k's into a vector
Vmat = [V1 ; V2 ; V3];
kmat = [k10 ; k12 ; k21 ; k13 ; k31 ; ke0];

V = pkmodel(infn, Tmat, Vmat, kmat);