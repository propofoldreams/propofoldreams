function [V, infn] = titrate(CpT,Vstate,Tmat,Vmat,kmat,maxinfnrate,effect)
%CpT is the new target compartment concentration
%Vstate is the final state of compartments before change
%Tmat is time matrix
%Vmat and kmat define the 3 compartment model
%maxinfnrate is the max infusion rate of pump
%effect = 1 for effect site target; = 0 for plasma target

%baseline parameters, setting everything up
V1 = Vmat(1);
V2 = Vmat(2);
V3 = Vmat(3);
k10 = kmat(1);
k12 = kmat(2);
k21 = kmat(3);
k13 = kmat(4);
k31 = kmat(5);
ke0 = kmat(6);

%figure out the time steps and set outputs up
siz = size(Tmat);

infn = zeros(siz(1),1);
V = zeros(siz(1),6);
V(:,1) = Tmat(:,1) + Vstate(1) - Tmat(1);
V(1,2:1:6) = Vstate(2:1:6);

%let's figure out what is the current status
if nargin == 6 %default is effect site targeting
    effect = 1;
end

if effect == 0
    conc = Vstate(5);
else
    %effect site target
    conc = Vstate(6);
    peakconc = peakinginfuse(Vmat,kmat,maxinfnrate);
end

%% MAIN BODY OF ALGO

% Let's write it for EFFECT SITE targeting first
if effect == 1
if CpT > conc
    %Uptitration - need to bolus + infuse

    %This part of the algorithm optimises the bolus infusion to get to CeT
    %as fast as possible iteratively
    for step = 1:1:200
        bolustime = (CpT-conc)/peakconc*(100+step-1) + Vstate(1); %bolus time initial GUESS
        for t = 2:1:300
            tstep = V(t,1)-V(t-1,1);
            dV1 = (k21*V(t-1,3)+k31*V(t-1,4)-V(t-1,2)*(k10+k12+k13))*tstep/60; %delta V1 compartment from redistribution
            if V(t,1) <= bolustime %during the bolus period
                infn(t) = maxinfnrate;
            elseif V(t,1) <= bolustime + 1 %borderline period
                sugginfn = 3600 * (CpT*V1-dV1-V(t-1,2)) / tstep;
                if sugginfn < 0
                    infn(t) = maxinfnrate * (bolustime + 1 - V(t,1));
                else
                    infn(t) = maxinfnrate * (bolustime + 1 - V(t,1)) + sugginfn * (V(t,1) - bolustime);
                end
            else %out of the bolus period
                sugginfn = 3600 * (CpT*V1-dV1-V(t-1,2)) / tstep;
                if sugginfn < 0
                    infn(t) = 0; %delta V1 compartment from infusion
                else
                    infn(t) = sugginfn;
                    if effect == 1 %i.e. we are in effect site targeting
                        if V(t,6) > CpT %continue pausing infn until Ce reaches CeT
                            infn(t) = 0;
                        end
                    end
                end
            end
            V(t,2) = V(t-1,2) + dV1 + infn(t)*tstep/3600; %iterative calc V1 drug
            V(t,3) = V(t-1,3) + (k12*V(t-1,2)-k21*V(t-1,3))*tstep/60; %iterative V2
            V(t,4) = V(t-1,4) + (k13*V(t-1,2)-k31*V(t-1,4))*tstep/60; %iterative V3
            V(t,5) = V(t,2)/V1;
            V(t,6) = V(t-1,6) + (V(t-1,5)-V(t-1,6))*ke0/60;
        end
        VAUCtest(step, :) = [bolustime sum((CpT - V(1:1:300,6)).^2)];
        if max(V(1:1:300,6)) > CpT
            break
        end
    end
    [~, d] = min(VAUCtest);
    bolustime = VAUCtest(d(2), 1);

    %reset everything
    infn = zeros(siz(1),1);
    V = zeros(siz(1),6);
    V(:,1) = Tmat(:,1) + Vstate(1) - Tmat(1);
    V(1,2:1:6) = Vstate(2:1:6);

    % Now that we figured out the best bolus, start the true infusion regime
    for t = 2:1:siz(1)
        tstep = V(t,1)-V(t-1,1);
        dV1 = (k21*V(t-1,3)+k31*V(t-1,4)-V(t-1,2)*(k10+k12+k13))*tstep/60; %delta V1 compartment from redistribution
        if V(t,1) <= bolustime %during the bolus period
            infn(t) = maxinfnrate;
        elseif V(t,1) <= bolustime + 1 %borderline period
            sugginfn = 3600 * (CpT*V1-dV1-V(t-1,2)) / tstep;
            if sugginfn < 0
                infn(t) = maxinfnrate * (bolustime + 1 - V(t,1));
            else
                infn(t) = maxinfnrate * (bolustime + 1 - V(t,1)) + sugginfn * (V(t,1) - bolustime);
            end
        else %out of the bolus period
            sugginfn = 3600 * (CpT*V1-dV1-V(t-1,2)) / tstep;
            if sugginfn < 0
                infn(t) = 0; %delta V1 compartment from infusion
            else
                infn(t) = sugginfn;
                if effect == 1 %i.e. we are in effect site targeting
                    if V(t,6) > CpT %continue pausing infn until Ce reaches CeT
                        infn(t) = 0;
                    end
                end
            end
        end
        V(t,2) = V(t-1,2) + dV1 + infn(t)*tstep/3600; %iterative calc V1 drug
        V(t,3) = V(t-1,3) + (k12*V(t-1,2)-k21*V(t-1,3))*tstep/60; %iterative V2
        V(t,4) = V(t-1,4) + (k13*V(t-1,2)-k31*V(t-1,4))*tstep/60; %iterative V3
        V(t,5) = V(t,2)/V1;
        V(t,6) = V(t-1,6) + (V(t-1,5)-V(t-1,6))*ke0/60;
    end

else
    %Downtitration or continue current
    tstep = Tmat(2) - Tmat(1);
    for t = 2:1:siz(1)
        dV1 = (k21*V(t-1,3)+k31*V(t-1,4)-V(t-1,2)*(k10+k12+k13))*tstep/60;
        sugginfn = 3600 * (CpT*V1-dV1-V(t-1,2)) / tstep;
        if sugginfn < 0
            infn(t) = 0; %pause infusion
        else
            infn(t) = sugginfn;
            if effect == 1 %i.e. we are in effect site targeting
                if V(t,6) > CpT %continue pausing infn until Ce reaches CeT
                    infn(t) = 0;
                end
            end
        end
        V(t,2) = V(t-1,2) + dV1 + infn(t)*tstep/3600; %iterative calc V1 drug
        V(t,3) = V(t-1,3) + (k12*V(t-1,2)-k21*V(t-1,3))*tstep/60; %iterative V2
        V(t,4) = V(t-1,4) + (k13*V(t-1,2)-k31*V(t-1,4))*tstep/60; %iterative V3
        V(t,5) = V(t,2)/V1;
        V(t,6) = V(t-1,6) + (V(t-1,5)-V(t-1,6))*ke0/60;
    end
end
end

% Now write it for PLASMA targeting
if effect == 0
    if CpT > conc
        %Uptitration - need to bolus + infuse

        bolustime = (CpT-conc)*V1/maxinfnrate + Vstate(1);
        for t = 2:1:300
            tstep = V(t,1)-V(t-1,1);
            dV1 = (k21*V(t-1,3)+k31*V(t-1,4)-V(t-1,2)*(k10+k12+k13))*tstep/60; %delta V1 compartment from redistribution
            if V(t,1) <= bolustime %during the bolus period
                infn(t) = maxinfnrate;
            elseif V(t,1) <= bolustime + 1 %borderline period
                sugginfn = 3600 * (CpT*V1-dV1-V(t-1,2)) / tstep;
                if sugginfn < 0
                    infn(t) = maxinfnrate * (bolustime + 1 - V(t,1));
                else
                    infn(t) = maxinfnrate * (bolustime + 1 - V(t,1)) + sugginfn * (V(t,1) - bolustime);
                end
            else %out of the bolus period
                % OK we are out of the bolus period, let's check what is
                % the final Cp reached here
                V(t,2) = V(t-1,2) + dV1 + infn(t)*tstep/3600;
                cptfinal = V(t,2)/V1; % this is the Cp that bolustime would reach - our first loop is done!
                break
            end
            V(t,2) = V(t-1,2) + dV1 + infn(t)*tstep/3600; %iterative calc V1 drug
            V(t,3) = V(t-1,3) + (k12*V(t-1,2)-k21*V(t-1,3))*tstep/60; %iterative V2
            V(t,4) = V(t-1,4) + (k13*V(t-1,2)-k31*V(t-1,4))*tstep/60; %iterative V3
            V(t,5) = V(t,2)/V1;
        end

        %reset everything
        infn = zeros(siz(1),1);
        V = zeros(siz(1),6);
        V(:,1) = Tmat(:,1) + Vstate(1) - Tmat(1);
        V(1,2:1:6) = Vstate(2:1:6);
        
        % scale the guesstime bolustime by cptfinal to figure out the best bolus time
        bolustime = (CpT-conc)*V1/maxinfnrate * CpT/cptfinal + Vstate(1);

        % Now that we figured out the best bolus, start the true infusion regime
        for t = 2:1:siz(1)
            dV1 = (k21*V(t-1,3)+k31*V(t-1,4)-V(t-1,2)*(k10+k12+k13))*tstep/60; %delta V1 compartment from redistribution
            if V(t,1) <= bolustime %during the bolus period
                infn(t) = maxinfnrate;
            elseif V(t,1) <= bolustime + 1 %borderline period
                sugginfn = 3600 * (CpT*V1-dV1-V(t-1,2)) / tstep;
                if sugginfn < 0
                    infn(t) = maxinfnrate * (bolustime + 1 - V(t,1));
                else
                    infn(t) = maxinfnrate * (bolustime + 1 - V(t,1)) + sugginfn * (V(t,1) - bolustime);
                end
            else %out of the bolus period
                sugginfn = 3600 * (CpT*V1-dV1-V(t-1,2)) / tstep;
                if sugginfn < 0
                    infn(t) = 0; %theoretically, we never get to here with plasma targeting
                else
                    infn(t) = sugginfn;
                end
            end
            V(t,2) = V(t-1,2) + dV1 + infn(t)*tstep/3600; %iterative calc V1 drug
            V(t,3) = V(t-1,3) + (k12*V(t-1,2)-k21*V(t-1,3))*tstep/60; %iterative V2
            V(t,4) = V(t-1,4) + (k13*V(t-1,2)-k31*V(t-1,4))*tstep/60; %iterative V3
            V(t,5) = V(t,2)/V1;
        end

    else
        %Downtitration or continue current
        tstep = Tmat(2) - Tmat(1);
        for t = 2:1:siz(1)
            dV1 = (k21*V(t-1,3)+k31*V(t-1,4)-V(t-1,2)*(k10+k12+k13))*tstep/60;
            sugginfn = 3600 * (CpT*V1-dV1-V(t-1,2)) / tstep;
            if sugginfn < 0
                infn(t) = 0; %pause infusion
            else
                infn(t) = sugginfn;
            end
            V(t,2) = V(t-1,2) + dV1 + infn(t)*tstep/3600; %iterative calc V1 drug
            V(t,3) = V(t-1,3) + (k12*V(t-1,2)-k21*V(t-1,3))*tstep/60; %iterative V2
            V(t,4) = V(t-1,4) + (k13*V(t-1,2)-k31*V(t-1,4))*tstep/60; %iterative V3
            V(t,5) = V(t,2)/V1;
        end
    end
end
