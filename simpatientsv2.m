%% Let's sim pts!

ptindex = 0;
propolconc = 10;
maxpumprate = 1200;
maxinfnrate = maxpumprate*propolconc;
outputmat = zeros(2768,9);

clf %We gonna do some plots!
hold("on")

tic

% Ladies first
sex = 0;
for weight = 40:10:200
    for bmi = 13:5:83
        height = round(sqrt(weight/bmi)*100);
        if height < 100 %check if the height is within limits and reasonable
            %skip this
        else
            if height > 210
                %skip this
            else
                %passed the weight and height criteria! Let's build the
                %matrix
                for age = 20:10:90
                    ptindex = ptindex + 1;
                    disp (['simulating patient ' num2str(ptindex)])
                    ptdata = [age ; weight ; height ; sex];
                    [ABW, Abol] = augmarsh(ptdata); %optimise the Marsh using PDreams algo
                    errorfn = regimetest(ptdata,ABW,Abol,maxinfnrate);
                    outputmat(ptindex, 1) = age;
                    outputmat(ptindex, 2) = weight;
                    outputmat(ptindex, 3) = height;
                    outputmat(ptindex, 4) = sex;
                    outputmat(ptindex, 5) = ABW;
                    outputmat(ptindex, 6) = Abol;
                    outputmat(ptindex, 7) = errorfn(1); %mdpe
                    outputmat(ptindex, 8) = errorfn(2); %mdape
                    outputmat(ptindex, 9) = errorfn(3); %maxape
                end
            end
        end
    end
end


% Now do the boys
disp ('--------------------------- Simulating XYs')
sex = 1;

for weight = 40:10:200
    for bmi = 13:5:83
        height = round(sqrt(weight/bmi)*100);
        if height < 100 %check if the height is within limits and reasonable
            %skip this
        else
            if height > 210
                %skip this
            else
                %passed the weight and height criteria! Let's build the
                %matrix
                for age = 20:10:90
                    ptindex = ptindex + 1;
                    disp (['simulating patient ' num2str(ptindex)])
                    ptdata = [age ; weight ; height ; sex];
                    [ABW, Abol] = augmarsh(ptdata); %optimise the Marsh using PDreams algo
                    errorfn = regimetest(ptdata,ABW,Abol,maxinfnrate);
                    outputmat(ptindex, 1) = age;
                    outputmat(ptindex, 2) = weight;
                    outputmat(ptindex, 3) = height;
                    outputmat(ptindex, 4) = sex;
                    outputmat(ptindex, 5) = ABW;
                    outputmat(ptindex, 6) = Abol;
                    outputmat(ptindex, 7) = errorfn(1); %mdpe
                    outputmat(ptindex, 8) = errorfn(2); %mdape
                    outputmat(ptindex, 9) = errorfn(3); %maxape
                end
            end
        end
    end
end

%% store result into a csv
%csvwrite("ptsim.csv",outputmat);

%% display results
disp ('plotting INDEX patients')
ptdata = [30 ; 200 ; 158 ; 0];
[ABW, Abol] = augmarsh(ptdata); %optimise the Marsh using PDreams algo
[~, Vgold, Vmarsh] = regimetest(ptdata,ABW,Abol,maxinfnrate);
plot (Vgold(:,1)/60, Vgold(:,6),'color', [0 0 0], 'linewidth', 2)
plot (Vmarsh(:,1)/60, Vmarsh(:,6),'r--', 'linewidth', 1.5)

ptdata = [30 ; 40 ; 163 ; 0];
[ABW, Abol] = augmarsh(ptdata); %optimise the Marsh using PDreams algo
[~, ~, Vmarsh] = regimetest(ptdata,ABW,Abol,maxinfnrate);
plot (Vmarsh(:,1)/60, Vmarsh(:,6),'g--', 'linewidth', 1.5)

ptdata = [90 ; 150 ; 158 ; 0];
[ABW, Abol] = augmarsh(ptdata); %optimise the Marsh using PDreams algo
[~, ~, Vmarsh] = regimetest(ptdata,ABW,Abol,maxinfnrate);
plot (Vmarsh(:,1)/60, Vmarsh(:,6),'b--', 'linewidth', 1.5)

set(gcf,'Position',[10 10 1618 1000])
set(gca,'FontWeight','bold') 
set(gca,'FontSize',18)
xticks([0 60 120 180 240])
yticks([0 1 2 3 4])


toc