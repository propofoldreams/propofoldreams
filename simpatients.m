%% Let's sim pts!

ptindex = 0;
maxinfnrate = 12000;
outputmat = zeros(2768,9);

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
toc