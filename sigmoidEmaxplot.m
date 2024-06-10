function sigmoidEmaxplot
ce50 = input ('Ce50 = ');
ce50shift = input ('Ce50 shift =');

ce = 0:0.01:6;
basebis = 95;
predbis = zeros(length(ce));
noshiftpredbis = predbis;

for i = 1:1:length(ce)
    if ce(i) - ce50shift <= 0
        predbis(i) = basebis;
    elseif ce(i) - ce50shift > ce50
        predbis(i) = basebis*((ce50^1.47))/(ce50^1.47+(ce(i)-ce50shift)^1.47);
    else
        predbis(i) = basebis*((ce50^1.89))/(ce50^1.89+(ce(i)-ce50shift)^1.89);
    end
end

for i = 1:1:length(ce)
    if ce(i) > ce50
        noshiftpredbis(i) = basebis*((ce50^1.47))/(ce50^1.47+(ce(i))^1.47);
    else
        noshiftpredbis(i) = basebis*((ce50^1.89))/(ce50^1.89+(ce(i))^1.89);
    end
end


plot(ce, predbis, 'g-', ce, noshiftpredbis, 'r-')