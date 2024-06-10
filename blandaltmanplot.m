function blandaltmanplot(data1, data2)

if nargin == 0
    % Example data
    data1 = [1, 2, 3, 4, 5];
    data2 = [1.1, 1.9, 3.2, 4.1, 4.8];
end

% Calculate means and differences
means = (data1 + data2) / 2;
differences = data2 - data1;

% Calculate average difference and limits of agreement
avg_diff = mean(differences);
loa = std(differences) * 1.96;

% Create the plot
figure;
scatter(means, differences);
hold on;
plot(means, avg_diff * ones(size(means)), 'k', 'LineWidth', 2); % Average line
plot(means, (avg_diff + loa) * ones(size(means)), 'r--'); % Upper limit of agreement
plot(means, (avg_diff - loa) * ones(size(means)), 'r--'); % Lower limit of agreement
hold off;

% Label the plot
xlabel('Mean of Measurements');
ylabel('Difference between Measurements');
title('Bland-Altman Plot');
