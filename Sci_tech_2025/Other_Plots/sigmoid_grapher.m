% Used to Recreate figures and Data presented in "Optimal Cislunar
% Trajectories with Continuous, High-Thrust Nuclear-Thermal Propulsion" at
% Scitech 2025. Consult README before running

savefig = false;

close all
% Define the sigmoid function
sigmoid = @(x) 1 ./ (1 + exp(-x));

% Generate x values
x = linspace(-10, 10, 100); % Adjust range and resolution as needed

% Compute y values
y = sigmoid(x);

% Plot the sigmoid function
figure;
plot(x, y, 'LineWidth', 2);
box on;
title('Sigmoid Function: $F(x)=\frac{1}{1+e^{-x}}$','Interpreter','latex',FontSize=18);

if savefig
    print(gcf, 'sigmoid', '-dsvg');
end