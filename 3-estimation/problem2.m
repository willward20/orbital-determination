%% Problem 2: The Batch Processor Algorithm
clear all 
clc

y = [1; 2; 1];
W = diag([2 1 1]);
H = [1; 1; 1];
xBar = 2;
WBar = 2;

xHats = [xBar];
xVars = [1/WBar];

for i = 1:6
    xHat = (H'*W*H + WBar)\(H'*W*y + WBar*xBar);
    xBar = xHat;
    WBar = 1/(H'*W*H + WBar);
    xHats = [xHats; xBar];
    xVars = [xVars; 1/WBar];
end
eHat = y - H*xHat

% Plot xHat
figure(1)
plot(xHats)
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
title('Estimate of $\hat{x}$ Over Time', ...
    fontweight='bold', Interpreter='latex', FontSize=16)
xlabel('Iterations', FontSize=14)
ylabel('$\hat{x}$', fontsize=18, fontweight='bold', Interpreter='latex')
grid on