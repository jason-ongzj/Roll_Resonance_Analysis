% Main script for running calculations. Note that RollResonance.m must
% be in the same folder as this script.

obj = RollResonance("CN_A.txt", "CL_Delta.txt", "CL_P.txt", "CMQ.txt");
obj.InitInputData("Dynamic_Roll.txt");
obj.InitializeVars();

eps = 0.00436;      % rad
gamma = 22.5;       % deg
static_AoA = 0.05;   % deg
cg_offset = 0.005;  % m
Toffset = 0.001;    % m
mu = 0;             % deg

fin_cant_angle = 0.35;

obj.SetRefQuantities(0.205, 0.03301, fin_cant_angle)
obj.AssumeQuantities(eps, gamma, Toffset, cg_offset, static_AoA, mu);
obj.FindNaturalFrequency();

obj.InitializeBC(100, -4.3138E-15, -0.003732135, 0, 2.144275746)
for i = 101:length(obj.Time)-1
    % Roll equation via Laplace Transform
    obj.CalcDynamicRoll(i, 100); 
    
    % Steady state equation
    obj.CalcSteadyStateRoll(i);
    
    if (i <= 2100)
       obj.CalcRoots(i);
       obj.CalcComplexPlane(i);
    end
    
    % Just to show calculations are running
    if(mod(i,100) ==0)
       fprintf("Iteration %d00, calculating...\n", i/100);
   end
end

x_limit = round(max(abs(obj.sideslip(101:2100))));
y_limit = round(max(abs(obj.AoA(101:2100))));

graph_limit = x_limit + 1;
if (y_limit > x_limit)
    graph_limit = y_limit + 1;
end

figure(1);
plot(obj.Time, obj.omega);
hold on;
plot(obj.Time, obj.P/360);
plot(obj.Time, obj.Pss/360);
plot(obj.Time, obj.PssCalc/360);
xlabel('Time (s)', 'fontsize', 12);
ylabel('Frequency (Hz)', 'fontsize', 12);
xlim([0 100]);
ylim([0 8]);
legend('Natural Frequency','Dynamic - Laplace Transform', ...
    'ASTOS Output', 'Steady State Equation');
grid on;
hold off;

figure(2);
plot(obj.POmegaRatio, obj.trimAoA/obj.staticAoA);
xlabel('P / \omega', 'fontsize', 12);
ylabel('\alpha_{TRIM} / \alpha_{STATIC}', 'fontsize', 12);
grid on
xlim([0 2]);

figure(3);
ax = gca;
% Calculated sideslip and AoA
plot(obj.sideslip(101:2100), obj.AoA(101:2100));
hold on;
ax.DataAspectRatio = [1 1 1];

% Astos Output
plot(obj.AstosSS(101:2100), obj.AstosAoA(101:2100));

% Trim Arm
plot(obj.trimArmReal(101:2100), obj.trimArmImag(101:2100));

plot([-graph_limit graph_limit], [0 0], 'k');
plot([0 0], [-graph_limit graph_limit], 'k');
grid on;
xlabel('Sideslip Angle (\circ)', 'fontsize', 12);
ylabel('Angle of Attack (\circ)', 'fontsize', 12);
legend('Tricyclic Theory', 'ASTOS Output', 'Trim Arm');
title("Static AoA: " + static_AoA + "\circ, CG Offset: " + obj.CGoffset*1000 + ...
    " mm, Thrust Offset: " + obj.Toffset*1000 + " mm", 'fontsize', 14);
x0=50;
y0=50;
width=800;
height=800;
set(gcf,'position',[x0,y0,width,height])
hold off;

figure(4);
plot(obj.POmegaRatio, obj.phase);
xlabel('P / \omega', 'fontsize', 12);
ylabel('Phase Angle (\circ)', 'fontsize', 12);
grid on
xlim([0 2]);

