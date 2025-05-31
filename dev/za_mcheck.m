%% * Preamble
clc
clear all

set(0,'defaultAxesTickLabelInterpreter','default');
set(0,'defaultTextInterpreter','latex');
set(0,'DefaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',13);

%% * Setup Function
% Analytical MMS Solution for a damped Duffing Oscillator of the form:
% \ddot{x} + 2*z0*w0*\dot{x} + w0^2*x + al*x^3 + F*cos(Om*t).
funduff = @(abO, p) struct('v', {[-p.z0*p.w0*abO(1)-p.F/(2*p.w0)*sin(abO(2)),
                                  -(abO(3)-p.w0)+3*p.al/(8*p.w0)*abO(1)^2-p.F/(2*p.w0)/abO(1)*cos(abO(2))], ...
                                 [-p.z0*p.w0 -p.F/(2*p.w0)*cos(abO(2)) 0;
                                  3*p.al/(4*p.w0)*abO(1)+p.F/(2*p.w0)/abO(1)^2*cos(abO(2)) p.F/(2*p.w0)/abO(1)*sin(abO(2)) -1]}).v;

% xi scales x such that x = xi*xtilde .
xi = 1e-3; 
pars = struct('z0', 0.5e-2, 'w0', 2., 'al', 0.1*xi^2, 'F', 0.1/xi);

%% * Continuation
Om0 = 0.85*pars.w0;
Om1 = 1.15*pars.w0;
ds = 0.025;

Ab0 = [1.0; 0];
Sopt = struct('stepmax', 1e3, 'Dscale', [xi;1.0;1.0], 'dynamicDscale', 1);

aboC = solve_and_continue_v2(Ab0, @(abO) funduff(abO, pars), Om0, Om1, ds, Sopt);
aboC(2,aboC(1,:)<0) = pi+aboC(2,aboC(1,:)<0);
aboC(1,aboC(1,:)<0) = -aboC(1,aboC(1,:)<0);
aboC(2,:) = wrapTo2Pi(aboC(2,:));

% Plotting
lspec = {'o-', 'MarkerFaceColor', 'w', 'LineWidth', 2};

fsz = 24;
figure(1)
poss=get(gcf, 'Position');
set(gcf, 'Color', 'white', 'Position', poss)
clf()
tiledlayout(6,2,'TileSpacing','compact', 'Padding','compact');

axs = [];

axs = [axs; nexttile([3, 1])];
plot(aboC(end,:), aboC(1,:), lspec{:})
ylabel('Response Amplitude $A$ (m)')
title('Solution $u(s)$')

axs = [axs; nexttile([2 1])];
plot(aboC(end,1:end-1), diff(aboC(1,:)), lspec{:});
ylabel('A-Secant (m)')
title('Secants $u(s_{p+1})-u(s_p)$')

axs = [axs; nexttile([2 1])];
plot(aboC(end,1:end-1), diff(aboC(2,:)), lspec{:});
ylabel('$\beta$-Secant (rad)')

axs = [axs; nexttile([3 1])];
plot(aboC(end,:), aboC(2,:), lspec{:})
xlabel('Excitation Frequency $\Omega$ (rad/s)')
ylabel('Response Phase $\beta$ (rad)')

axs = [axs; nexttile([2 1])];
plot(aboC(end,1:end-1), diff(aboC(2,:)), lspec{:});
ylabel('$\Omega$-Secant (rad)')

xlabel('Excitation Frequency $\Omega$ (rad/s)')

for i=1:length(axs)
    axes(axs(i))
    set(gca, 'FontSize', fsz)
    grid on;grid minor;grid minor;box on;
    axis tight;
end
