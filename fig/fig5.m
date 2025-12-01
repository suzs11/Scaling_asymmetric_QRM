clear all; clc;

nnc =28; nnc_m =2;  nnn =25; 
fig = gcf; fig.Position =[50, 10, 530, 750]; %50, 70, 950, 250
TNR ='Times New Roman'; %font name and size
psize = [5, 5, 20, 17]; lw = 1.5; ms = 13;

%sizes of the picture frames
w =0.37;  pw1 =0.105;  pw2 =0.60;
h =0.37;  ph2 =0.11;   ph1 =0.60; phh1 =ph1 +h; phh2 =ph2 +h;
position11 =[pw1, ph1, w, h];  position12 =[pw2, ph1, w, h];
position21 =[pw1, ph2, w, h];  position22 =[pw2, ph2, w, h];
% Load data
av0 = load('./dat/eta_gapF.dat');
av1 = load('./dat/ec_eta2old.dat');
av2 = load('./dat/mmeta.dat');

% Extract theta
teta = av1(1:3:end, 1);
heta = av2(1:2:end, 1);

% Initialize arrays
Eav = zeros(length(av1(1:3:end, 2)), 1);
Eco = zeros(length(av1(1:3:end, 3)), 1);
mav = zeros(length(av2(1:2:end, 2)), 1);
mco = zeros(length(av2(1:2:end, 2)), 1);

% Assign values
Eav = av1(1:3:end, 2);
Eco = av1(1:3:end, 3);
mav = av2(1:2:end, 2);
mco = av2(1:2:end, 3);

% Fit lines
p1 = polyfit(log(av0(:,1)), log(abs(av0(:,2))), 1);
p2 = polyfit(log(teta), Eco, 1);
p3 = polyfit(log(heta), log(abs(mav)), 1);
p4 = polyfit(log(heta), log(abs(mco)), 1);

l1 = polyval(p1, log(av0(:,1)));
l2 = polyval(p2, log(teta));
l3 = polyval(p3, log(heta));
l4 = polyval(p4, log(heta));

fprintf('Slope of log|Delta E| = %.4f\n', p1(1));
fprintf('Slope of C = %.4f\n', p2(1));
fprintf('Slope of log(|m|) = %.4f\n', p3(1));
fprintf('Slope of log(|X|) = %.4f\n', p4(1));

x = 4:0.1:17; 
% Plotting
subplot('Position', position11); ax1 = gca; 
% plot(ax1, log(teta), log(abs(Eav)), '--o','linewidth',lw,'markersize',ms, 'Color', [166  64  45]./255, 'MarkerFaceColor', [166 64 45]./255); 
plot(ax1, log(av0(1:2:end,1)), log(abs(av0(1:2:end,2))), '--o','linewidth',lw,'markersize',ms, 'Color', [166  64  45]./255, 'MarkerFaceColor', [166 64 45]./255); 
% plot(log(teta), l1, 'r');
% hold on;
% plot(x, -4/3*x, '--','linewidth',lw, 'Color', [166  64  45]./255, 'MarkerFaceColor', [166 64 45]./255)
xlabel(ax1,'$\log\eta$', 'Interpreter', 'latex','Fontsize', nnc);
ylabel(ax1, '$\log |\Delta E|$', 'Interpreter', 'latex','Fontsize',nnc);
% text(11, -10, sprintf('$d=%.4f$', p1(1)), 'Interpreter', 'latex','Fontsize',36);
text(8.0, -8, '$-\tilde{d}/\nu_R=-4/3$', 'Interpreter', 'latex','Fontsize',nnc);
xlim(ax1, [min(log(av0(:,1)))-0.6, max(log(av0(:,1)))-0.2]);
ylim(ax1, [min(log(abs(av0(:,2))))+0.2, max(log(abs(av0(:,2))))+1.0]);
xticks([5 10 15]);
%yticks([-15, -10, -5]);
ytickformat('%3d')
set(ax1,'FontName','Times','FontSize',nnc)
 
subplot('Position', position12); ax2 = gca; 
plot(ax2, log(teta), Eco, 'p','linewidth',lw,'markersize',ms, 'Color', [16  139 150]./255, 'MarkerFaceColor', [16 13 150]./255); 
hold on;
%plot(x, 0*x-0.27, '--','linewidth',lw, 'Color', [16  139 150]./255, 'MarkerFaceColor', [16  139 150]./255)
plot(ax2, log(teta), l2, '--','linewidth',lw,'Color',[16 139 150]./255);
xlabel(ax2, '$\log\eta$', 'Interpreter', 'latex','Fontsize', nnc);
ylabel(ax2, '$C_v$', 'Interpreter', 'latex', 'Fontsize',nnc);
% text(8, 0, ['$\alpha=',num2str(p2(1),'%.4f'),'$'], 'Interpreter', 'latex','Fontsize',48);
text(8, 0, '$\alpha=0$', 'Interpreter', 'latex','Fontsize',nnc);
xlim(ax2, [min(log(teta))-0.6, max(log(teta))+0.6]);
ylim(ax2, [-1.1, 0.6]);
xticks([5 10 15]);
yticks([-1.0 -0.5 0 0.5]);
ytickformat('%.1f')
set(ax2, 'FontName','Times','FontSize',nnc)


subplot('Position', position21); ax3 = gca; 
plot(ax3, log(heta), log(abs(mav)), '--s','linewidth',lw,'markersize',ms, 'Color', [219  155  52]./255, 'MarkerFaceColor', [219 155 52]./255);
% plot(log(heta), l3, 'r');
% hold on;
% plot(x, -1*x/3-6.1, '--','linewidth',lw, 'Color', [219  155  52]./255, 'MarkerFaceColor', [219  155  52]./255)
% text(7, -8, ['$-\beta/\nu=-',num2str(p3(1),'%.4f'),'$'], 'Interpreter', 'latex','Fontsize',48);
text(7, -8, '$-\beta/\nu_R=-1/3$', 'Interpreter', 'latex','Fontsize',nnc);
xlabel(ax3, '$\log\eta$', 'Interpreter', 'latex','Fontsize', nnc);
ylabel(ax3, '$\log|m|$',  'Interpreter', 'latex','Fontsize', nnc);
xlim(ax3, [min(log(heta))-0.6, max(log(heta))+0.6]);
ylim(ax3, [min(log(abs(mav)))-0.3, max(log(abs(mav)))+0.3]);
xticks([5 10 15]);
ytickformat('%3d');
set(ax3, 'FontName','Times','FontSize',nnc)
 
subplot('Position', position22); ax4 = gca; 
plot(ax4, log(heta), log(abs(mco)), '--d','linewidth',lw,'markersize',ms, 'Color', [255 0 0]./255, 'MarkerFaceColor', [255  0  0]./255); 
%text(9, 3.6, ['$\gamma/\nu=',num2str(p4(1),'%.4f'),'$'], 'Interpreter', 'latex','Fontsize',48);
text(9, 3.6, '$\gamma/\nu_R=2/3$', 'Interpreter', 'latex','Fontsize',nnc);
% plot(log(heta), l4, 'k');
% hold on
%plot(x, 2*x/3-0.9, '--','linewidth',lw, 'Color', [255 0 0]./255, 'MarkerFaceColor', [255 0 0]./255)
xlabel(ax4, '$\log\eta$', 'Interpreter', 'latex','Fontsize', nnc);
ylabel(ax4, '$\log|\chi_g|$', 'Interpreter', 'latex','Fontsize',nnc);
xlim(ax4, [min(log(heta))-0.6, max(log(heta))+0.6]);
ylim(ax4, [min(log(abs(mco)))-0.6, max(log(abs(mco)))+0.6]);
xticks([5 10 15]);
ytickformat(ax4, '%5d')
set(ax4,'FontName','Times','FontSize',nnc)
print('fig5b_.eps', '-depsc');
