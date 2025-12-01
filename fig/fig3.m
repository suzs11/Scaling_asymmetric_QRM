clear all; clc;
% load data
av1 = load('./dat/xpN1.dat');
av2 = load('./dat/xpN2.dat');
av3 = load('./dat/xpN3.dat');
av4 = load('./dat/xpN4.dat');
av5 = load('./dat/xpN5.dat');
teta = av1(:,1);
xav = zeros(length(av1(:,1)),5);
pav = zeros(length(av1(:,1)),5);
Nav = zeros(length(av1(:,1)),5);
xav(:,1) = av1(:,2); xav(:,2) = av2(:,2); xav(:,3) = av3(:,2);
pav(:,1) = av1(:,3); pav(:,2) = av2(:,3); pav(:,3) = av3(:,3);
Nav(:,1) = av1(:,4); Nav(:,2) = av2(:,4); Nav(:,3) = av3(:,4);
xav(:,4) = av4(:,2); xav(:,5) = av5(:,2);
pav(:,4) = av4(:,3); pav(:,5) = av5(:,3);
Nav(:,4) = av4(:,4); Nav(:,5) = av5(:,4);

% load data
avh1 = load('./dat/xpNh1.dat');
avh2 = load('./dat/xpNh2.dat');
avh3 = load('./dat/xpNh3.dat');
avh4 = load('./dat/xpNh4.dat');
avh5 = load('./dat/xpNh5.dat');
heta = avh1(:,1);
xavh = zeros(length(avh1(:,1)),5);
pavh = zeros(length(avh1(:,1)),5);
Navh = zeros(length(avh1(:,1)),5);
xavh(:,1) = avh1(:,2); xavh(:,2) = avh2(:,2); xavh(:,3) = avh3(:,2);
pavh(:,1) = avh1(:,3); pavh(:,2) = avh2(:,3); pavh(:,3) = avh3(:,3);
Navh(:,1) = avh1(:,4); Navh(:,2) = avh2(:,4); Navh(:,3) = avh3(:,4);
xavh(:,4) = avh4(:,2); xavh(:,5) = avh5(:,2);
pavh(:,4) = avh4(:,3); pavh(:,5) = avh5(:,3);
Navh(:,4) = avh4(:,4); Navh(:,5) = avh5(:,4);

etl = 2.^[16:20];
Rc = 1.0;
f1 = figure(1);
plot((teta-Rc)*etl(1)^(2/3),xav(:,1)*etl(1)^(-1/3),'^','linewidth',3,'markersize',15,'Color',[166 64 54]./255,...
    'DisplayName','$\eta=2^{16}$')
hold on
plot((teta-Rc)*etl(2)^(2/3),xav(:,2)*etl(2)^(-1/3),'s','linewidth',3,'markersize',15,'Color',[255  0  0]./255,...
    'DisplayName','$\eta=2^{17}$')
hold on
plot((teta-Rc)*etl(3)^(2/3),xav(:,3)*etl(3)^(-1/3),'o','linewidth',3,'markersize',15,'Color',[16 139 150]./255,...
    'DisplayName','$\eta=2^{18}$')
hold on
plot((teta-Rc)*etl(4)^(2/3),xav(:,4)*etl(4)^(-1/3),'d','linewidth',3,'markersize',15,'Color',[16 139 150]./255,...
    'DisplayName','$\eta=2^{19}$')
hold on
plot((teta-Rc)*etl(5)^(2/3),xav(:,5)*etl(5)^(-1/3),'p','linewidth',3,'markersize',15,'Color',[255  0  0]./255,...
    'DisplayName','$\eta=2^{20}$')
xlabel('$t\eta^{1/\nu_{R}}$','interpreter','latex','Fontsize',55)
% ylabel('$\langle \hat{x}^2\rangle\eta^{2/3}$','interpreter','latex','Fontsize',55)
ylabel('$\langle \hat{x}^2\rangle\eta^{\beta_{x^2}/\nu_R}$','interpreter','latex','Fontsize',55)
xlim([-2.0,2.0]);
ylim([0,3.8]);
xticks([-2 0 2]);
yticks([0  2]);
ytickformat('%.1f')
xtickformat('%.1f')
hl = legend('Location','west');
set(hl,'Interpreter','latex')
legend boxoff;
set(gca,'FontName','Times','FontSize',50)
print(f1,sprintf('xav1.eps'),'-depsc')

f2 = figure(2);
plot((teta-Rc)*etl(1)^(2/3),pav(:,1)*etl(1)^(1/3),'^','linewidth',3,'markersize',15,'Color',[166 64 54]./255,...
    'DisplayName','$\eta=2^{16}$')                  
hold on                                             
plot((teta-Rc)*etl(2)^(2/3),pav(:,2)*etl(2)^(1/3),'s','linewidth',3,'markersize',15,'Color',[255  0  0]./255,...
    'DisplayName','$\eta=2^{17}$')                  
hold on                                             
plot((teta-Rc)*etl(3)^(2/3),pav(:,3)*etl(3)^(1/3),'o','linewidth',3,'markersize',15,'Color',[16 139 150]./255,...
    'DisplayName','$\eta=2^{18}$')                  
hold on                                             
plot((teta-Rc)*etl(4)^(2/3),pav(:,4)*etl(4)^(1/3),'d','linewidth',3,'markersize',15,'Color',[16 139 150]./255,...
    'DisplayName','$\eta=2^{19}$')                  
hold on                                             
plot((teta-Rc)*etl(5)^(2/3),pav(:,5)*etl(5)^(1/3),'p','linewidth',3,'markersize',15,'Color',[255  0  0]./255,...
    'DisplayName','$\eta=2^{20}$')
xlabel('$t\eta^{1/\nu_{R}}$','interpreter','latex','Fontsize',55)
% ylabel('$\langle \hat{p}^2\rangle\eta^{-2/3}$','interpreter','latex','Fontsize',55)
ylabel('$\langle \hat{p}^2\rangle\eta^{-\beta_{p^2}/\nu_R}$','interpreter','latex','Fontsize',55)
xlim([-2.0,2.0]);
ylim([0.33,1.3]);
xticks([-2 0 2]);
yticks([0.5 1.0]);
ytickformat('%.1f')
xtickformat('%.1f')
% hl = legend('Location','southwest');
% set(hl,'Interpreter','latex')
% legend boxoff;
set(gca,'FontName','Times','FontSize',50)
print(f2,sprintf('pav1.eps'),'-depsc')

f3 = figure(3); % 166 64 54
plot((teta-Rc)*etl(1)^(2/3),Nav(:,1)*etl(1)^(-1/3),'^','linewidth',3,'markersize',15,'Color',[166 64 54]./255,...
    'DisplayName','$\eta=2^{16}$')                    
hold on                                               
plot((teta-Rc)*etl(2)^(2/3),Nav(:,2)*etl(2)^(-1/3),'s','linewidth',3,'markersize',15,'Color',[255  0  0]./255,...
    'DisplayName','$\eta=2^{17}$')                    
hold on                                               
plot((teta-Rc)*etl(3)^(2/3),Nav(:,3)*etl(3)^(-1/3),'o','linewidth',3,'markersize',15,'Color',[16 139 150]./255,...
    'DisplayName','$\eta=2^{18}$')                     
hold on                                               
plot((teta-Rc)*etl(4)^(2/3),Nav(:,4)*etl(4)^(-1/3),'d','linewidth',3,'markersize',15,'Color',[16 139 150]./255,...
    'DisplayName','$\eta=2^{19}$')                    
hold on                                               
plot((teta-Rc)*etl(5)^(2/3),Nav(:,5)*etl(5)^(-1/3),'p','linewidth',3,'markersize',15,'Color',[255  0  0]./255,...
    'DisplayName','$\eta=2^{20}$')
xlabel('$t\eta^{1/\nu_{R}}$','interpreter','latex','Fontsize',55)
ylabel('$\langle \hat{N}\rangle\eta^{-\beta_{N}/\nu_R}$','interpreter','latex','Fontsize',55)
xlim([-2.0,2.0]);
ylim([0,1.9]);
xticks([-2 0 2]);
yticks([0  1]);
ytickformat('%.1f')
xtickformat('%.1f')
% hl = legend('Location','west');
% set(hl,'Interpreter','latex')
% legend boxoff;
set(gca,'FontName','Times','FontSize',50)
print(f3,sprintf('Nav1.eps'),'-depsc')

f4 = figure(4);
plot(heta*etl(1),xavh(:,1)*etl(1)^(-1/3),'^','linewidth',3,'markersize',15,'Color',[166 64 54]./255,...
    'DisplayName','$\eta=2^{16}$')
hold on
plot(heta*etl(2),xavh(:,2)*etl(2)^(-1/3),'s','linewidth',3,'markersize',15,'Color',[255  0  0]./255,...
    'DisplayName','$\eta=2^{17}$')
hold on
plot(heta*etl(3),xavh(:,3)*etl(3)^(-1/3),'o','linewidth',3,'markersize',15,'Color',[16 139 150]./255,...
    'DisplayName','$\eta=2^{18}$')
hold on
plot(heta*etl(4),xavh(:,4)*etl(4)^(-1/3),'d','linewidth',3,'markersize',15,'Color',[16 139 150]./255,...
    'DisplayName','$\eta=2^{19}$')
hold on
plot(heta*etl(5),xavh(:,5)*etl(5)^(-1/3),'p','linewidth',3,'markersize',15,'Color',[255  0  0]./255,...
    'DisplayName','$\eta=2^{20}$')
xlabel('$h\eta^{1/\nu_h}$','interpreter','latex','Fontsize',55)
% ylabel('$\langle \hat{x}^2\rangle\eta^{2/3}$','interpreter','latex','Fontsize',55)
ylabel('$\langle \hat{x}^2\rangle\eta^{\beta_{x^2}/\nu_R}$','interpreter','latex','Fontsize',55)
xlim([-2.0,2.0]);
ylim([0.41,1.1]);
xticks([-2 0 2]);
yticks([0.5 1.0]);
ytickformat('%.1f')
xtickformat('%.1f')
% hl = legend('Location','west');
% set(hl,'Interpreter','latex')
% legend boxoff;
set(gca,'FontName','Times','FontSize',50)
print(f4,sprintf('xavh1.eps'),'-depsc')

f5 = figure(5);
plot(heta*etl(1),pavh(:,1)*etl(1)^(1/3),'^','linewidth',3,'markersize',15,'Color',[166 64 54]./255,...
    'DisplayName','$\eta=2^{16}$')
hold on
plot(heta*etl(2),pavh(:,2)*etl(2)^(1/3),'s','linewidth',3,'markersize',15,'Color',[255  0  0]./255,...
    'DisplayName','$\eta=2^{17}$')
hold on
plot(heta*etl(3),pavh(:,3)*etl(3)^(1/3),'h','linewidth',3,'markersize',15,'Color',[16 139 150]./255,...
    'DisplayName','$\eta=2^{18}$')
hold on
plot(heta*etl(4),pavh(:,4)*etl(4)^(1/3),'d','linewidth',3,'markersize',15,'Color',[16 139 150]./255,...
    'DisplayName','$\eta=2^{19}$')
hold on
plot(heta*etl(5),pavh(:,5)*etl(5)^(1/3),'p','linewidth',3,'markersize',15,'Color',[255  0  0]./255,...
    'DisplayName','$\eta=2^{20}$')
xlabel('$h\eta^{1/\nu_h}$','interpreter','latex','Fontsize',55)
% ylabel('$\langle \hat{p}^2\rangle\eta^{-2/3}$','interpreter','latex','Fontsize',55)
ylabel('$\langle \hat{p}^2\rangle\eta^{-\beta_{p^2}/\nu_R}$','interpreter','latex','Fontsize',55)
xlim([-2.0,2.0]);
ylim([0.54, 0.88]);
xticks([-2 0 2]);
yticks([0.6 0.8]);
ytickformat('%.1f')
xtickformat('%.1f')
% hl = legend('Location','west');
% set(hl,'Interpreter','latex')
% legend boxoff;
set(gca,'FontName','Times','FontSize',50)
print(f5,sprintf('pavh1.eps'),'-depsc')

f6 = figure(6);
plot(heta*etl(1),Navh(:,1)*etl(1)^(-1/3),'^','linewidth',3,'markersize',15,'Color',[166 64 54]./255,...
    'DisplayName','$\eta=2^{16}$')
hold on
plot(heta*etl(2),Navh(:,2)*etl(2)^(-1/3),'s','linewidth',3,'markersize',15,'Color',[255  0  0]./255,...
    'DisplayName','$\eta=2^{17}$')
hold on
plot(heta*etl(3),Navh(:,3)*etl(3)^(-1/3),'h','linewidth',3,'markersize',15,'Color',[16 139 150]./255,...
    'DisplayName','$\eta=2^{18}$')
hold on
plot(heta*etl(4),Navh(:,4)*etl(4)^(-1/3),'d','linewidth',3,'markersize',15,'Color',[16 139 150]./255,...
    'DisplayName','$\eta=2^{19}$')
hold on
plot(heta*etl(5),Navh(:,5)*etl(5)^(-1/3),'p','linewidth',3,'markersize',15,'Color',[255  0  0]./255,...
    'DisplayName','$\eta=2^{20}$')
xlabel('$h\eta^{1/\nu_h}$','interpreter','latex','Fontsize',55)
ylabel('$\langle \hat{N}\rangle\eta^{-\beta_{N}/\nu_R}$','interpreter','latex','Fontsize',55)
xlim([-2.0,2.0]);
ylim([0.2, 0.55]);
xticks([-2 0 2]);
yticks([0.2 0.4]);
ytickformat('%.1f')
xtickformat('%.1f')
% hl = legend('Location','west');
% set(hl,'Interpreter','latex')
% legend boxoff;
set(gca,'FontName','Times','FontSize',50)
print(f6,sprintf('Navh1.eps'),'-depsc')
