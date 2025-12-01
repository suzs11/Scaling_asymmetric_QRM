clear all; clc;
% load data
avh1 = load('./dat/xpNBF1.dat');
avh2 = load('./dat/xpNBF2.dat');
avh3 = load('./dat/xpNBF3.dat');
avh4 = load('./dat/xpNBF4.dat');
avh5 = load('./dat/xpNBF5.dat');
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

etl = 2.^[16:20]; Rc = 1.0;

f1 = figure(1);
plot(heta*etl(1),Navh(:,1)*etl(1)^(-1.05/3),'^','linewidth',3,'markersize',15, ...
    'Color',[166 64 54]./255,'DisplayName','$\eta=2^{16}$'); hold on;
plot(heta*etl(2),Navh(:,2)*etl(2)^(-1.05/3),'s','linewidth',3,'markersize',15, ...
    'Color',[255  0  0]./255,'DisplayName','$\eta=2^{17}$'); hold on;
plot(heta*etl(3),Navh(:,3)*etl(3)^(-1.05/3),'o','linewidth',3,'markersize',15, ...
    'Color',[16 139 150]./255, 'DisplayName','$\eta=2^{18}$'); hold on;
plot(heta*etl(4),Navh(:,4)*etl(4)^(-1.05/3),'d','linewidth',3,'markersize',15, ...
    'Color',[16 139 150]./255, 'DisplayName','$\eta=2^{19}$'); hold on;
plot(heta*etl(5),Navh(:,5)*etl(5)^(-1.05/3),'p','linewidth',3,'markersize',15, ...
    'Color',[255 0 0]./255, 'DisplayName','$\eta=2^{20}$');
xlabel('$h\eta^{1/\nu_h}$','interpreter','latex','Fontsize',40)
ylabel('$\langle \hat{N}\rangle\eta^{-\beta_{\hat{N}}/\nu_R}$','interpreter','latex','Fontsize',40)
xlim([-2.0,2.0]);
ylim([0.11, 0.202]);
% xticks([-2 0 2]);
yticks([0.12, 0.16, 0.20]);
xtickformat('%.1f')
ytickformat('%.2f')
hl = legend('Location','north');
set(hl,'Interpreter','latex')
legend boxoff;
set(gca,'FontName','Times','FontSize',38)
print(f1,sprintf('fig4b.eps'),'-depsc')
