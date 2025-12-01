clear all;

nnc =38; nnc_m =2;  nnn =25; 
% fig = gcf; fig.Position =[50, 10, 860, 550]; %50, 70, 950, 250
TNR ='Times New Roman'; %font name and size
psize = [5, 5, 20, 17];

%sizes of the picture frames
w =0.405;  pw1 =0.070;  pw2 =0.568;
h =0.368;  ph2 =0.123;  ph1 =0.612; phh1 =ph1 +h; phh2 =ph2 +h;
position11 =[pw1, ph1, w, h];  position12 =[pw2, ph1, w, h];
position21 =[pw1, ph2, w, h];  position22 =[pw2, ph2, w, h];

% data
av1 = load('./dat/xpN_ln_eta1.dat');
av2 = load('./dat/xpN_ln_eta2.dat');
av3 = load('./dat/xpN_ln_eta3.dat');
av4a = load('./m_Rlog.dat');


leta = av1(:,1);
Nav  = zeros(length(av1(:,1)),3); 

xav(:,1) = av1(:,2);
xav(:,2) = av2(:,2);
xav(:,3) = av3(:,2);

pav(:,1) = av1(:,3);
pav(:,2) = av2(:,3);
pav(:,3) = av3(:,3);

Nav(:,1) = av1(:,4);
Nav(:,2) = av2(:,4);
Nav(:,3) = av3(:,4);

av4(:,:) = av4a(1:9,:); 
mav(:,1:3) = av4(:,2:4);

coef = [0.9999, 1.0000, 1.0001];
p1 = polyfit(log(leta), log(Nav(:,2)),1);
l1 = polyval(p1, log(leta));
fprintf('Slope of N = %.4f\n',p1(1));

p2 = polyfit(log(leta), log(xav(:,2)),1);
l2 = polyval(p2, log(leta));
fprintf('Slope of x = %.4f\n',p1(2));

p3 = polyfit(log(leta), log(pav(:,2).*leta),1);
l3 = polyval(p3, log(leta));
fprintf('Slope of p = %.4f\n',p1(1));

p4 = polyfit(log(av4(:,1)), log(abs(mav(:,2))),1);
l4 = polyval(p4, log(leta));
fprintf('Slope of m = %.4f\n',p4(1));


%----------------------------------------------------------------
% subfigure (b)
f1=figure(1); set(f1, 'Units', 'centimeters', 'Position', psize);
% subplot('position', position12);
ax1 = gca;
plot(ax1, log(leta),log(xav(:,1)./leta),'-^','linewidth',4,'markersize',20,'Color',[166 64 54]./255,...
    'DisplayName',['$g=',num2str(coef(1)),'$'])
hold on
plot(ax1, log(leta),log(xav(:,2)./leta),'-s','linewidth',4,'markersize',20,'Color',[16 139 150]./255,...
    'DisplayName','$g=1.0000$')
hold on
plot(ax1, log(leta),log(xav(:,3)./leta),'-p','linewidth',4,'markersize',20,'Color',[255 0 0]./255,...
    'DisplayName',['$g=',num2str(coef(3)),'$'])
xlabel(ax1, '$\log \eta$','interpreter','latex','Fontsize',nnc)
ylabel(ax1, '$\log \langle\hat{x}^2\rangle$','interpreter','latex','Fontsize',nnc)
xlim(ax1,[5.3, 14.2]);
%ylim(ax1,[-11.3, -4]);
set(ax1, 'xtick', [6, 10, 14]);
%set(ax1, 'ytick', [-11, -8, -5]);
ytickformat(ax1, '%3.0f')
%xtickformat('%.1f')
%hl = legend('Location','northwest');
%set(hl,'Interpreter','latex')
%legend boxoff;
% ha1 = annotation('textbox',[pw2+0.01, phh1-0.2, 0.01, 0.01], ...
%    'String','(b)', 'FontName', TNR, 'FitBoxToText','on');
% ha1.LineStyle = 'none'; ha1.Color ='k'; ha1.FontSize =nnc; ha1.Interpreter ='latex';
set(ax1,'FontName','Times','FontSize',nnc)
print(f1,sprintf('fb.eps'),'-depsc')

%----------------------------------------------------------------
% subfigure (c)
% subplot('position', position21);
f2 = figure(2); set(f2, 'Units', 'centimeters', 'Position', psize);
ax2 = gca;
plot(ax2, log(leta),log(pav(:,1).*leta),'-^','linewidth',4,'markersize',20,'Color',[166 64 54]./255,...
    'DisplayName',['$g=',num2str(coef(1)),'g_c$'])
hold on
plot(ax2, log(leta),log(pav(:,2).*leta),'-s','linewidth',4,'markersize',20,'Color',[16 139 150]./255,...
    'DisplayName','$g=1.0000g_c$')
hold on
plot(ax2, log(leta),log(pav(:,3).*leta),'-p','linewidth',4,'markersize',20,'Color',[255 0 0]./255,...
    'DisplayName',['$g=',num2str(coef(3)),'g_c$'])
xlabel(ax2, '$\log \eta$','interpreter','latex','Fontsize',nnc)
ylabel(ax2, '$\log \langle \hat{p}^2\rangle$','interpreter','latex','Fontsize',nnc)
xlim(ax2,[5.3, 14.2]);
ylim(ax2,[2.6, 9.5]);
set(ax2, 'xtick', [6, 10, 14]);
set(ax2, 'ytick', [3, 6, 9]);
ytickformat(ax2, '%2.0f')
%xtickformat('%.1f')
%hl = legend('Location','northwest');
%set(hl,'Interpreter','latex')
%legend boxoff;
% ha2 = annotation('textbox',[pw1+0.01, phh2-0.2, 0.01, 0.01], ...
%    'String','(c)', 'FontName', TNR, 'FitBoxToText','on');
% ha2.LineStyle = 'none'; ha2.Color ='k'; ha2.FontSize =nnc; ha2.Interpreter ='latex';
set(ax2,'FontName','Times','FontSize',nnc)
print(f2,sprintf('fc.eps'),'-depsc')
%----------------------------------------------------------------
% subfigure (a)
% subplot('position', position11);
f3 = figure(3); set(f3, 'Units', 'centimeters', 'Position', psize);
ax3 = gca;
plot(ax3, log(leta),log(Nav(:,1)),'-^','linewidth',4,'markersize',20,'Color',[166 64 54]./255,...
    'DisplayName',['$g=',num2str(coef(1)),'$'])
hold on
plot(ax3, log(leta),log(Nav(:,2)),'-s','linewidth',4,'markersize',20,'Color',[16 139 150]./255,...
    'DisplayName','$g=1.0000$')
hold on
plot(ax3, log(leta),log(Nav(:,3)),'-p','linewidth',4,'markersize',20,'Color',[255 0 0]./255,...
    'DisplayName',['$g=',num2str(coef(3)),'$'])
xlabel(ax3, '$\log \eta$','interpreter','latex','Fontsize',nnc)
ylabel(ax3, '$\log \langle \hat{N}\rangle$','interpreter','latex','Fontsize',nnc)
xlim(ax3,[5.3, 14.2]);
ylim(ax3,[-0.1,4.3]);
set(ax3, 'xtick', [6, 10, 14]);
set(ax3, 'ytick', [0, 2, 4]);
ytickformat(ax3, '%2.0f');
%xtickformat('%.1f')
hl = legend('Location','northwest');
set(hl,'Interpreter','latex')
legend boxoff;
% ha3 = annotation('textbox',[pw1+0.01, phh1-0.2, 0.01, 0.01], ...
%     'String','(a)', 'FontName', TNR, 'FitBoxToText','on');
% ha3.LineStyle = 'none'; ha3.Color ='k'; ha3.FontSize =nnc; ha3.Interpreter ='latex';
set(ax3,'FontName','Times','FontSize',nnc);
print(f3,sprintf('fa.eps'),'-depsc');
%----------------------------------------------------------------
% subfigure (d)
% subplot('position', position22);
f4 = figure(4); set(f4, 'Units', 'centimeters', 'Position', psize);
ax4 = gca;
plot(ax4, log(av4(:,1)),log(abs(mav(:,1))),'-^','linewidth',4,'markersize',20,'Color',[166 64 54]./255,...
    'DisplayName',['$g=',num2str(coef(1)),'g_c$'])
hold on;
plot(ax4, log(av4(:,1)),log(abs(mav(:,2))),'-s','linewidth',4,'markersize',20,'Color',[16 139 150]./255,...
    'DisplayName','$g=1.0000g_c$')
hold on;
plot(ax4, log(av4(:,1)),log(abs(mav(:,3))),'-p','linewidth',4,'markersize',20,'Color',[255 0 0]./255,...
    'DisplayName',['$g=',num2str(coef(3)),'g_c$'])
xlabel(ax4, '$\log \eta$','interpreter','latex','Fontsize',nnc)
ylabel(ax4, '$\log |m|$','interpreter','latex','Fontsize',nnc)
xlim(ax4,[5.3, 14.2]);
ylim(ax4,[-13.2, -8.6]);
set(ax4, 'xtick', [6, 10, 14]);
set(ax4, 'ytick', [-13, -11, -9]);
ytickformat(ax4,'%3.0f');
% xtickformat('%.1f')
% hl = legend('Location','northwest');
% set(hl,'Interpreter','latex')
% legend boxoff;
% ha4 = annotation('textbox',[pw2+0.01, phh2-0.2, 0.01, 0.01], ...
%    'String','(d)', 'FontName', TNR, 'FitBoxToText','on');
% ha4.LineStyle = 'none'; ha4.Color ='k'; ha4.FontSize =nnc; ha4.Interpreter ='latex';
set(ax4,'FontName','Times','FontSize',nnc)
print(f4,sprintf('fd.eps'),'-depsc')

% print(sprintf('fig2s.eps'),'-depsc')
