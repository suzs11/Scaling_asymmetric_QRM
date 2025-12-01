clear; clc;
d  = 2500;
hb = 1;
Wf = 1;
h  = 0;
A  = sparse(diag(sqrt(1:d-1),1));
Ad = A';
AdA = Ad*A;
Sz = sparse([1,0;0,-1]);
Sx = sparse([0,1;1, 0]);
Sp = sparse([0,1;0, 0]);
Sm = sparse([0,0;1, 0]);
gs = sparse([0;1]);    
es = sparse([1;0]);    
I_a = speye(2);
I_f = speye(d);
An  = kron(I_a, A);
Adn = kron(I_a, Ad);
x2  =  0.5*(Adn+An)*(Adn+An);
p2  = -0.5*(Adn-An)*(Adn-An);
N2  = kron(I_a,AdA);
ln  = 6:18;
leta = exp(ln);
coef = [0.9999, 1.0, 1.0001];
xav = zeros(length(ln),length(coef));
pav = zeros(length(ln),length(coef));
Nav = zeros(length(ln),length(coef));
% p = gcp('nocreate');
% % if there is no parallel pool, create one
% if isempty(p)
%    profileName = 'local';
%    numWorkers  = 20;
%    pc = parcluster(profileName);
%    pc.NumWorkers = numWorkers;
%    poolobj = parpool(pc, numWorkers);
% end
for cl = 1:3
  for l = 1:length(ln)
    eta = leta(l);
    Hatom  = 0.5*hb*kron(Sz, I_f);
    Hfield = hb*kron(I_a, Ad*A)/eta;
    gt = 0.5*coef(cl)/sqrt(eta);
    Hint = hb*gt*(kron(Sx, A + Ad)) + 0.5*h*kron(Sx,I_f);
    H = Hatom + Hfield + Hint;
    Ne = 0.05*d;
    [eigv, eige] = eigs(sparse(H),Ne,'sa');
    x2a = x2./eta; p2a = p2.*eta;
    xav(l,cl) = eigv(:,1)'*x2a*eigv(:,1);
    pav(l,cl) = eigv(:,1)'*p2a*eigv(:,1);
    Nav(l,cl) = eigv(:,1)'*N2*eigv(:,1);
  end
  %plot(teta,pav,'linewidth',4)
  %hold on;
  %xlabel('$t\eta^{1/\nu}$','interpreter','latex')
  %ylabel('$\langle x^2\rangle\eta^{-1/3}$','interpreter','latex')
end
% p = ('nocreate');
% if isempty(p) == false
%    delete(poolobj);
% end

%disp(xav(end))
%fprintf("********************\n");
%fprintf("The p is the :%f\n",pav(end));
for cl = 1:3
  fa = fopen(['./xpN_ln_eta',num2str(cl),'.dat'],'wt');
  for l = 1:length(ln)
    fprintf(fa,'%12.8f\t',leta(l),xav(l,cl),pav(l,cl),Nav(l,cl));
    fprintf(fa,'\n');
  end
  fclose(fa);
end

% Fit a straight line
p1 = polyfit(log(leta), log(Nav(:,2))',1);
p2 = polyfit(log(leta), log(xav(:,2))',1);
p3 = polyfit(log(leta), log(pav(:,2))',1);

% Plot the original data and the fitted line
l1 = polyval(p1, log(leta));
l2 = polyval(p2, log(leta));
l3 = polyval(p3, log(leta));

fprintf('Slope of N = %.4f\n',p1(1));
fprintf('Slope of x = %.4f\n',p2(1));
fprintf('Slope of p = %.4f\n',p3(1));


f1 = figure(1);
plot(log(leta),log(Nav(:,1)),'-^','linewidth',3,'markersize',15,...
    'MarkerFaceColor',[166 64 54]./255,'Color',[166 64 54]./255,...
    'DisplayName',['$g=',num2str(coef(1)),'g_c$'])
hold on
str1= ['$\beta_N = ',num2str(p1(1)),'$'];
%annotation('arrow','X',[0.6 0.5],'Y',[0.8 0.1],'LineWidth',1.5,'LineStyle','--');
text(11,3,str1,'Interpreter','latex','FontSize',24);
plot(log(leta),log(Nav(:,2)),'s','linewidth',3,'markersize',15,...
    'MarkerFaceColor',[16 139 150]./255,'Color',[16 139 150]./255,...
    'DisplayName',['$g=',num2str(coef(2)),'g_c$'])
hold on
plot(log(leta),l1,'-','linewidth',3,'DisplayName','Fitted Line',...
    'Color',[219 155 52]./255)
hold on
plot(log(leta),log(Nav(:,3)),'-p','linewidth',3,'markersize',15,...
    'MarkerFaceColor',[255 0 0]./255,'Color',[255 0 0]./255,...
    'DisplayName',['$g=',num2str(coef(3)),'g_c$'])
xlabel('$\log \eta$','interpreter','latex','Fontsize',24)
ylabel('$\log \langle N\rangle$','interpreter','latex','Fontsize',24)
%xlim([-2.0,2.0]);
%ylim([0,4]);
%ytickformat('%.0f')
%xtickformat('%.1f')
hl = legend('Location','west');
set(hl,'Interpreter','latex')
legend boxoff;
set(gca,'FontName','Times','FontSize',24)
fpath = './fig/NvEta.eps';
% print(f1,sprintf('NvEta.eps'),'-depsc')
%print(f1,fpath,'-depsc')

f2 = figure(2);
plot(log(leta),log(xav(:,1)),'-^','linewidth',3,'markersize',15,...
    'MarkerFaceColor',[166 64 54]./255,'Color',[166 64 54]./255,...
    'DisplayName',['$g=',num2str(coef(1)),'g_c$'])
hold on
str2= ['$\beta_x = ',num2str(p2(1)),'$'];
%annotation('arrow','X',[0.6 0.5],'Y',[0.8 0.1],'LineWidth',1.5,'LineStyle','--');
text(12,4,str2,'Interpreter','latex','FontSize',36);
plot(log(leta),log(xav(:,2)),'s','linewidth',3,'markersize',15,...
    'MarkerFaceColor',[16 139 150]./255,'Color',[16 139 150]./255,...
    'DisplayName',['$g=',num2str(coef(2)),'g_c$'])
hold on
plot(log(leta),l2,'-','linewidth',3,'DisplayName','Fitted Line',...
    'Color',[219 155 52]./255)
hold on
plot(log(leta),log(xav(:,3)),'-p','linewidth',3,'markersize',15,...
    'MarkerFaceColor',[255 0 0]./255,'Color',[255 0 0]./255,...
    'DisplayName',['$g=',num2str(coef(3)),'g_c$'])
xlabel('$\log \eta$','interpreter','latex','Fontsize',24)
ylabel('$\log \langle x^2\rangle$','interpreter','latex','Fontsize',24)
%xlim([-2.0,2.0]);
%ylim([0,4]);
%ytickformat('%.0f')
%xtickformat('%.1f')
hl = legend('Location','west');
set(hl,'Interpreter','latex')
legend boxoff;
set(gca,'FontName','Times','FontSize',24)
fpath = './fig/XvEta.eps';
% print(f2,sprintf('XvEta.eps'),'-depsc')
%print(f2,fpath,'-depsc')

f3 = figure(3);
plot(log(leta),log(pav(:,1)),'-^','linewidth',3,'markersize',15,...
    'MarkerFaceColor',[166 64 54]./255,'Color',[166 64 54]./255,...
    'DisplayName',['$g=',num2str(coef(1)),'g_c$'])
hold on
str3= ['$\beta_p = ',num2str(p3(1)),'$'];
%annotation('arrow','X',[0.6 0.5],'Y',[0.8 0.1],'LineWidth',1.5,'LineStyle','--');
text(10,-3,str3,'Interpreter','latex','FontSize',24);
plot(log(leta),log(pav(:,2)),'s','linewidth',3,'markersize',15,...
    'MarkerFaceColor',[16 139 150]./255,'Color',[16 139 150]./255,...
    'DisplayName',['$g=',num2str(coef(2)),'g_c$'])
hold on
plot(log(leta),l3,'-','linewidth',3,'DisplayName','Fitted Line',...
    'Color',[219 155 52]./255)
hold on
plot(log(leta),log(pav(:,3)),'-p','linewidth',3,'markersize',15,...
    'MarkerFaceColor',[255 0 0]./255,'Color',[255 0 0]./255,...
    'DisplayName',['$g=',num2str(coef(3)),'g_c$'])
xlabel('$\log \eta$','interpreter','latex','Fontsize',24)
ylabel('$\log \langle p^2\rangle$','interpreter','latex','Fontsize',24)
%xlim([-2.0,2.0]);
%ylim([0,4]);
%ytickformat('%.0f')
%xtickformat('%.1f')
hl = legend('Location','west');
set(hl,'Interpreter','latex')
legend boxoff;
set(gca,'FontName','Times','FontSize',24)
fpath='./fig/PvEta.eps';
% print(f3,sprintf('PvEta.eps'),'-depsc')
%print(f3,fpath,'-depsc')
