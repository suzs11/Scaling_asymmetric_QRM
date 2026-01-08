clear;
clc;
d = 600;
hbar = 1;
%eta  = 2^(12);
etl  =2.^[16:20];
Wf   = 1;
h    = 0;
A  = diag(sqrt(1:d-1),1);
Ad = A';
AdA = Ad*A;
Sz = [1,0;0,-1];
Sx = [0,1;1, 0];
Sp = [0,1;0, 0];
Sm = [0,0;1, 0];
gs = [0;1];
es = [1;0];
I_a = eye(2);
I_f = eye(d);
An  = kron(I_a, A);
Adn = kron(I_a, Ad);
x2  =  0.5*(An+Adn)*(An+Adn);
p2  = -0.5*(Adn-An)*(Adn-An);
N2  = kron(I_a,AdA);
R = linspace(1-2^(-9), 1+2^(-9), 80);
xav = zeros(length(etl),length(R));
pav = zeros(length(etl),length(R));
Nav = zeros(length(etl),length(R));
for l = 1:length(etl)
    eta = etl(l);
    Hatom  = 0.5*hbar*kron(Sz, I_f);
    Hfield = hbar*kron(I_a, AdA)/eta;
    for tl = 1:length(R)
        g = 0.5*R(tl)/sqrt(eta);
        Hint = hbar*g*(kron(Sx, A + Ad)) + 0.5*h*kron(Sx,I_f);
        H = Hatom + Hfield + Hint;
        Ne = 0.1*d;
        [eigv, eige] = eigs(sparse(H),Ne, 'sa');
        xav(l,tl) = eigv(:,1)'*x2*eigv(:,1);
        pav(l,tl) = eigv(:,1)'*p2*eigv(:,1);
        Nav(l,tl) = eigv(:,1)'*N2*eigv(:,1);
    end
end

for l = 1 : length(etl)
    fx = fopen(['../fig/dat/xpN',num2str(l),'.dat'],'wt');
    for tl = 1:length(R)
        fprintf(fx,'%12.8f\t', R(tl), xav(l,tl),pav(l,tl),Nav(l,tl));
        fprintf(fx,'\n');
    end
end

%disp(xav(end))
%fprintf("********************\n");
%fprintf("The p is the :%f\n",pav(end));
Rc = 1.0;
f1 = figure(1);
plot((R-Rc)*etl(1)^(2/3),xav(1,:)*etl(1)^(-1/3),'^','linewidth',3,'markersize',15,'Color',[166 64 54]./255,...
    'DisplayName','$\eta=2^{16}$')
hold on
plot((R-Rc)*etl(2)^(2/3),xav(2,:)*etl(2)^(-1/3),'s','linewidth',3,'markersize',15,'Color',[16 139 150]./255,...
    'DisplayName','$\eta=2^{17}$')
hold on
plot((R-Rc)*etl(3)^(2/3),xav(3,:)*etl(3)^(-1/3),'o','linewidth',3,'markersize',15,'Color',[255 0 0]./255,...
    'DisplayName','$\eta=2^{18}$')
plot((R-Rc)*etl(4)^(2/3),xav(4,:)*etl(4)^(-1/3),'h','linewidth',3,'markersize',15,'Color',[219 155 52]./255,...
    'DisplayName','$\eta=2^{19}$')
plot((R-Rc)*etl(5)^(2/3),xav(5,:)*etl(5)^(-1/3),'p','linewidth',3,'markersize',15,'Color',[255 0 0]./255,...
    'DisplayName','$\eta=2^{20}$')
xlabel('$t\eta^{1/\nu}$','interpreter','latex','Fontsize',50)
ylabel('$\langle x^2\rangle\eta^{-1/3}$','interpreter','latex','Fontsize',50)
xlim([-2.0,2.0]);
ylim([0,4]);
ytickformat('%.1f')
xtickformat('%.1f')
hl = legend('Location','west');
set(hl,'Interpreter','latex')
legend boxoff;
set(gca,'FontName','Times','FontSize',48)
print(f1,sprintf('xav.eps'),'-depsc')

f2 = figure(2);
plot((R-Rc)*etl(1)^(2/3),pav(1,:)*etl(1)^(1/3),'^','linewidth',3,'markersize',15,'Color',[166 64 54]./255,...
    'DisplayName','$\eta=2^{16}$')
hold on
plot((R-Rc)*etl(2)^(2/3),pav(2,:)*etl(2)^(1/3),'s','linewidth',3,'markersize',15,'Color',[16 139 150]./255,...
    'DisplayName','$\eta=2^{17}$')
hold on
plot((R-Rc)*etl(3)^(2/3),pav(3,:)*etl(3)^(1/3),'o','linewidth',3,'markersize',15,'Color',[255 0 0]./255,...
    'DisplayName','$\eta=2^{18}$')
plot((R-Rc)*etl(4)^(2/3),pav(4,:)*etl(4)^(1/3),'h','linewidth',3,'markersize',15,'Color',[219 155 52]./255, ...
    'DisplayName','$\eta=2^{19}$')
plot((R-Rc)*etl(5)^(2/3),pav(5,:)*etl(5)^(1/3),'p','linewidth',3,'markersize',15,'Color',[255 0 0]./255, ...
    'DisplayName','$\eta=2^{20}$')
xlabel('$t\eta^{1/\nu}$','interpreter','latex','Fontsize',50)
ylabel('$\langle p^2\rangle\eta^{1/3}$','interpreter','latex','Fontsize',50)
xlim([-2.0,2.0]);
ylim([0.3,1.3]);
ytickformat('%.1f')
xtickformat('%.1f')
% hl = legend('Location','southwest');
% set(hl,'Interpreter','latex')
% legend boxoff;
set(gca,'FontName','Times','FontSize',38)
print(f2,sprintf('pav.eps'),'-depsc')

f3 = figure(3);
plot((R-Rc)*etl(1)^(2/3),Nav(1,:)*etl(1)^(-1/3),'^','linewidth',3,'markersize',15,'Color',[166 64 54]./255,...
    'DisplayName','$\eta=2^{16}$')
hold on
plot((R-Rc)*etl(2)^(2/3),Nav(2,:)*etl(2)^(-1/3),'s','linewidth',3,'markersize',15,'Color',[16 139 150]./255,...
    'DisplayName','$\eta=2^{17}$')
hold on
plot((R-Rc)*etl(3)^(2/3),Nav(3,:)*etl(3)^(-1/3),'o','linewidth',3,'markersize',15,'Color',[255 0 0]./255,...
    'DisplayName','$\eta=2^{18}$')
plot((R-Rc)*etl(4)^(2/3),Nav(4,:)*etl(4)^(-1/3),'h','linewidth',3,'markersize',15,'Color',[219 155 52]./255, ...
    'DisplayName','$\eta=2^{19}$')
plot((R-Rc)*etl(5)^(2/3),Nav(5,:)*etl(5)^(-1/3),'p','linewidth',3,'markersize',15,'Color',[255 0 0]./255, ...
    'DisplayName','$\eta=2^{20}$')
xlabel('$t\eta^{1/\nu}$','interpreter','latex','Fontsize',50)
ylabel('$\langle N\rangle\eta^{-\beta_{N}/\nu}$','interpreter','latex','Fontsize',50)
xlim([-2.0,2.0]);
%ylim([0.3,1.3]);
ytickformat('%.1f')
xtickformat('%.1f')
% hl = legend('Location','west');
% set(hl,'Interpreter','latex')
% legend boxoff;
set(gca,'FontName','Times','FontSize',38)
print(f3,sprintf('Nav.eps'),'-depsc')
