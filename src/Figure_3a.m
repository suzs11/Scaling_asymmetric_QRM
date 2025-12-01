clear;
clc;
d = 600;
hbar = 1;
% eta  = 2^(12);
et  = [8, 12, 16, 20];
etl = 2.^et;
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
teta = -2:0.1:2;
xav = zeros(length(etl),length(teta));
pav = zeros(length(etl),length(teta));
Nav = zeros(length(etl),length(teta));
p = gcp('nocreate'); 
% if there is no parallel pool, create one
if isempty(p)
   profileName = 'local';
   numWorkers = 10;
   pc = parcluster(profileName);
   pc.NumWorkers = numWorkers;
   poolobj = parpool(pc, numWorkers);
end
for l = 1:length(etl)
    eta = etl(l);
    W0  = eta;
    Hatom  = 0.5*hbar*eta*kron(Sz, I_f);
    Hfield = hbar*Wf*kron(I_a, Ad*A);
    parfor tl = 1:length(teta)
        g = 0.5*(teta(tl)/(eta^(2/3))+1)*sqrt(eta);
        Hint = hbar*g*(kron(Sx, A)+kron(Sx, Ad))+0.5*h*eta*kron(Sx,I_f);
        H = Hatom + Hfield + Hint;
        [eigv, eige] = eig(H);
        xav(l,tl) = eigv(:,1)'*x2*eigv(:,1);
        pav(l,tl) = eigv(:,1)'*p2*eigv(:,1);
        Nav(l,tl) = eigv(:,1)'*N2*eigv(:,1);
        % xav(l,tl) = xav(l,tl)*eta^(-1/3);
        % pav(l,tl) = pav(l,tl)*eta^(1/3);
        % Nav(l,tl) = Nav(l,tl)*eta^(-1/3);
    end
    %plot(teta,pav,'linewidth',4)
    %hold on;
    %xlabel('$t\eta^{1/\nu}$','interpreter','latex')
    %ylabel('$\langle x^2\rangle\eta^{-1/3}$','interpreter','latex')
end
p = ('nocreate');
if isempty(p)==false
   delete(poolobj);
end

for l = 1 : length(etl)
    fx = fopen(['./fig/dat/xpNa',num2str(l),'.dat'],'wt');
    for tl = 1:length(teta)
        fprintf(fx,'%12.8f\t', teta(tl), xav(l,tl),pav(l,tl),Nav(l,tl));
        fprintf(fx,'\n');
    end
end
%disp(xav(end))
%fprintf("********************\n");
%fprintf("The p is the :%f\n",pav(end));
% figsize = [8, 5]; % Double column with small size
% figsize = [16, 10]; % Double column with small size
f1 = figure(1);
plot(teta,xav(1,:)*eta^(-1/3),'^','linewidth',3,'markersize',10,'Color',[166 64 54]./255,...
    'DisplayName',['$\eta=2^{',num2str(et(1)),'}$'])
hold on                                            
plot(teta,xav(2,:)*eta^(-1/3),'s','linewidth',3,'markersize',10,'Color',[16 139 150]./255,...
    'DisplayName',['$\eta=2^{',num2str(et(2)),'}$'])
hold on                                            
plot(teta,xav(3,:)*eta^(-1/3),'p','linewidth',3,'markersize',10,'Color',[255 0 0]./255,...
    'DisplayName',['$\eta=2^{',num2str(et(3)),'}$'])
xlabel('$t\eta^{1/\nu}$','interpreter','latex','Fontsize',26)
ylabel('$\langle x^2\rangle\eta^{-1/3}$','interpreter','latex','Fontsize',26)
xlim([-2.0,2.0]);
ylim([0,4]);
ytickformat('%.1f')
xtickformat('%.1f')
hl = legend('Location','west');
set(hl,'Interpreter','latex')
legend boxoff;
set(gca,'FontName','Times','FontSize',26)
% set(f1,'Unit','Centimeters','position',[20,10,figsize])
fpath='./fig/xav.eps';
% print(f1,sprintf('xav.eps'),'-depsc')
print(f1,fpath,'-depsc')

f2 = figure(2);
plot(teta,pav(1,:)*eta^(1/3),'^','linewidth',3,'markersize',10,'Color',[166 64 54]./255,...
    'DisplayName',['$\eta=2^{',num2str(et(1)),'}$'])
hold on
plot(teta,pav(2,:)*eta^(1/3),'s','linewidth',3,'markersize',10,'Color',[16 139 150]./255,...
    'DisplayName',['$\eta=2^{',num2str(et(2)),'}$'])
hold on
plot(teta,pav(3,:)*eta^(1/3),'p','linewidth',3,'markersize',10,'Color',[255 0 0]./255,...
    'DisplayName',['$\eta=2^{',num2str(et(3)),'}$'])
xlabel('$t\eta^{1/\nu}$','interpreter','latex','Fontsize',26)
ylabel('$\langle p^2\rangle\eta^{1/3}$','interpreter','latex','Fontsize',26)
xlim([-2.0,2.0]);
ylim([0.3,1.3]);
ytickformat('%.1f')
xtickformat('%.1f')
hl = legend('Location','southwest');
set(hl,'Interpreter','latex')
legend boxoff;
set(gca,'FontName','Times','FontSize',26)
% set(f2,'Unit','Centimeters','position',[20,10,figsize])
fpath = './fig/pav.eps';
% print(f2,sprintf('pav.eps'),'-depsc')
print(f2,fpath,'-depsc')

f3 = figure(3);
plot(teta,Nav(1,:)*eta^(-1/3),'^','linewidth',3,'markersize',10,'Color',[166 64 54]./255,...
    'DisplayName',['$\eta=2^{',num2str(et(1)),'}$'])
hold on
plot(teta,Nav(2,:)*eta^(-1/3),'s','linewidth',3,'markersize',10,'Color',[16 139 150]./255,...
    'DisplayName',['$\eta=2^{',num2str(et(2)),'}$'])
hold on
plot(teta,Nav(3,:)*eta^(-1/3),'p','linewidth',3,'markersize',10,'Color',[255 0 0]./255,...
    'DisplayName',['$\eta=2^{',num2str(et(3)),'}$'])
xlabel('$t\eta^{1/\nu}$','interpreter','latex','Fontsize',26)
ylabel('$\langle N\rangle\eta^{-\beta_{N}/\nu}$','interpreter','latex','Fontsize',26)
xlim([-2.0,2.0]);
%ylim([0.3,1.3]);
ytickformat('%.1f')
xtickformat('%.1f')
hl = legend('Location','west');
set(hl,'Interpreter','latex')
legend boxoff;
set(gca,'FontName','Times','FontSize',26)
% set(f3,'Unit','Centimeters','position',[20,10,figsize])
fpath = './fig/Nav.eps';
% print(f3,sprintf('Nav.eps'),'-depsc')
print(f3,fpath,'-depsc')
