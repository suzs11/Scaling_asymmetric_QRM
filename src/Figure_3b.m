clear;
clc;
d = 1000;
hbar = 1;
% eta  = 2^(12);
% etl  = [2^(18), 2^(19), 2^(20)];
et  = [8, 12, 16, 20];
etl = 2.^et;
Wf  = 1;
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
heta = -2:0.1:2;
xavh = zeros(length(etl),length(heta));
pavh = zeros(length(etl),length(heta));
Navh = zeros(length(etl),length(heta));
nuh = 1.0;
p = gcp('nocreate'); 
% if there is no parallel pool, create one
if isempty(p)
   profileName = 'local';
   numWorkers = 9;
   pc = parcluster(profileName);
   pc.NumWorkers = numWorkers;
   poolobj = parpool(pc, numWorkers);
end
for l = 1:length(etl)
    eta = etl(l);
    W0   = eta;
    Hatom  = 0.5*hbar*eta*kron(Sz, I_f);
    Hfield = hbar*Wf*kron(I_a, Ad*A);
    parfor hl = 1:length(heta)
        g = 0.5*sqrt(eta);
        Hint = hbar*g*(kron(Sx, A)+kron(Sx, Ad))...
               +0.5*(heta(hl)*eta^(1-1/nuh))*kron(Sx,I_f);
        H = Hatom + Hfield + Hint;
        [eigv, eige] = eig(H);
        xavh(l,hl) = eigv(:,1)'*x2*eigv(:,1);
        pavh(l,hl) = eigv(:,1)'*p2*eigv(:,1);
        Navh(l,hl) = eigv(:,1)'*N2*eigv(:,1);
        % xavh(l,hl) = xavh(l,hl)*eta^(-1/3);
        % pavh(l,hl) = pavh(l,hl)*eta^(1/3);
        % Navh(l,hl) = Navh(l,hl)*eta^(-0.33);
    end
    %plot(teta,pav,'linewidth',4)
    %hold on;
    %xlabel('$t\eta^{1/\nu}$','interpreter','latex')
    %ylabel('$\langle x^2\rangle\eta^{-1/3}$','interpreter','latex')
end
p = ('nocreate')
if isempty(p)==false
   delete(poolobj)
end
for l=1:length(etl)
    fx = fopen(['./fig/dat/xpNha',num2str(l),'.dat'],'wt')
    for hl = 1:length(heta)
        fprintf(fx,'%12.8f\t', heta(hl), xavh(l,hl),pavh(l,hl),Navh(l,hl))
        fprintf(fx,'\n')
    end
end

%disp(xav(end))
%fprintf("********************\n");
%fprintf("The p is the :%f\n",pav(end));
% figsize = [8,5];
f1 = figure(1);
plot(heta,xavh(1,:)*eta^(-1/3),'^','linewidth',3,'markersize',10,'Color',[166 64 54]./255,...
    'DisplayName',['$\eta=2^{',num2str(et(1)),'}$'])
hold on
plot(heta,xavh(2,:)*eta^(-1/3),'s','linewidth',3,'markersize',10,'Color',[16 139 150]./255,...
    'DisplayName',['$\eta=2^{',num2str(et(2)),'}$'])
hold on
plot(heta,xavh(3,:)*eta^(-1/3),'p','linewidth',3,'markersize',10,'Color',[255 0 0]./255,...
    'DisplayName',['$\eta=2^{',num2str(et(3)),'}$'])
xlabel('$h\eta^{1/\nu_h}$','interpreter','latex','Fontsize',26)
ylabel('$\langle x^2\rangle\eta^{-1/3}$','interpreter','latex','Fontsize',26)
xlim([-2.0,2.0]);
%ylim([0,4]);
ytickformat('%.1f')
xtickformat('%.1f')
hl = legend('Location','southwest');
set(hl,'Interpreter','latex')
legend boxoff;
set(gca,'FontName','Times','FontSize',26)
% set(f1,'Unit','Centimeters','position',[20,10,figsize])
fpath = './fig/xavh.eps'
% print(f1,sprintf('xavh.eps'),'-depsc')
print(f1,fpath,'-depsc')

f2 = figure(2);
plot(heta,pavh(1,:)*eta^(1/3),'^','linewidth',3,'markersize',10,'Color',[166 64 54]./255,...
    'DisplayName',['$\eta=2^{',num2str(et(1)),'}$'])
hold on
plot(heta,pavh(2,:)*eta^(1/3),'s','linewidth',3,'markersize',10,'Color',[16 139 150]./255,...
    'DisplayName',['$\eta=2^{',num2str(et(2)),'}$'])
hold on
plot(heta,pavh(3,:)*eta^(1/3),'p','linewidth',3,'markersize',10,'Color',[255 0 0]./255,...
    'DisplayName',['$\eta=2^{',num2str(3),'}$'])
xlabel('$h\eta^{1/\nu_h}$','interpreter','latex','Fontsize',26)
ylabel('$\langle p^2\rangle\eta^{1/3}$','interpreter','latex','Fontsize',26)
xlim([-2.0,2.0]);
%ylim([0.3,1.3]);
ytickformat('%.1f')
xtickformat('%.1f')
hl = legend('Location','southwest');
set(hl,'Interpreter','latex')
legend boxoff;
set(gca,'FontName','Times','FontSize',26)
% set(f2,'Unit','Centimeters','position',[20,10,figsize])
fpath = './fig/pavh.eps'
% print(f2,sprintf('pavh.eps'),'-depsc')
print(f2,fpath,'-depsc')

f3 = figure(3);
plot(heta,Navh(1,:)*eta^(-1/3),'^','linewidth',3,'markersize',10,'Color',[166 64 54]./255,...
    'DisplayName',['$\eta=2^{',num2str(et(1)),'}$'])
hold on
plot(heta,Navh(2,:)*eta^(-1/3),'s','linewidth',3,'markersize',10,'Color',[16 139 150]./255,...
    'DisplayName',['$\eta=2^{',num2str(et(2)),'}$'])
hold on
plot(heta,Navh(3,:)*eta^(-1/3),'p','linewidth',3,'markersize',10,'Color',[255 0 0]./255,...
    'DisplayName',['$\eta=2^{',num2str(et(3)),'}$'])
xlabel('$h\eta^{1/\nu_h}$','interpreter','latex','Fontsize',26)
ylabel('$\langle N\rangle\eta^{-\beta_{N}/\nu}$','interpreter','latex','Fontsize',26)
xlim([-2.0,2.0]);
ylim([0.15,0.42]);
ytickformat('%.2f')
xtickformat('%.1f')
hl = legend('Location','southwest');
set(hl,'Interpreter','latex')
legend boxoff;
set(gca,'FontName','Times','FontSize',26)
% set(f3,'Unit','Centimeters','position',[20,10,figsize])
fpath = './fig/Navh.eps'
%print(f3,sprintf('Navh.eps'),'-depsc')
print(f3,fpath,'-depsc')
