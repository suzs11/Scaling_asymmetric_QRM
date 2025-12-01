clear;
clc;
format long;
d = 1000;
hbar = 1;
eta  = [16:20];
etl  = 2.^eta;
Wf   = 1;
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
heta = linspace(-2^(-15), 2^(-15), 80);
xavh = zeros(length(etl),length(heta));
pavh = zeros(length(etl),length(heta));
Navh = zeros(length(etl),length(heta));
nuh = 1.0;
p = gcp('nocreate');
% if there is no parallel pool, create one
if isempty(p)
    profileName = 'local';
    numWorkers  = 40;
    pc = parcluster(profileName);
    pc.NumWorkers = numWorkers;
    poolobj = parpool(pc, numWorkers);
end
for l = 1:3
    eta = etl(l);
    Hatom  = 0.5*hbar*kron(Sz, I_f);
    Hfield = hbar*Wf*kron(I_a, Ad*A)/eta;
    parfor hl = 1:length(heta)
        g = (0.5/sqrt(eta))*(-1.0/(eta^(2/3))+1);
        Hint = hbar*g*kron(Sx, A + Ad) ...
               +0.5*(heta(hl)*eta^(-1/nuh))*kron(Sx,I_f);
        H = Hatom + Hfield + Hint;
        [eigv, eige] = eig(H);
        xavh(l,hl) = eigv(:,1)'*x2*eigv(:,1);
        pavh(l,hl) = eigv(:,1)'*p2*eigv(:,1);
        Navh(l,hl) = eigv(:,1)'*N2*eigv(:,1);
    end
end

%disp(xav(end))
%fprintf("********************\n");
%fprintf("The p is the :%f\n",pav(end));
p = ('nocreate');
if isempty(p) == false
    delete(poolobj);
end

for l = 1:3
    fx = fopen(['./fig/dat/xpNF',num2str(l),'.dat'],'wt');
    for hl = 1:length(heta)
        fprintf(fx,'%18.9f\t',heta(hl),xavh(l,hl),pavh(l,hl),Navh(l,hl));
        fprintf(fx,'\n');
    end
    fclose(fx);
end

f1 = figure(1);
plot(heta*etl(1), Navh(1,:)*etl(1)^(-1/3),'^','linewidth',3,'markersize',15,'Color',[166 64 54]./255,...
    'DisplayName','$\eta=2^{16}$')
hold on
plot(heta*etl(2), Navh(2,:)*etl(2)^(-1/3),'s','linewidth',3,'markersize',15,'Color',[16 139 150]./255,...
    'DisplayName','$\eta=2^{17}$')
hold on
plot(heta*etl(3), Navh(3,:)*etl(3)^(-1/3),'o','linewidth',3,'markersize',15,'Color',[255 0 0]./255,...
    'DisplayName','$\eta=2^{18}$')
plot(heta*etl(4), Navh(4,:)*etl(4)^(-1/3),'h','linewidth',3,'markersize',15,'Color',[219 155 52]./255,...
    'DisplayName','$\eta=2^{19}$')
plot(heta*etl(5), Navh(5,:)*etl(5)^(-1/3),'p','linewidth',3,'markersize',15,'Color',[255 0 0]./255,...
    'DisplayName','$\eta=2^{20}$')
xlabel('$h\eta^{1/\nu_h}$','interpreter','latex','Fontsize',40)
ylabel('$\langle N\rangle\eta^{-\beta_{N}/\nu}$','interpreter','latex','Fontsize',40)
xlim([-2.0,2.0]);
% ylim([0.15,0.42]);
ytickformat('%.2f')
xtickformat('%.1f')
hl = legend('Location','southwest');
set(hl,'Interpreter','latex')
legend boxoff;
set(gca,'FontName','Times','FontSize',38)
fpath = './fig/Nfig4.eps';
% print(f1,sprintf('NFig4.eps'),'-depsc')
print(f1, fpath, '-depsc')
