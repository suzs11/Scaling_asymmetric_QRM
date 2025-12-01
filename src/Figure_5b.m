clear all; clc;
format long;
d  = 15000;
hb = 1;
Wf = 1;
A  = sparse(diag(sqrt(1:d-1),1));
Ad = A';  AdA = Ad*A;

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
ln  = 6:24;
leta = exp(ln);
h   = [-0.02, -0.015, -0.01, -0.005, 0.0, 0.005, 0.01];
% dh  = h(2)-h(1);
Eh  = zeros(length(ln),length(h));
Eav = zeros(length(ln),1);
Eco = zeros(length(ln),1);
mav = zeros(length(ln),1);
mco = zeros(length(ln),1);

p = gcp('nocreate');
% if there is no parallel pool, create one
if isempty(p)
   profileName = 'local';
   numWorkers  = 20;
   pc = parcluster(profileName);
   pc.NumWorkers = numWorkers;
   poolobj = parpool(pc, numWorkers);
end
for i = 1:length(h)
    parfor l = 1:length(ln)
        eta = leta(l);
        Hatom  = 0.5*hb*kron(Sz, I_f);
        Hfield = hb*kron(I_a, Ad*A)/eta;
        gt = 0.5/sqrt(eta);
        Hint =  hb*gt*(kron(Sx, A + Ad)) + 0.5*h(i)*kron(Sx,I_f)/eta;
        H = Hatom + Hfield + Hint;
        Ne = 0.05*d;
        [eigv, eige] = eigs(sparse(H), Ne, 'sa');
        Eh(l,i) = eige(1,1);
    end
end
p = ('nocreate');
if isempty(p) == false
   delete(poolobj);
end

%disp(xav(end))
%fprintf("********************\n");
%fprintf("The p is the :%f\n",pav(end));
fx = fopen('./fig/dat/m_eta.dat','wt');
fy = fopen('./fig/dat/Ehl.dat','wt');
for l = 1:length(ln)
    dh = (h(2) - h(1)) / leta(l);
    mav(l) = (Eh(l,5) - Eh(l,4))/dh;
    mco(l) = (Eh(l,4) - 2*Eh(l,5) + Eh(l,6))/(dh*dh);
    % mav(l) = (Eh(1)-8*Eh(2)+8*Eh(4)-Eh(5))/(12*dh);
    % mco(l) = (Eh(l,1)-16*Eh(l,2)+30*Eh(l,3)-16*Eh(l,4)+Eh(l,5))/(12*(dh^2));
    fprintf(fx,'%15.8f\t',leta(l), mav(l), mco(l));
    fprintf(fx,'\n');
    fprintf(fy,'%15.8f\t',leta(l), Eh(l,:));
    fprintf(fy,'\n');
end
fclose(fx);
fclose(fy);

p1 = polyfit(log(leta),log(abs(mav))',1);
l1 = polyval(p1, log(leta));
fprintf('Slope of N = %.4f\n',p1(1));

f1 = figure(1);
plot(log(leta),log(abs(mav)),'^','linewidth',3,'MarkerFaceColor',...
    [166, 64 54]./255,'markersize',15,'Color',[166 64 54]./255)
hold on
plot(log(leta), l1,'r-','linewidth',4)
xlabel('$\log \eta$','interpreter','latex','Fontsize',24)
ylabel('$\log |m|$','interpreter','latex','Fontsize',24)
xlim([min(log(leta)),max(log(abs(leta)))]);
ylim([min(log(abs(mav))),max(log(abs(mav)))]);
%ytickformat('%.0f')
%xtickformat('%.1f')
set(gca,'FontName','Times','FontSize',24)
fpath = './fig/logm.eps';
%print(f1,sprintf('logm1.eps'),'-depsc')
print(f1, fpath,'-depsc')

p2 = polyfit(log(leta),log(abs(mco))',1);
l2 = polyval(p2, log(leta));
fprintf('Slope of N = %.4f\n',p2(1));

f2 = figure(2);
plot(log(leta),log(abs(mco)),'^','linewidth',3,'MarkerFaceColor',...
    [166, 64 54]./255,'markersize',15,'Color',[166 64 54]./255)
hold on
plot(log(leta),l2,'-','linewidth',3,'DisplayName','Fitted Line',...
    'Color',[219 155 52]./255)
xlabel('$\log \eta$','interpreter','latex','Fontsize',24)
ylabel('$\log |\chi|$','interpreter','latex','Fontsize',24)
xlim([min(log(leta)), max(log(leta))]);
ylim([min(log(abs(mco))), max(log(abs(mco)))]);
%ytickformat('%.0f')
%xtickformat('%.1f')
set(gca,'FontName','Times','FontSize',24)
fpath = './fig/logx.eps';
%print(f2,sprintf('logX1.eps'),'-depsc')
print(f2,fpath,'-depsc')
