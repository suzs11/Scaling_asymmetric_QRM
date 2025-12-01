clear; clc;

nnc =22; nnc_m =2;  nnn =25; 
% fig = gcf; fig.Position =[50, 10, 560, 435]; %50, 70, 950, 250
fig = gcf; fig.Position =[50, 10, 860, 750]; %50, 70, 950, 250
TNR ='Times New Roman'; %font name and size
psize = [5, 5, 20, 17]; lw = 3; ms = 10;

%sizes of the picture frames
w =0.39;  pw1 =0.10;  pw2 =0.58;
h =0.39;  ph2 =0.105;  ph1 =0.60; phh1 =ph1 +h; phh2 =ph2 +h;
position11 =[pw1, ph1, w, h];  position12 =[pw2, ph1, w, h];
position21 =[pw1, ph2, w, h];  position22 =[pw2, ph2, w+0.024, h];

% Define the range of R and h values
nn = 200;
%R_vec = linspace(0,  2, nn);    % R from 0 to 2
%h_vec = linspace(-1, 1, nn);    % h from -1 to 1

R_vec =  0.0:0.005:2;
h_vec = -1.0:0.01:1.2;
nh = length(h_vec);
nR = length(R_vec);

% Create a meshgrid for R and h
[R, h] = meshgrid(R_vec, h_vec);

% Define the coefficients of the cubic equation in yg
a = R.^4;
c = (1 - R.^2);
d = -R.*h/sqrt(2);

%---------------------
 Egx = zeros(size(R));
 X_min = zeros(size(R));
 x_value = -1:0.0005:1;
 E_values = zeros(1,length(x_value));
 for i = 1:nh
  for j = 1:nR
    for k = 1:length(x_value)
      E_values(k) = Ex(R(i,j), h(i,j), x_value(k));
    end
    [Egx(i,j), ind] = min(E_values);
    X_min(i, j) = x_value(ind);
  end
 end  
%------------------------ 
 f8 = figure(8);  ax8 = gca;
 i1=25; i2 = 50; i3 = 70; i5 = 25; i6 = 50; i7 = 100; i8 = 150;
 %i1= 100; i2 = 102; i3 = 89; i4 = 111; i5 = 1; i6 = 201;
 plot(ax8, h_vec, Egx(:,i1),  'k-', 'linewidth', 5, ...
      'DisplayName',['$R=',num2str(R_vec(i1)),'$']);
 hold on;
 plot(ax8, h_vec, Egx(:,i2),  'c--', 'linewidth', 5, ...
      'DisplayName',['$R=',num2str(R_vec(i2)),'$']);
 hold on;
 plot(ax8, h_vec, Egx(:,i3),  'm-', 'linewidth', 5, ...
      'DisplayName',['$R=',num2str(R_vec(i3)),'$']);
 hold on;
 plot(ax8, h_vec, Egx(:,i6),  'g--', 'linewidth', 5, ...
      'DisplayName',['$R=',num2str(R_vec(i6)),'$']);
 xlabel(ax8, '$h$', 'interpreter', 'latex');
 ylabel(ax8, '$E_g$', 'interpreter', 'latex');
 h8 = legend('Location','best');
 set(h8,'Interpreter','latex')
 legend boxoff;
 set(ax8,'FontName','Times','FontSize',32);
 print(f8, 'EXg.eps', '-depsc');     

%-------------------------------------------------------
 i11=101; i21 = 201; i31 = 221; i41 = 301; i51= 401;
 f9 = figure(9);  ax9 = gca;
 plot(ax9, h_vec(1:end-1), X_min(1:end-1,i11), 'linewidth', 5, ...
     'DisplayName',['$R=',num2str(R_vec(i11)),'$']);
 hold on;
 plot(ax9, h_vec(1:end-1), X_min(1:end-1,i21), 'linewidth', 5, ...
     'DisplayName',['$R=',num2str(R_vec(i21)),'$']);
 hold on;
 plot(ax9, h_vec(1:end-1), X_min(1:end-1,i31), 'linewidth', 5, ...
     'DisplayName',['$R=',num2str(R_vec(i31)),'$']);
 hold on;
 plot(ax9, h_vec(1:end-1), X_min(1:end-1,i41), 'linewidth', 5, ...
     'DisplayName',['$R=',num2str(R_vec(i41)),'$']);
 hold on;
 plot(ax9, h_vec(1:end-1), X_min(1:end-1,i51), 'linewidth', 5, ...
     'DisplayName',['$R=',num2str(R_vec(i51)),'$']);
 xlabel(ax9, '$h$', 'interpreter', 'latex');
 ylabel(ax9, '$x_{g}$', 'interpreter', 'latex');
 xlim(ax9, [-1, 1]);
 hl2 = legend('Location','best');
 set(hl2,'Interpreter','latex')
 legend boxoff;
 set(ax9,'FontName','Times','FontSize',32);
 print(f9, 'DEXG.eps', '-depsc');      
%------------------------------------------------------
 % Eg = -0.5 + 0.5.*(1-R.^2).*(yg.^2) + 0.25.*(R.^4).*(yg.^4)-(R.*h.*yg)./sqrt(2);
 Eg = Egx;
 i1=25; i2 = 50; i3 = 70; i4 = 95;
 h_vec(i1);
 % Create a 2D image using imagesc
 f5 = figure(5);  ax5 = gca;
 plot(ax5, R_vec, Eg(i1,:), 'linewidth',  5, ...
     'DisplayName',['$h=',num2str(h(i1)),'$']);
 hold on;
 plot(ax5, R_vec, Eg(i2,:), 'linewidth', 5, ...
     'DisplayName',['$h=',num2str(h(i2)),'$']);
 hold on;
 plot(ax5, R_vec, Eg(i3,:), 'linewidth', 5, ...
     'DisplayName',['$h=',num2str(h(i3)),'$']);
 hold on;
 plot(ax5, R_vec, Eg(i4,:), 'linewidth', 5, ...
     'DisplayName',['$h=',num2str(h(i4)),'$']);
 xlabel(ax5, '$R$', 'interpreter', 'latex');
 ylabel(ax5, '$E_g$', 'interpreter', 'latex');
 hl = legend('Location','best');
 set(hl,'Interpreter','latex')
 legend boxoff;
 set(ax5,'FontName','Times','FontSize',32);
 print(f5, 'EgR.eps', '-depsc');
 %------------------------------------------------------- 
 f6 = figure(6);  ax6 = gca;
 plot(ax6, h_vec, Eg(:,i1), 'linewidth',  5, ...
     'DisplayName',['$R=',num2str(R_vec(i1)),'$']);
 hold on;
 plot(ax6, h_vec, Eg(:,i2), 'linewidth', 5, ...
     'DisplayName',['$R=',num2str(R_vec(i2)),'$']);
 hold on;
 plot(ax6, h_vec, Eg(:,i3), 'linewidth', 5, ...
     'DisplayName',['$R=',num2str(R_vec(i3)),'$']);
 hold on;
 plot(ax6, h_vec, Eg(:,i4), 'linewidth', 5, ...
     'DisplayName',['$R=',num2str(R_vec(i4)),'$']);
 xlabel(ax6, '$h$', 'interpreter', 'latex');
 ylabel(ax6, '$E_g$', 'interpreter', 'latex');
 hl = legend('Location','southeast');
 set(hl,'Interpreter','latex')
 legend boxoff;
 set(ax6,'FontName','Times','FontSize',32);
 print(f6, 'EgH.eps', '-depsc');  
%-------------------------------------------------------
 DER  = zeros(nh, nR-1);
 for i = 1 : nh
  for j = 1 : nR-1
    DER(i,j) = (Eg(i, j+1) - Eg(i, j))/(R(i, j+1)-R(i, j));
  end
 end
 %i1=25; i2 = 50; i3 = 70; i4 = 95;
 i1= 100; i2 = 102; i3 = 91; i4 = 111; i5 = 1; i6 = 201; i7 = 101;
 h_vec(i1);
 % Create a 2D image using imagesc
 % f1 = figure(1);  ax1 = gca;
 subplot('position', position11); ax1 = gca;
 plot(ax1, R_vec(1:end-1), DER(i7,:), '-',  'linewidth', lw, ...
     'DisplayName',['$h=',num2str(h(i7)),'$']);
 hold on;
 % plot(ax1, R(i1,1:end-1),  DER(i1,:), 'r-.', 'linewidth', 5, ...
 %      'DisplayName',['$h=',num2str(h(i1)),'$'],'HandleVisibility', 'off');
 % hold on;
 plot(ax1, R_vec(1:end-1), DER(i2,:), '-.', 'linewidth', lw, ...
     'DisplayName',['$h=',num2str(h(i2)),'$']);
 hold on;
 % plot(ax1, R_vec(1:end-1), DER(i3,:), 'b:', 'linewidth', 5, ...
 %     'DisplayName',['$h=',num2str(h(i3)),'$'],'HandleVisibility', 'off');
 % hold on;
 plot(ax1, R_vec(1:end-1), DER(i4,:), ':', 'linewidth', lw, ...
     'DisplayName',['$h=',num2str(h(i4)),'$']);
 hold on;
 % plot(ax1, R_vec(1:end-1), DER(i5,:), 'c--', 'linewidth', 5, ...
 %     'DisplayName',['$h=\pm',num2str(h(i5)),'$'],'HandleVisibility', 'off');
 % hold on;
 plot(ax1, R_vec(1:end-1), DER(i6,:), '--', 'linewidth', lw, ...
     'DisplayName',['$h=',num2str(h(i6)),'$']);
 
 xlabel(ax1, '$R$', 'interpreter', 'latex');
 ylabel(ax1, '$\partial E_g/\partial R$', 'interpreter', 'latex');
 set(ax1, 'xtick',[0, 1, 2]);
 %ytickformat(ax1, '%.1f');
 hl = legend('Location','southeast');
 set(hl,'Interpreter','latex')
 legend boxoff;
 set(ax1,'FontName','Times','FontSize',nnc);
 % print(f1, 'DERtest.eps', '-depsc');   
% %-------------------------------------------------------
% % R1 = R(:,1:end-1); h1 = h(:,1:end-1);
% % yg1 = yg(1:end,1:end-1);
% % DER = -yg1.^2.*R1 - (1-R1.^2).*yg1.*dyg + 3.*R1.^2.*(yg1.^4) ...
% %       + 4.*(R1.^3).*(yg1.^2).*dyg - h1.*dyg./sqrt(2);
% % imagesc(R_vec(1:end-1), h_vec, DER);
% %-------------------------------------------------------
% 
%-------------------------------------------------------
 DEH  = zeros(nh-1, nR);
 for j = 1 : nR
  for i = 1 : nh-1
    DEH(i,j) = -(Eg(i+1, j) - Eg(i, j))/(h(i+1, j)-h(i, j));
  end
 end 
 i11=101; i21 = 201; i31 = 221; i41 = 301; i51= 401; 
 h_vec(i1);
 % Create a 2D image using imagesc
 % f2 = figure(2);  
 subplot('position', position21); ax2 = gca;
 plot(ax2, h_vec(1:end-1), DEH(:,i51), 'linewidth', lw, ...
     'DisplayName',['$R=',num2str(R_vec(i51)),'$']);
 hold on;
 plot(ax2, h_vec(1:end-1), DEH(:,i41), ':','linewidth', lw, ...
     'DisplayName',['$R=',num2str(R_vec(i41)),'$']);
 hold on;
 plot(ax2, h_vec(1:end-1), DEH(:,i31), '-.', 'linewidth', lw, ...
     'DisplayName',['$R=',num2str(R_vec(i31)),'$']);
 hold on;
 plot(ax2, h_vec(1:end-1), DEH(:,i21), '--','linewidth', lw, ...
     'DisplayName',['$R=',num2str(R_vec(i21)),'$']);
 hold on;
 plot(ax2, h_vec(1:end-1), DEH(:,i11), 'linewidth', lw, ...
     'DisplayName',['$R=',num2str(R_vec(i11)),'$']);
 hold on;
 line([0, 0],[-1, 1],'Color','black','LineStyle','--','linewidth',3,'HandleVisibility','off');
 % uistack(hp, 'bottom');
 xlabel(ax2, '$h$', 'interpreter', 'latex');
 ylabel(ax2, '$-\partial E_g/\partial h$', 'interpreter', 'latex');
 xlim(ax2, [-1, 1]);
 set(ax2, 'xtick',[-1, 0, 1]);
 set(ax2, 'ytick',[-1, 0, 1]);
 ytickformat(ax2, '%5.0f');
 hl2 = legend('Location','best');
 set(hl2,'Interpreter','latex')
 legend boxoff;
 set(ax2,'FontName','Times','FontSize',nnc);
 % print(f2, 'DEHtest.eps', '-depsc'); 

%-------------------------------------------------------
 DER2 = zeros(nh, nR-2);
 for i = 1 : nh
   for j = 1 : nR-2
    DER2(i,j) = (Eg(i, j+2) - 2.*Eg(i, j+1) + Eg(i,j))/((R(i, j+1)-R(i, j)).^2);
   end
 end
 i1= 100; i2 = 102; i3 = 91; i4 = 111; i5 = 1; i6 = 201; i7 = 101;
 h_vec(i1);
 % Create a 2D image using imagesc
 %f3 = figure(3);  ax3 = gca;
 subplot('position', position12); ax3 = gca;
 plot(ax3, R_vec(2:end-1), DER2(i7,:), 'linewidth', lw, ...
     'DisplayName',['$h=',num2str(h(i7)),'$']);
 hold on;
 % plot(ax3, R_vec(2:end-1), DER2(i1,:), 'r-.','linewidth',  5, ...
 %     'DisplayName',['$h=',num2str(h(i1)),'$'],'HandleVisibility', 'off');
 % hold on;
 plot(ax3, R_vec(2:end-1), DER2(i2,:), '-.','linewidth', lw, ...
     'DisplayName',['$h=',num2str(h(i2)),'$']);
 hold on;
 % plot(ax3, R_vec(2:end-1), DER2(i3,:), 'b:','linewidth', 5, ...
 %     'DisplayName',['$h=',num2str(h(i3)),'$'],'HandleVisibility', 'off');
 % hold on;
 plot(ax3, R_vec(2:end-1), DER2(i4,:), ':','linewidth', lw, ...
     'DisplayName',['$h=',num2str(h(i4)),'$']);
 hold on;
 % plot(ax3, R_vec(2:end-1), DER2(i5,:), 'c--','linewidth', 5, ...
 %     'DisplayName',['$h=',num2str(h(i5)),'$'],'HandleVisibility', 'off');
 % hold on;
 plot(ax3, R_vec(2:end-1), DER2(i6,:), '--','linewidth', lw, ...
     'DisplayName',['$h=',num2str(h(i6)),'$']);
 hold on;
 line([1, 1],[-3, 2],'Color','black','LineStyle','--','linewidth',3,'HandleVisibility','off');
 xlabel(ax3, '$R$', 'interpreter', 'latex');
 ylabel(ax3, '$\partial^2 E_g/\partial R^2$', 'interpreter', 'latex');
 set(ax3, 'xtick', [0, 1, 2]);
 hl = legend('Location','best');
 set(hl,'Interpreter','latex')
 legend boxoff;
 set(ax3,'FontName','Times','FontSize',nnc);
 % print(f3, 'DER2test.eps', '-depsc'); 
%-------------------------------------------------------
%  DEH2  = zeros(nn-2, nn);
%  for i = 1 : nn - 2
%   for j = 1 : nn
%     DEH2(i,j) = -(Eg(i+2, j) - 2.*Eg(i+1, j) + Eg(i,j))/((h(i+1, j)-h(i, j)).^2);
%   end
%  end
%  i1=25; i2 = 50; i3 = 70; i4 = 95;
%  h_vec(i1);
%  % Create a 2D image using imagesc
%  f4 = figure(4);  ax4 = gca;
%  plot(ax4, h_vec(2:end-1), DEH2(:,i1), 'linewidth', 5, ...
%      'DisplayName',['$R=',num2str(R_vec(i1)),'$']);
%  hold on;
%  plot(ax4, h_vec(2:end-1), DEH2(:,i2), 'linewidth', 5, ...
%      'DisplayName',['$R=',num2str(R_vec(i2)),'$']);
%  hold on;
%  plot(ax4, h_vec(2:end-1), DEH2(:,i3), 'linewidth', 5, ...
%      'DisplayName',['$R=',num2str(R_vec(i3)),'$']);
%  hold on;
%  plot(ax4, h_vec(2:end-1), DEH2(:,i4), 'linewidth', 5, ...
%      'DisplayName',['$R=',num2str(R_vec(i4)),'$']);
%  xlabel(ax4, '$h$', 'interpreter', 'latex');
%  ylabel(ax4, '$-\partial^2 E_g/\partial h^2$', 'interpreter', 'latex');
%  hl = legend('Location','best');
%  set(hl,'Interpreter','latex')
%  legend boxoff;
%  set(ax4,'FontName','Times','FontSize',32);
%  print(f4, 'DEH2test.eps', '-depsc');
% %---------------------------------------------------------
%set(gca, 'YDir', 'normal');     % Make y-axis normal (bottom to top)
% hold on;
% plot(R_vec(1:lenR), nh1(1:lenR), 'k-', 'linewidth',3);
% hold on;
% plot(R_vec(1:lenR), nh2(1:lenR), 'k-', 'linewidth',3);
% hold on;

%----------------------------------------------------------
% plot(R_vec(lenR:end), nh3(lenR:end), 'c-', 'linewidth',4);
% hold on;
% plot(R_vec(1:lenR), nh3(1:lenR), 'c--', 'linewidth',4);
% hold on;
% % plot(R_vec(lenR-1)-0.005, 0, 'k.', 'MarkerSize', 60); 
% plot(1, 0, 'b.', 'MarkerSize', 60,'Color',[0, 0.47, 0.8]); 
% hold on
% plot(2:0.01:2.2, zeros(1,21), 'k-', 'linewidth',4);
%---------------------------------------------------------
%xlabel('$R$', 'interpreter', 'latex');
%ylabel('$h$', 'interpreter', 'latex');
%xticks([0,1,2]);
%yticks([-1,0,1]);
%ytickformat('%.0f')
%xtickformat('%.0f')

%cb = colorbar;
%colormap('hot');
%caxis([-1,1]);
%cb.YTick = [-1.0, 0.0, 1.0];
%cb.YTickLabel = {'-1.0','0.0','1.0'};
%ylabel(cb, '$m$', 'Interpreter','latex', 'Rotation', 0);
%cb_label = cb.Label;
%cb_label.Position(1) = 0.7;  
%cb_label.Position(2) = 1.135; 
%set(get(cb,'ylabel'),'string','$x_g$','fontsize',20);
%delete(cb);

%set(gca,'FontName','Times','FontSize',32,'TickLength', [0 0]);
% print(f2, 'Yg_imagesc.eps', '-depsc');
 subplot('position', position22); ax = gca;
 imagesc(R_vec, h_vec(1:end-1), DEH);      % h is y-axis, R is x-axis
 set(ax, 'YDir', 'normal');        % Make y-axis normal (bottom to top)
 xlabel(ax, '$R$', 'interpreter', 'latex');
 ylabel(ax, '$h$', 'interpreter', 'latex');
 xlim(ax, [0,2]);
 ylim(ax, [-1, 1]);
 xticks([0,1,2]);
 yticks([-1,0,1]);
 ytickformat('%.0f')
 xtickformat('%.0f')
 
 cb = colorbar;
 colormap('jet');
 caxis([-1,1]);
 cb.YTick = [-1, 0, 1];
 cb.YTickLabel = {'-1','0','1'};
 ylabel(cb, '$m$', 'Interpreter','latex', 'Rotation', 0);
 cb_label = cb.Label;
 cb_label.Position(1) = 0.6;  
 cb_label.Position(2) = 1.20; 
 %set(get(cb,'ylabel'),'string','$x_g$','fontsize',20);
 %delete(cb);
 
 set(ax,'FontName','Times','FontSize',nnc,'TickLength', [0 0]);
 %print(f1, 'Xgtest.eps', '-depsc');
%-----------------------------
 function Ex = Ex(R, h, x)
  % Ex  = -1/2 + (1 - R^2)./2 .* x.^2 + 1/4 * R^4.* x.^4 - sqrt(2)/2 * R * h * x;
  Ex = -1/2 + ((1 - R^2)/2) * x.^2 + 1/4 * R^4 * x.^4 - (sqrt(2)/2) * R * h * x;
 end

 function DEh = DEh(R, h, x)
  DEh = - sqrt(2)/2 * R * x;
 end
 
 function DEr = DEr(R, h, x)
  DEr = - R.*x.^2 + R.^3.*x.^4 - sqrt(2)/2 * h * x;
 end
 
 function DEr2 = DEr2(R, h, x)
  DEr2 = - x.^2 + 3*R.^2.*x.^4;
 end
%-----------------------------
