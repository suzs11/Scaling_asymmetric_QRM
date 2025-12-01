clear all; clc;
format long;

% ==================================================
% System Parameters
% ==================================================
d  = 1000;     % Dimension of field Hilbert space
hb = 1;        % Reduced Planck's constant (set to 1)
Wf = 1;        % Field frequency
h  = 0;        % External field strength

% ==================================================
% Operator Definitions
% ==================================================

% Field operators (harmonic oscillator)
A  = sparse(diag(sqrt(1:d-1),1));     % Annihilation operator
Ad = A';                              % Creation operator  
AdA = Ad*A;                           % Number operator

% Atomic operators (two-level system)
Sz = sparse([1,0;0,-1]);              % Pauli Z matrix
Sx = sparse([0,1;1, 0]);              % Pauli X matrix
Sp = sparse([0,1;0, 0]);              % Raising operator
Sm = sparse([0,0;1, 0]);              % Lowering operator

% Atomic states
gs = sparse([0;1]);                   % Ground state |g>
es = sparse([1;0]);                   % Excited state |e>

% Identity operators
I_a = speye(2);    % Identity for atom
I_f = speye(d);    % Identity for field

% ==================================================
% Tensor Product Operators
% ==================================================
An  = kron(I_a, A);      % A ⊗ I
Adn = kron(I_a, Ad);     % A† ⊗ I
x2  = 0.5*(An+Adn)*(An+Adn);    % Field quadrature X²
p2  = -0.5*(Adn-An)*(Adn-An);   % Field quadrature P²
N2  = kron(I_a,AdA);            % Total number operator

% ==================================================
% Parameter Setup
% ==================================================
ln   = 6:24;                % Logarithmic range for eta
% leta = exp(ln);           % Eta values (exponential scale)
leta = [50, 500, 1000, 2000, 3000, 4000];               % Alternative: power law scale

% Initialize arrays for results
Eav = zeros(1, length(ln));    % Average energy
Eco = zeros(1, length(ln));    % Energy curvature
mav = zeros(1, length(ln));    % Average magnetization (unused)
mco = zeros(1, length(ln));    % Magnetization curvature (unused)

% Time perturbation parameters for finite difference
t  = [-0.0002, -0.0001, 0.0, 0.0001, 0.0002];
dt = t(2)-t(1);           % Time step
Et = zeros(length(t), length(ln));  % Energy vs time and eta

% ==================================================
% Parallel Computing Setup
% ==================================================
% p = gcp('nocreate');
% % Create parallel pool if none exists
% if isempty(p)
%    profileName = 'local';
%    numWorkers  = 20;
%    pc = parcluster(profileName);
%    pc.NumWorkers = numWorkers;
%    poolobj = parpool(pc, numWorkers);
% end

% ==================================================
% Main Computation Loop
% ==================================================
% Calculate ground state energy for different eta and time perturbations
for i = 1:length(t)
  for l = 1:length(leta)
    eta = leta(l);
    
    % Hamiltonian components
    Hatom  = 0.5*hb*kron(Sz, I_f);                    % Atomic Hamiltonian
    Hfield = hb*Wf*kron(I_a, Ad*A)/eta;               % Field Hamiltonian
    g = 0.5*(t(i)+1)/sqrt(eta);                       % Coupling strength
    Hint = hb*g*(kron(Sx, A)+kron(Sx, Ad)) + 0.5*h*kron(Sx,I_f); % Interaction
    
    % Total Hamiltonian
    H = Hatom + Hfield + Hint;
    
    % Number of eigenvalues to compute (increases with ln(l))
    Ne = 10*ln(l);
    
    % Compute smallest algebraic eigenvalues
    [eigv, eige] = eigs(sparse(H), Ne, 'sa');
    Et(i,l) = eige(1,1);    % Ground state energy
  end
end

% ==================================================
% Cleanup Parallel Pool
% ==================================================
% p = gcp('nocreate');
% if isempty(p) == false
%    delete(poolobj);
% end

% ==================================================
% Data Analysis and File Output
% ==================================================
% fx = fopen('./fig/dat/ec_eta3.dat', 'wt');
fx = fopen('Eg_etax.dat', 'wt');
for l = 1:length(leta)
    % Calculate average energy (shifted by 0.5)
    Eav(l) = Et(3,l) + 0.5;
    
    % Calculate energy curvature using finite difference
    %Eco(l) = (Et(2,l) - 2*Et(3,l) + Et(4,l)) / (dt*dt);
    
    % Alternative higher-order finite difference (commented)
    % Eco(l) = (Et(1,l) - 16*Et(2,l) + 30*Et(3,l) - ...
    %          16*Et(4,l) + Et(5,l))/(12*(dt^2));
    
    % Write data to file
    fprintf(fx, '%20.8f\t', leta(l), Eav(l), Eco(l));
    fprintf(fx, '\n');
end
fclose(fx);

% ==================================================
% Figure 1: Log-Log Plot of Energy vs Eta
% ==================================================
% Fit power law to energy scaling
p1 = polyfit(log(leta), log(abs(Eav)), 1);
l1 = polyval(p1, log(leta));
fprintf('Slope of energy scaling = %.4f\n', p1(1));

% Create and customize plot
f1 = figure(1);
plot(log(leta), log(abs(Eav)), '^', 'linewidth', 2, 'markersize', 10, ...
     'MarkerFaceColor', [166 64 54]./255, 'Color', [166 64 54]./255)
hold on
plot(log(leta), l1, '-', 'linewidth', 2, 'DisplayName', 'Fitted Line', ...
     'Color', [219 155 52]./255)

% Add slope annotation
str1 = ['$-d/\nu = ', num2str(p1(1)), '$'];
text(8, -8, str1, 'Interpreter', 'latex', 'FontSize', 24);

% Axis labels and formatting
xlabel('$\log \eta$', 'interpreter', 'latex', 'Fontsize', 24)
ylabel('$\log|E+1/2|$', 'interpreter', 'latex', 'Fontsize', 24)
xlim([min(log(leta)), max(log(leta))]);
ylim([min(log(abs(Eav))), max(log(abs(Eav)))]);
set(gca, 'FontName', 'Times', 'FontSize', 24)

% Save figure
fpath = './fig/Egtt1.eps';
print(f1, fpath, '-depsc')

% ==================================================
% Figure 2: Curvature vs Log Eta
% ==================================================
% Fit linear relationship to curvature
% p2 = polyfit(ln, Eco', 1);
% l2 = polyval(p2, log(leta));
% fprintf('Slope of curvature = %.4f\n', p2(1));
% 
% % Create and customize plot
% f2 = figure(2);
% plot(log(leta), Eco, '^', 'linewidth', 2, 'markersize', 10, ...
%      'MarkerFaceColor', [166 64 54]./255, 'Color', [166 64 54]./255)
% hold on
% plot(log(leta), l2, '-', 'linewidth', 2, 'DisplayName', 'Fitted Line', ...
%      'Color', [219 155 52]./255)
% 
% % Add slope annotation
% str2 = ['$\alpha = ', num2str(p2(1)), '$'];
% text(7, 0.7, str2, 'Interpreter', 'latex', 'FontSize', 24);
% 
% % Axis labels and formatting
% xlabel('$\log \eta$', 'interpreter', 'latex', 'Fontsize', 24)
% ylabel('$C$', 'interpreter', 'latex', 'Fontsize', 24)
% xlim([min(log(leta)), max(log(leta))]);
% ylim([min(Eco), max(Eco)]);
% ytickformat('%.1f')
% set(gca, 'FontName', 'Times', 'FontSize', 24)
% 
% % Save figure
% fpath = './fig/C.eps';
% print(f2, fpath, '-depsc')
