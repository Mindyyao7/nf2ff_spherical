% 
% NF2FF Conversion using Spherical Scanner
%
clear 
close all
clc

disp('************************************************')
disp('          Near-To-Far-Field Conversion')
disp('          Spherical Scanner')
disp('************************************************')

addpath('misc_functions')
addpath('plot_functions')
addpath('transformation_functions')
% %% Load Data
% disp('Load Data...')
% f = 5e9;
% 
% % setup = '../measurement/SphericalScan_HornAntenna/';
% setup = 'fake_antenna';
% 
% data_ff = readtable([setup 'Farfield/farfield (f=' num2str(f*1e-9) ') [1].txt']);
% 
% scans = dir(setup);
% scan_names = {scans.name};
% scan_names(ismember(scan_names,{'.','..','Farfield','Scan-NTheta-NPhi.txt'})) = []; %removing something
% 
% data_nf= cellfun(@(scan_names) readtable([setup,scan_names,'/NearFieldProbeResults' num2str(f*1e-9) 'GHz.txt']),scan_names,'UniformOutput',false);
% 
% % Rearrange farfield table
% data_ff = table(data_ff.Var1*pi/180,data_ff.Var2*pi/180,data_ff.Var3,data_ff.Var4,data_ff.Var6);
% data_ff.Properties.VariableNames = {'theta' 'phi' 'Eabs' 'Ethetaabs' 'Ephiabs'};

% Rearrange nearfield table
% data_nf = cellfun(@rearrangeTables,data_nf,'UniformOutput',false);
% data_nf = cellfun(@rotateSphericalNFData,data_nf,'UniformOutput',false);
% 
% % Select measurements to process
% data_nf = data_nf([2]);% maybe there are multiple 
% scan_names = scan_names([2]);
% 
% disp('Done!')


% NF2FF transformation

%% Initialization: dtheta, dphi, 
disp('NF2FF Transformation...')
% ff position sampling 
dtheta_ff=3;
dphi_ff =3;
phi_range = (0+dphi_ff:dphi_ff:2*pi-dphi_ff);
theta_range = (0+dtheta_ff:dtheta_ff:pi-dphi_ff);
r_ff = 10000; 
% nf parameters
A = 0.06; % sphere of measurment 
r0 = 0.005; % minimal sphere radius 5mm 
c0 = 0.020; % minimal cylinder radius 
f = 60*10^9; % f in Hz

%% HFSS data: upper hemisphere
csvFilePath = '20240301_OCAarray_HFSS_Circular_NF_Broadside.csv';
data = readtable(csvFilePath);
data.Properties.VariableNames{'Theta_deg_'}='theta';
data.Properties.VariableNames{'Phi_deg_'}='phi';
data.Properties.VariableNames{'dB_NearEPhi_'}='Ephi_mag';
data.Properties.VariableNames{'ang_deg_NearEPhi__deg_'}='Ephi_phase';
data.Properties.VariableNames{'dB_NearETheta_'}='Etheta_mag';
data.Properties.VariableNames{'ang_deg_NearETheta__deg_'}='Etheta_phase';
data.frequency = 60*ones(size(data.theta));

Etheta_li = power(10,(data.Etheta_mag/20));
Ephi_li = power(10,(data.Ephi_mag/20));
Etheta_complex = Etheta_li.*exp(1i*data.Etheta_phase/180*pi);
Ephi_complex = Ephi_li.*exp(1i*data.Ephi_phase/180*pi);
r = A*ones(size(data.phi)); % need to be changed if frequency change
Er = 1*ones(length(data.phi), 1);
Eabs = sqrt(Ephi_complex.^2+Etheta_complex.^2+Er.^2);
data_nf = table(r, data.theta/180*pi, data.phi/180*pi, [Er, Etheta_complex, Ephi_complex], Eabs);
data_nf.Properties.VariableNames = {'r' 'theta' 'phi' 'E' 'Eabs'};
data_nf = data_nf(704:1387, :);

% %% Plot NF
% [theta, phi] = meshgrid(unique(data_nf.theta), unique(data_nf.phi));
% %% rearrange
% for i = 1:1:37
%     for j = 1:1:19
%         Etheta(i,j) = data_nf.E(19*(i-1)+j,2);
%     end
% end
% 
% for a = 1:1:37
%     for b = 1:1:19
%         Ephi(a,b) = data_nf.E(19*(a-1)+b, 3);
%     end
% end

% %% NF at a glance
% figure;
% subplot(1,2,1);
% % Etheta = reshape(10*log10(abs(Etheta_complex)),);
% surf(theta, phi, 10*log10(abs(Etheta)));
% title(sprintf('f = %f GHz',f/1000000000));
% xlabel('\theta (rad)');
% ylabel('\phi (rad)');
% zlabel('|Etheta| (dBi)');
% view(-37.5,30);
% shading flat;
% colorbar;
% 
% subplot(1,2,2);
% % Etheta = reshape(10*log10(abs(Etheta_complex)),);
% surf(theta, phi, 10*log10(abs(Ephi)));
% title(sprintf('f = %f GHz',f/1000000000));
% xlabel('\theta (rad)');
% ylabel('\phi (rad)');
% zlabel('|Ephi| (dBi)');
% view(-37.5,30);
% shading flat;
% colorbar;

% %% fake data
% f = 60*10^9;
% R = 0.005; %radius for the antenna
% k = 2*pi*f/physconst('LightSpeed');
% n = (k*R+0.045*(k*R)^(1/3)*(60)); % 11 ish
% dtheta_nf = 18;
% dphi_nf = 36;
% theta = (0:dtheta_nf:180)*pi/180; % length will be 10
% theta = reshape(theta, [length(theta), 1]);
% phi = (0:dphi_nf:360)*pi/180; % length 10\
% phi = reshape(phi, [length(phi), 1]);
% r = 0.06*ones(length(phi), 1);
% 
% N = length(phi); % Number of random complex numbers
% Etheta = 5*rand(N, 1) + 5i * rand(N, 1);
% 
% % Ephi complementary to Etheta
% Ephi = (rand(N, 1) + 1i * rand(N, 1));
% Er = 1*ones(length(Ephi), 1);
% Eabs = sqrt(Ephi.^2+Etheta.^2+Er.^2);
% data_nf = table(r, theta, phi, [Er, Etheta, Ephi], Eabs);
% data_nf.Properties.VariableNames = {'r' 'theta' 'phi' 'E' 'Eabs'};
% % data_nf_struct = table2struct(data_nf_table);
% % data_nf = num2cell(data_nf_table);


%% Transnformation Start
% Wave Number
lambda = physconst('LightSpeed')/f;
k0 = 2*pi/lambda;
% Radius of Measurement Sphere
A = mean(data_nf.r);

theta = data_nf.theta;
phi = data_nf.phi;
Etheta = data_nf.E(:,2);
Ephi = data_nf.E(:,3);

% Step Sizes
delta_theta = diff(unique(round(theta,2)));
delta_theta = delta_theta(1);
delta_phi = diff(unique(round(phi,2)));
delta_phi = delta_phi(1);

% Theta/Phi values of Far-field points to calculate
[theta_ff,phi_ff] = meshgrid(theta_range,phi_range);
theta_ff = reshape(theta_ff,numel(theta_ff),1);
phi_ff = reshape(phi_ff,numel(phi_ff),1);

% Number of spherical wave modes
N = round(k0*r0 + 5);
m = round(k0*c0 + 5);
% Preallocate space for increased speed
Etheta_ff = zeros(size(theta_ff));
Ephi_ff = zeros(size(phi_ff));

idx = 1;
for s = 1:2
    for n = 1:N
        M = repmat((-n:n), length(theta), 1);
        % Compute Spherical Expansion Coefficients
        [~, ftheta, fphi, Y] = sphericalVectorWaveFunction(s, M, n, A, theta, phi, k0);
        ftheta_t = conj(ftheta) .* sin(theta);
        fphi_t = conj(fphi) .* sin(theta);
        q = (1 ./ Y .* (ftheta_t + Ephi' * fphi_t))' * delta_theta * delta_phi;

        % Compute Far-Field 
        M = repmat((-n:n), length(theta_ff), 1);
        [xtheta,xphi] = sphericalVectorWaveFunctionFarField(s,M,n,r_ff,theta_ff,phi_ff,k0);

        Etheta_ff = Etheta_ff + sum(repmat(q, 1, length(theta_ff))' .* xtheta, 2);
        Ephi_ff = Ephi_ff + sum(repmat(q, 1, length(theta_ff))' .* xphi, 2);

        idx = idx + 1;
    end
end


% Create Results Table
s = numel(Etheta_ff);
data_nf2ff = table(zeros(s, 1), zeros(s, 1), zeros(s, 1), zeros(s, 1), zeros(s, 1));
data_nf2ff.Properties.VariableNames = {'theta', 'phi', 'Etheta', 'Ephi', 'Eabs'};
data_nf2ff.theta = reshape(theta_ff, s, 1);
data_nf2ff.phi = reshape(phi_ff, s, 1);
data_nf2ff.Etheta = reshape(Etheta_ff, s, 1);
data_nf2ff.Ephi = reshape(Ephi_ff, s, 1);
data_nf2ff.Eabs = reshape(sqrt(Etheta_ff.^2 + Ephi_ff.^2), s, 1);

disp('Done!')

% data_nf2ff = cellfun(@(data_nf) nf2ff_spherical_manual(data_nf,f,theta_range,phi_range),data_nf,'Uniformoutput',false);

% [data_nf2ff] = nf2ff_spherical_manual(data_nf, f, theta_range, phi_range);


%% nf2ff
% applies the function nf2ff to each cell of data_nf, 
% % theta_range, phi_range, for sampling at FF
% data_nf2ff = cellfun(@(data_nf) nf2ff_spherical_manual(data_nf,f,theta_range,phi_range),data_nf,'Uniformoutput',false);
% disp('Done!')
%% graphs
% disp('Done!')
% %% Plots
% disp('Plotting...')
% close all
% 
% normalized = true;
% logarithmic = true;
% 
% % Phi=0 cut
% figure('name','Far-Field Cuts,Phi=0°','numbertitle','off',...
%         'units','normalized','outerposition',[0 0 1 1]);
% plotFFPhiCut(data_ff,0,normalized,logarithmic)
% cellfun(@(data_nf2ff) plotNFPhiCutCylindrical(data_nf2ff,0,normalized,logarithmic),data_nf2ff)
% grid on
% xlabel('Theta [°]')
% if logarithmic == true
%     ylim([-50 0])
%     ylabel('E-Field Pattern [dB]')
% elseif normalized == true
%     ylim([0 1])
%     ylabel('E-Field Pattern [-]')
% end
% ylabel('E-Field Pattern [-]')
% title('Far-Field Cut Phi=0°')
% legend(['Far-Field',scan_names])
% 
% 
% figure('name','Far-Field Error,Phi=0°','numbertitle','off',...
%         'units','normalized','outerposition',[0 0 1 1]);
% cellfun(@(data_nf2ff) plotDiffPhiCutCylindrical(data_nf2ff,data_ff,0,theta_range),data_nf2ff)
% grid on
% xlabel('Theta [°]')
% ylabel('Difference to Reference Far-Field [dB]')
% title('Difference to Reference Far-Field, Phi=0°')
% legend(scan_names)
% 
% 
% % Theta=90 cut
% figure('name','Far-Field Cuts,Theta=90°','numbertitle','off',...
%         'units','normalized','outerposition',[0 0 1 1]);
% plotFFThetaCut(data_ff,pi/2,normalized,logarithmic)
% cellfun(@(data_nf2ff) plotNFThetaCutCylindrical(data_nf2ff,pi/2,normalized,logarithmic),data_nf2ff)
% grid on
% xlabel('Phi [°]')
% if logarithmic == true
%     ylim([-50 0])
%     ylabel('E-Field Pattern [dB]')
% elseif normalized == true
%     ylim([0 1])
%     ylabel('E-Field Pattern [-]')
% end
% ylabel('E-Field Pattern [-]')
% title('Far-Field Cut Theta=90°')
% legend(['Far-Field',scan_names])
% 
% 
% figure('name','Far-Field Error,Theta=90°','numbertitle','off',...
%         'units','normalized','outerposition',[0 0 1 1]);
% cellfun(@(data_nf2ff) plotDiffThetaCutCylindrical(data_nf2ff,data_ff,pi/2),data_nf2ff)
% grid on
% xlabel('Phi [°]')
% ylabel('Difference to Reference Far-Field [dB]')
% title('Difference to Reference Far-Field, Theta=90°')
% legend(scan_names)
% 
% 
% disp('Done!')
% 
% 



