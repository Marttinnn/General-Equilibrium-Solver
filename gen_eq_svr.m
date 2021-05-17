clc, clear all, close all;
addpath c:\dynare\4.5.6\matlab;
%cd

%% Cell arrays for IRFs and Welfare measure
ineterestrules=2; %Number of interest rate rules (+1 for PEG)
calibrations=1; %Number of different calibrations
irfperiods=20; %Impulse Response Function periods

IRF_ALL=cell(calibrations,ineterestrules+1);
U=table(zeros(ineterestrules+1,1), zeros(ineterestrules+1,1), 'VariableNames', {'Conditional' 'Unconditional'}); %Table for saving welfare results (of each iteration / different MP)
%U.Properties.RowNames = {'CPI' 'CPIY' 'PEG'};
%% Welfare calculations
i=1; %Policy regimes
j=1; %Calibrations
phi_pis=[1.5, 1.5]; %Different weights on inflation (Interest rate rule for MP)
phi_ys=[0, 0.5]; %Different weights on output gap (How strongly monetary authority reacts)

for i=1:(ineterestrules+1)
if i<=ineterestrules
    phi_pi=phi_pis(i); %Taking different weights from phi_pis vector
    phi_y=phi_ys(i);
    save MP_weights.mat phi_pi phi_y %Saving the weights into a file that will be loaded in dynare
    dynare SOE_CITR.mod; %noclearall; 
elseif i==(ineterestrules+1) %SOE with currency peg
    dynare SOE_PEG.mod 
end

    % Non-stationary variables
    PH=ones(1,21); %Domestic price level
    P = ones(1,21); %CPI level
    NER = ones(1,21); %Nominal exchange rate level
    
    % Stationary variables
    PIH = 1+oo_.irfs.Pi_H_eps_A; %Domestic inflation
    PI = 1+oo_.irfs.Pi_eps_A; % CPI inflation
    Delta_NER = 1+oo_.irfs.delta_eps_eps_A; %Change in NER

for ii=(2:(irfperiods+1)) %Evolution of the non-stationary variables (given inflation, changes in NER)
    PH(ii) =  PH(ii-1)*PIH(ii-1);
    P(ii) = P(ii-1)*PI(ii-1);
    NER(ii) = NER(ii-1)*Delta_NER(ii-1);
end
    
    %To start from zero
    PH=PH-1;
    P=P-1;
    NER=NER-1;

    %Extracting IRF from Dynare results
irf = array2table(zeros(irfperiods,14), 'VariableNames', {'PIH' 'Y_gap' 'PI' 'ToT' ...
    'NER' 'IR' 'PH' 'CPI' 'REER' 'MC' 'N' 'C'...
    'w' 'nx'});
irf{:,1} = transpose(oo_.irfs.Pi_H_eps_A);
irf{:,2} = transpose(oo_.irfs.y_gap_eps_A);
irf{:,3} = transpose(oo_.irfs.Pi_eps_A);
irf{:,4} = transpose(oo_.irfs.S_eps_A);
irf{:,5} = transpose(NER(2:end));
irf{:,6} = transpose(oo_.irfs.R_eps_A);
irf{:,7} = transpose(PH(2:end));
irf{:,8} = transpose(P(2:end));
irf{:,9} = transpose(oo_.irfs.q_eps_A);
irf{:,10} = transpose(oo_.irfs.mc_eps_A);
irf{:,11} = transpose(oo_.irfs.N_eps_A);
irf{:,12} = transpose(oo_.irfs.C_eps_A);
irf{:,13} = transpose(oo_.irfs.w_eps_A);
irf{:,14} = transpose(oo_.irfs.nx_eps_A);

    %Saving IRFs (j~different calibration, i~different MP rule)
IRF_ALL{j,i}=irf; 

%%
U_DR= oo_.dr.inv_order_var(1); % Declaration order -> Decision Rule (DR) order // Utility is declared 1st -> oo_.dr.inv_order_var(1) gives the DR order of the utility
U{i,1}=oo_.dr.ys(1) + 0.5*oo_.dr.ghs2(U_DR); % oo_.dr.ys(Declaration Order of Utility) + 0.5*oo_.dr.ghs2(DR Order of Utility) // The stochastic steady state (Dynare Manual 4.13.4)
U{i,2}=oo_.mean(1); % Unconditional mean 
end%%


%% Plot
CITR=IRF_ALL{1,1};
CITR_Y=IRF_ALL{1,2};
PEG=IRF_ALL{1,3};

figure
for ii=1:8
subplot(4,2,ii)
plot(CITR{:,ii}, 'black', 'LineWidth', 1, 'LineStyle', '-')
hold on
plot(CITR_Y{:,ii}, 'red', 'LineWidth', 1, 'LineStyle', '-')
plot(PEG{:,ii}, 'blue', 'LineWidth', 1, 'LineStyle', '-')
hold off
grid on

if ii==1
    title('Domestic inflation')
elseif ii==2
    title('Output gap')
elseif ii==3
    title('CPI inflation')
elseif ii==4
    title('TOT') %Terms of Trade
elseif ii==5
    title('NER') %Nominal Exchange Rate
elseif ii==6
    title('R') %Nominal interest rate // Policy Rate
elseif ii==7
    title('PH')
elseif ii==8
    title('P')
end
end
legend('CITR', 'CITRY', 'PEG')

figure
for ii=1:6
subplot(4,2,ii)
plot(CITR{:,ii+8}, 'black', 'LineWidth', 1, 'LineStyle', '-')
hold on
plot(CITR_Y{:,ii}, 'red', 'LineWidth', 1, 'LineStyle', '-')
plot(PEG{:,ii}, 'blue', 'LineWidth', 1, 'LineStyle', '-')
hold off
grid on

if ii==1
    title('REER') %Real Exchange Rate
elseif ii==2
    title('MC real')
elseif ii==3
    title('N') %Labor
elseif ii==4
    title('C') %Consumption
elseif ii==5
    title('w') %Real wage
elseif ii==6
    title('nx') %Net Export
end
end
legend('CITR', 'CITR_Y', 'PEG')

