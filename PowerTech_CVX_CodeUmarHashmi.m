% Created by Md Umar Hashmi
% Date: August 2019
% Place: Funchal, Madeira

% Cite the paper for referencing:
% @article{hashmi2019energy,
%   title={Energy Storage in Madeira, Portugal: Co-optimizing for Arbitrage, Self-Sufficiency, Peak Shaving and Energy Backup},
%   author={Hashmi, Md Umar and Pereira, Lucas and Bu{\v{s}}i{\'c}, Ana},
%   journal={arXiv preprint arXiv:1904.00463},
%   year={2019}
% }

% Cite SMILE H2020 Project for using the consumer load data in the folder

clear
close all
clc
load('Madeira15minTOUrates.mat')

load('ActiveLoadDataMadeira.mat')
load('ActivePVgenDataMadeira.mat')

load('probability_mat.mat')

pts = 96; %%number of points considered

P_no = P_g + P_pv;      % Net load with solar and inelastic load

price=pricedata(1:pts,3);           % ToU with 3 levels
price_nominal=pricedata(1:pts,1);   % Fixed electricity price
e_ch=0.95;                          % Charging efficiency
e_dis =0.95;                        % Discharging efficiency
del_max = 4000;                     % Maximum charging rate
del_min = -del_max;                 % Minimum discharging rate
b_initial = 1000;                   % Initial battery capacity
B_0 = b_initial*ones(pts,1);
b_rated = 2000;                     % rated battery capacity
B_max = b_rated*ones(pts,1);
B_min = 200*ones(pts,1);            % Minimum battery capacity
h=0.25;                             % sampling time of 1 minute
x_upper= del_max*ones(pts,1);
x_lower= del_min*ones(pts,1);

M = tril(ones(pts,pts));
t_incident = 24;

% Peak Power Contracts and corresponding charge per day
peak_power_contract = [3450   4600   5750   6900   10350  13800  17250  20700];
peak_charge_perDay =  [0.1611 0.2096 0.2560 0.3040 0.4478 0.5902 0.7326 0.8751];


count=1;                            % Denotes the PPC selected
lambda1 = 0.001;                  % Penalty for maintaining back-up accordance to probability of power failure

%% Convex optimization formulation for battery performing (a) arbitrage, (b) self-sufficiency, (c) peak demand shaving and (d) energy backup

cvx_begin
variables x_ch(pts,1) x_ds(pts,1) b_cap_mat(pts,1)
variables b_p_bat(pts,1) theta(pts,1) load_total(pts,1)
minimize sum(price'*theta) - lambda1*sum(prob'*(b_cap_mat))            %trace((x_ch/e_ch-x_ds*e_dis)*price')
subject to
    zeros(pts,1)<= x_ch <= x_upper;    %%charging
    zeros(pts,1)<= x_ds <= -x_lower;    %% discharging
    b_p_bat == x_ch-x_ds;
    load_total == x_ch/e_ch-x_ds*e_dis + P_g;
    theta >= zeros(pts,1);
    theta >= load_total;
    b_cap_mat == B_0 + M*b_p_bat*h;
    b_cap_mat(t_incident) >= 1600;
    B_min <= b_cap_mat <= B_max;  %%battery capacity
    load_total <= peak_power_contract(count);
cvx_end
profit_only_arbitrage = (sum(price_nominal'*P_g) - sum(price'*load_total))*h/1000 ;
B = B_0 + M*b_p_bat*h;


%% Plots 

t = 0: 1/4: 24-1/4;

figure;
subplot(211)
plot(t,P_g); hold on; plot(t,load_total)
legend('Nominal active power without storage', 'Active power seen by the energy meter with storage')
subplot(212)
plot(t,price)
legend('3 level ToU')


self_sufficiency_solar = 1 - sum(max(0, P_g))/sum(max(0, P_no));

self_sufficiency_esd = 1 - sum(max(0, load_total))/sum(max(0, P_no));


figure;
subplot(211)
plot(t,B/2000); 
legend('State of charge of the battery')

subplot(212)
plot(t, prob)
legend('Probability of power failure')



profit_fixed_cost = peak_charge_perDay(5) - peak_charge_perDay(count);

Total_profit = profit_fixed_cost + profit_only_arbitrage
figure;plot(t,P_g); hold on; plot(t, load_total); hold on; plot(t, load_total - P_g)
legend('Nominal active power without storage', 'Active power seen by the energy meter with storage', 'Storage active power')
