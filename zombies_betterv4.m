%20.320 Project 2
%Modeling Train to Busan zombie outbreak
%2020 11 02

%RUN THIS SECTION FIRST
close all
clear all
clc

%Initial values
S_o = 51249999;     %suceptible; #people; population of South Korea in 2016
R_o = 0;            %removed (dead zombies)
D_o = 0;            %naturally dead
Z_o = 1;            %zombies
Q_o = 0;            %quarantined
N = S_o + Q_o + Z_o;%total population

%Rate constants
Vb  = 0;           %birth rate; zero on short timescale 
Ksd = 0;           %natural death rate; zero on short timescale
Ksz = 10/24;       %hours-1, rate of collision/zombification of suceptible
Kdz = 0;           %zombification rate of naturally dead; set to zero, as no D population
Kzr = Ksz*1/5000;  %hours-1; zombie destruction rate; rate of collision x probability someone can kill a zombie
Kqs = 0;           %rate of exiting quarantine; assume no leaving quarantine on short timescale
Ksq = 1/120;       %hours-1, rate of entering quarantine

%Multiplication factors for changing rate constants
Alpha = 0.01;      %accounting for increased skill at avoiding zombies over time
Beta = 0.01;       %accounting for increased skill at killing zombies over time

tspan = linspace(0,120); %hours

options = odeset('AbsTol', 1e-8, 'RelTol', 1e-8);%set accuracy 

%% No quarantine
%Basic zombie outbreak model with no quarantine of healthy individuals

y0 = [S_o, D_o, Z_o, R_o];%initial values
k = [Vb,Ksd, Ksz, Kdz, Kzr]; %rate constants

[t, p] = ode23s(@zombies, tspan, y0 , options, k);

%Linear plots
tiledlayout(2,2);
nexttile
plot(t, p(:,1), '-b') %Suceptible
hold on
plot(t, p(:,3), '-r') %Zombies
xlabel('Time (hours)')
ylabel('# individuals')
title('No quarantine')
legend(["Susceptible", "Zombies"])
hold off

nexttile
plot(t, p(:,4), '-k') %Removed
hold on
% plot(t, p(:,2), '-c') %Naturally dead
xlabel('Time (hours)')
ylabel('# individuals')
legend('Removed (dead zombies)')
title('No quarantine')
hold off

%Fraction of population plots, semilogy
%Assume little change in N over timescale
nexttile
semilogy(t, p(:,1)./N, '-b') %Suceptible
hold on
semilogy(t, p(:,3)./N, '-r') %Zombies
semilogy(t, p(:,4)./N, '-k') %Removed
% semilogy(t, p(:,2)./N, '-c') %Dead
xlabel('Time (hours)')
ylabel('Fraction of population')
title('No quarantine')
legend('Susceptible', 'Zombies','Removed (dead zombies)')
hold off

% Plot total population
% Note: variation in total due to variants in accuracy of ODE solver
% figure
% plot(t, p(:,1)+p(:,2)+p(:,3)+p(:,4))

%% Quarantine 
%Zombie outbreak model with quarantine of healthy individuals

y0 = [S_o, D_o, Z_o, R_o, Q_o]; %initial values
k = [Vb,Ksd, Ksz, Kdz, Kzr, Kqs, Ksq]; %rate constants

[t_q, p_q] = ode15s(@zombies_quarantine, tspan, y0 , options, k);

%Linear plots
tiledlayout(2,2);
nexttile
plot(t_q, p_q(:,1), 'b') %suceptible
hold on
plot(t_q, p_q(:,3), 'r') %zombies
plot(t_q, p_q(:,5), '-g')%quarantined
xlabel('Time (hours)')
ylabel('# individuals')
title('Quarantine')
legend("Susceptible", "Zombies", "Quarantined")
hold off

nexttile
plot(t_q, p_q(:,4), 'k') %removed
hold on 
% plot(t_q, p_q(:,2), 'c')%naturally dead
xlabel('Time (hours)')
ylim([0,8000]);
ylabel('# individuals')
title('Quarantine')
legend("Removed (dead zombies)")
hold off

%Fraction of population plots, semilogy
%Assume little change in N over timescale
nexttile
semilogy(t_q, p_q(:,1)./N, 'b') %suceptible
hold on
semilogy(t_q, p_q(:,3)./N, 'r') %zombies
semilogy(t_q, p_q(:,5)./N, '-g') %quarantined
semilogy(t_q, p_q(:,4)./N, 'k') %removed
% semilogy(t_q, p_q(:,2)./N, 'c') %naturally dead
xlabel('Time (hours)')
ylabel('Fraction of population')
title('Quarantine')
legend("Susceptible", "Zombies", "Quarantined", "Removed (dead zombies)")
hold off

%% No quarantine, changing rates
%Basic zombie outbreak model with no quarantine of healthy individuals and
%increased skill at zombie avoidance/killing over time

y0 = [S_o, D_o, Z_o, R_o, Ksz, Kzr]; %initial values
k = [Vb,Ksd, Kdz, Kzr, Alpha, Beta]; %rates

[t, p_a] = ode23s(@zombies_changing_rates,tspan,y0, options,k);

%Linear plots
tiledlayout(3,2)
nexttile
plot(t, p_a(:,1), '-b') %suceptible
hold on
plot(t, p_a(:,3), '-r') %zombies
% plot(t, p_a(:,2), '-c') %naturally dead
xlabel('Time (hours)')
ylabel('# individuals')
title('No quarantine, changing rates')
legend("Susceptible", "Zombies")
hold off

nexttile
plot(t, p_a(:,4), '-k') %removed
xlabel('Time (hours)')
ylabel('# individuals')
ylim([0, 22000])
title('No quarantine, changing rates')
legend("Removed (dead zombies)")

nexttile
plot(t, p_a(:,5), 'm')%Ksz
xlabel('Time (hours)')
ylabel('Rate (hours-1)')
title('No quarantine, changing rates')
legend("Ksz")

nexttile
plot(t, p_a(:,6), 'Color', [0, 0.75, 0.75]) %Kzr
xlabel('Time (hours)')
ylabel('Rate (hours-1)')
title('No quarantine, changing rates')
legend("Kzr")

%Fraction of total population, semilogy
nexttile
semilogy(t, p_a(:,1)./N, '-b') %suceptible
hold on
semilogy(t, p_a(:,4)./N, '-k') %removed
semilogy(t, p_a(:,3)./N, '-r') %zombies
xlabel('Time (hours)')
ylabel('Fraction of population')
title('Changing rates')
legend("Susceptible", "Zombies", "Removed (dead zombies)")
hold off

%% Quarantine, changing rates
%Basic zombie outbreak model with quarantine of healthy individuals and
%increased skill at zombie avoidance/killing over time

y0 = [S_o, D_o, Z_o, R_o, Q_o, Ksz, Kzr]; %initial values
k = [Vb, Ksd, Kdz, Kqs, Ksq, Alpha, Beta]; %rates

[t_q, p_qa] = ode15s(@zombies_changing_rates_quarantine, tspan, y0 , options, k);

%Linear plots
tiledlayout(3,2);
nexttile
plot(t_q, p_qa(:,1), 'b') %suceptible
hold on
plot(t_q, p_qa(:,3), 'r') %zombies
plot(t_q, p_qa(:,5), '-g')%quarantined
xlabel('Time (hours)')
ylabel('# individuals')
title('Quarantine, changing rates')
legend("Susceptible", "Zombies", "Quarantined")
hold off

nexttile
plot(t_q, p_qa(:,4), 'k') %removed
hold on 
% plot(t_q, p_qa(:,2), 'c')%naturally dead
xlabel('Time (hours)')
ylabel('# individuals')
title('Quarantine, changing rates')
legend("Removed (dead zombies)")
hold off

nexttile
plot(t, p_qa(:,6), 'm')%Ksz
xlabel('Time (hours)')
ylabel('Rate (hours-1)')
title('Quarantine, changing rates')
legend("Ksz")

nexttile
plot(t, p_qa(:,7), 'Color', [0, 0.75, 0.75]) %Kzr
xlabel('Time (hours)')
ylabel('Rate (hours-1)')
title('Quarantine, changing rates')
legend("Kzr")

%Fraction of population plots, semilogy
%Assume little change in N over timescale
nexttile
semilogy(t_q, p_qa(:,1)./N, 'b') %suceptible
hold on
semilogy(t_q, p_qa(:,3)./N, 'r') %zombies
semilogy(t_q, p_qa(:,5)./N, '-g') %quarantined
semilogy(t_q, p_qa(:,4)./N, 'k') %removed
% semilogy(t_q, p_qa(:,2)./N, 'c') %naturally dead
xlabel('Time (hours)')
ylabel('Fraction of population')
title('Quarantine, changing rates')
legend("Susceptible", "Zombies", "Quarantined", "Removed (dead zombies)")
hold off

%% Functions

% No quarantine
function dydt = zombies(~,y,k)
% System of ODEs describing a zombie outbreak without quarantine

% S is susceptible population, Z is zombie population, R is removed population
% Inputs:
% y: a 1x4 [S, D, Z, R]
% k: a 1x5 [Vb, Ksd, Ksz, Kdz, Kzr] 
%                            Vb  - birth rate
%                            Ksd - background natural death rate
%                            Ksz - rate of collision/zombification of suceptible
%                            Kdz - zombification rate of naturally dead
%                            Kzr - zombie destruction rate
%                            
% Output (output is dydt):
% A 1X4  [dS/dt, dD/dt, dZ/dt, dR/dt]
 
dydt = zeros(4,1);
S = y(1);
D = y(2);
Z = y(3);
R = y(4);

Vb = k(1);
Ksd = k(2);
Ksz = k(3);
Kdz = k(4);
Kzr = k(5);

dydt(1,:) = Vb - Ksd.*S - Ksz.*S.*Z./(S+Z);
dydt(2,:) = Ksd.*S - Kdz.*D;
dydt(3,:) = Ksz.*S.*Z./(S+Z) + Kdz.*D - Kzr.*S.*Z./(S+Z);
dydt(4,:) = Kzr.*Z.*S./(S+Z);

end

% Quarantine of healthy individuals
function dydt = zombies_quarantine(~,y,k)
% System of ODEs describing a zombie outbreak with quarantine

% S is susceptible population, Z is zombie population, R is removed population
% Inputs:
% y: a 1x5 [S, D, Z, R, Q]
% k: a 1x7 [Kb, Ksd, Ksz, Kdz, Kzr, Kqs, Ksq] 
%                            Vb  - birth rate
%                            Ksd - background natural death rate
%                            Ksz - rate of collision/zombification of suceptible
%                            Kdz - zombification rate of naturally dead
%                            Kzr - zombie destruction rate
%                            Kqs - rate of leaving quarantine
%                            Ksq - rate of entering quarantine
%                            
% Output (output is dydt):
% A 1X5  [dS/dt, dD/dt dZ/dt, dR/dt, dQ/dt]
 

dydt = zeros(4,1);
S = y(1);
D = y(2);
Z = y(3);
R = y(4);
Q = y(5);

Vb = k(1);
Ksd = k(2);
Ksz = k(3);
Kdz = k(4);
Kzr = k(5);
Kqs = k(6);
Ksq = k(7);

dydt(1,:) = Vb - Ksd.*S - Ksz.*S.*Z./(S+Z+Q) - Ksq.*S + Kqs.*Q;
dydt(2,:) = Ksd.*S - Kdz.*D;
dydt(3,:) = Ksz.*S.*Z./(S+Z+Q) + Kdz.*D - Kzr.*S.*Z./(S+Z+Q);
dydt(4,:) = Kzr.*Z.*S./(S+Z+Q);
dydt(5,:) = Ksq.*S - Kqs.*Q;

end 

%No quarantine, with changing rates
function dydt = zombies_changing_rates(~,y,k)
% System of ODEs describing a zombie outbreak with no quarantine and
% increased skill at killing/avoiding zombies over time

% S is susceptible population, Z is zombie population, R is removed population
% Inputs:
% y: a 1x6 [S, D, Z, R, Ksz, Kzr]
% k: a 1x5 [Vb, Ksd, Kdz, Alpha, Beta] 
%                            Vb  - birth rate
%                            Ksd - background natural death rate
%                            Ksz - rate of collision/zombification of suceptible
%                            Kdz - zombification rate of naturally dead
%                            Kzr - zombie destruction rate
%                            Alpha - factor of decrease in Ksz
%                            Beta - factor of increase in Kzr
%                            
% Output (output is dydt):
% A 1X6  [dS/dt, dD/dt, dZ/dt, dR/dt, dKsz/dt, dKzr/dt]
 
dydt = zeros(6,1);
S = y(1);
D = y(2);
Z = y(3);
R = y(4);
Ksz_mod = y(5);
Kzr_mod = y(6);

Vb = k(1);
Ksd = k(2);
Kdz = k(3);
Alpha = k(4);
Beta = k(5);

dydt(1,:) = Vb - Ksd.*S - Ksz_mod.*S.*Z./(S+Z);
dydt(2,:) = Ksd.*S - Kdz.*D;
dydt(3,:) = Ksz_mod.*S.*Z./(S+Z) + Kdz.*D - Kzr_mod.*S.*Z./(S+Z);
dydt(4,:) = Kzr_mod.*Z.*S./(S+Z);
dydt(5,:) = -Alpha.*Ksz_mod;
dydt(6,:) = Beta.*Kzr_mod;

end

%Quarantine, with changing rates
function dydt = zombies_changing_rates_quarantine(~,y,k)
% System of ODEs describing a zombie outbreak with quarantine and
% increased skill at killing/avoiding zombies over time
% S is susceptible population, Z is zombie population, R is removed
% population, Q is quarantined population

% Inputs:
% y: a 1x7 [S, D, Z, R, Q, Ksz, Kzr]
% k: a 1x7 [Vb, Ksd, Kdz, Kqs, Ksq, Alpha, Beta] 
%                            Vb  - birth rate
%                            Ksd - background natural death rate
%                            Ksz - rate of collision/zombification of suceptible
%                            Kdz - zombification rate of naturally dead
%                            Kzr - zombie destruction rate
%                            Kqs - rate of leaving quarantine
%                            Ksq - rate of entering quarantine
%                            Alpha - factor of decrease in Ksz
%                            Beta - factor of increase in Kzr
%                            
% Output (output is dydt):
% A 1X7  [dS/dt, dD/dt, dZ/dt, dR/dt, dQ/dt, dKsz/dt, dKzr/dt]
 
dydt = zeros(7,1);
S = y(1);
D = y(2);
Z = y(3);
R = y(4);
Q = y(5);
Ksz_mod = y(6);
Kzr_mod = y(7);

Vb = k(1);
Ksd = k(2);
Kdz = k(3);
Kqs = k(4);
Ksq = k(5);
Alpha = k(6);
Beta = k(7);

dydt(1,:) = Vb - Ksd.*S - Ksz_mod.*S.*Z./(S+Z+Q) - Ksq.*S + Kqs.*Q;
dydt(2,:) = Ksd.*S - Kdz.*D;
dydt(3,:) = Ksz_mod.*S.*Z./(S+Z+Q) + Kdz.*D - Kzr_mod.*S.*Z./(S+Z+Q);
dydt(4,:) = Kzr_mod.*Z.*S./(S+Z+Q);
dydt(5,:) = Ksq.*S - Kqs.*Q;
dydt(6,:) = -Alpha.*Ksz_mod;
dydt(7,:) = Beta.*Kzr_mod;

end