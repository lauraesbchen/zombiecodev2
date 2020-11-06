%20.320 Project 2
%Modeling Train to Busan zombie outbreak
%2020 11 06

%RUN THIS SECTION FIRST
close all
clear all
clc

%Initial values
S_o = 51249999;     %suceptible; #people; population of South Korea in 2016
R_o = 0;            %removed (dead zombies)
D_o = 0;            %naturally dead
Z_o = 1;            %zombies
Q_o = 0;            %quarantined population
N = S_o + Q_o + Z_o;%total population

%Rate constants
Kc = 13.4/24;      %average rate of collision, individuals/hour/individual = hour-1
Vb  = 0;           %birth rate; zero on short timescale 
Ksd = 0;           %natural death rate; zero on short timescale
Ksz = Kc*8/10;     %hours-1, rate of collision x probability of bite
Kdz = 0;           %zombification rate of naturally dead; set to zero, as no D population
Kzr = Kc*1/5000;   %hours-1; zombie destruction rate; rate of collision x probability someone can kill a zombie
Kqs = 0;           %rate of exiting quarantine; assume no leaving quarantine on short timescale
Ksq = 1/120;       %hours-1, rate of entering quarantine

Re = Ksz/Kzr;      %effective reproduction number (R0 with no intervention)

%Multiplication factors for changing rate constants
Alpha = 0.006;     %accounting for increased skill at avoiding zombies over time, Ksz
Beta = 0.006;      %accounting for increased skill at killing zombies over time, Kzr

tspan = linspace(0,120); %hours

options = odeset('AbsTol', 1e-8, 'RelTol', 1e-8);%set accuracy 

%% No Quarantine
%Basic zombie outbreak model with no quarantine of healthy individuals

y0 = [S_o, D_o, Z_o, R_o];%initial values
k = [Vb,Ksd, Ksz, Kdz, Kzr]; %rate constants

[t, p] = ode15s(@zombies, tspan, y0 , options, k);

%Linear plots
figure
plot(t, p(:,1), '-b') %suceptible
hold on
plot(t, p(:,3), '-r') %zombies
xlabel('Time (hours)')
ylabel('# individuals')
title('No Quarantine: S and Z vs t when Z0 = 0')
legend(["Susceptible", "Zombies"])
hold off

figure
plot(t, p(:,4), '-k') %removed
hold on
% plot(t, p(:,2), '-c') %naturally dead
xlabel('Time (hours)')
ylabel('# individuals')
ylim([0,14500])
legend('Removed (dead zombies)')
title('No Quarantine: R vs t')
hold off

%Fraction of population plots, semilogy
%Assume little change in N over timescale
figure
semilogy(t, p(:,1)./N, '-b') %suceptible
hold on
semilogy(t, p(:,3)./N, '-r') %zombies
semilogy(t, p(:,4)./N, '-k') %removed
% semilogy(t, p(:,2)./N, '-c') %dead
xlabel('Time (hours)')
ylabel('Fraction of population')
title('No Quarantine: Fraction of Population vs t')
legend('Susceptible', 'Zombies','Removed (dead zombies)')
hold off

%% Quarantine 
%Zombie outbreak model with quarantine of healthy individuals

y0 = [S_o, D_o, Z_o, R_o, Q_o]; %initial values
k = [Vb,Ksd, Ksz, Kdz, Kzr, Kqs, Ksq]; %rate constants

[t_q, p_q] = ode15s(@zombies_quarantine, tspan, y0 , options, k);

%Linear plots
figure
plot(t_q, p_q(:,1), 'b') %suceptible
hold on
plot(t_q, p_q(:,3), 'r') %zombies
plot(t_q, p_q(:,5), '-g')%quarantined
xlabel('Time (hours)')
ylabel('# individuals')
title('Quarantine: S, Z, and Q vs t')
legend("Susceptible", "Zombies", "Quarantined")
hold off

figure
plot(t_q, p_q(:,4), 'k') %removed
% plot(t_q, p_q(:,2), 'c')%naturally dead
xlabel('Time (hours)')
ylabel('# individuals')
ylim([0, 10000]);
title('Quarantine: R vs t')
legend("Removed (dead zombies)")
hold off

%Fraction of population plots, semilogy
%Assume little change in N over timescale
figure
semilogy(t_q, p_q(:,1)./N, 'b') %suceptible
hold on
semilogy(t_q, p_q(:,3)./N, 'r') %zombies
semilogy(t_q, p_q(:,5)./N, '-g') %quarantined
semilogy(t_q, p_q(:,4)./N, 'k') %removed
% semilogy(t_q, p_q(:,2)./N, 'c') %naturally dead
xlabel('Time (hours)')
ylabel('Fraction of population')
title('Quarantine: Fraction of Population vs t')
legend("Susceptible", "Zombies", "Quarantined", "Removed (dead zombies)")
hold off

%% Increased Skill
%Zombie outbreak model with quarantine of healthy individuals and
%increased skill at zombie avoidance/killing over time. Ksz and Kzr change
%with time. 

y0 = [S_o, D_o, Z_o, R_o, Q_o, Ksz, Kzr]; %initial values
k = [Vb, Ksd, Kdz, Kqs, Ksq, Alpha, Beta]; %rates

[t_q, p_qa] = ode15s(@zombies_changing_rates_quarantine, tspan, y0 , options, k);

%Linear plots
figure
plot(t_q, p_qa(:,1), 'b') %suceptible
hold on
plot(t_q, p_qa(:,3), 'r') %zombies
plot(t_q, p_qa(:,5), '-g')%quarantined
xlabel('Time (hours)')
ylabel('# individuals')
title('Increased Skill: S, Z, and Q vs t')
legend("Susceptible", "Zombies", "Quarantined")
hold off

figure
plot(t_q, p_qa(:,4), 'k') %removed
hold on 
% plot(t_q, p_qa(:,2), 'c')%naturally dead
xlabel('Time (hours)')
ylabel('# individuals')
ylim([0, 17500]);
title('Increased Skill: R vs t')
legend("Removed (dead zombies)")
hold off

figure
plot(t_q, p_qa(:,6), 'm')%Ksz
xlabel('Time (hours)')
ylabel('Rate (hours-1)')
title('Increased Skill: Ksz vs t')
legend("Ksz")

figure
plot(t_q, p_qa(:,7), 'Color', [0, 0.75, 0.75]) %Kzr
xlabel('Time (hours)')
ylabel('Rate (hours-1)')
title("Increased Skill: Kzr vs t")
legend("Kzr")

%Fraction of population plots, semilogy
%Assume little change in N over timescale
figure
semilogy(t_q, p_qa(:,1)./N, 'b') %suceptible
hold on
semilogy(t_q, p_qa(:,3)./N, 'r') %zombies
semilogy(t_q, p_qa(:,5)./N, '-g') %quarantined
semilogy(t_q, p_qa(:,4)./N, 'k') %removed
% semilogy(t_q, p_qa(:,2)./N, 'c') %naturally dead
xlabel('Time (hours)')
ylabel('Fraction of population')
title('Increased Skill: Fraction of Population vs t')
legend("Susceptible", "Zombies", "Quarantined", "Removed (dead zombies)")
hold off

%Plot Re over time
figure
plot(t_q, p_qa(:,6)./p_qa(:,7)); %Ksz/Kzr
xlabel('Time (hours)')
ylabel('Re effective reproduction number')
legend("Re")
title("Increased Skill: Re vs t")

%% Functions

% No Quarantine
function dydt = zombies(~,y,k)
% System of ODEs describing a zombie outbreak without quarantine

% S is susceptible population, D is naturally dead population, Z is zombie population, R is removed population
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

N = S+Z;

dydt(1,:) = Vb - Ksd.*S - Ksz.*S.*Z./(N);
dydt(2,:) = Ksd.*S - Kdz.*D;
dydt(3,:) = Ksz.*S.*Z./(N) + Kdz.*D - Kzr.*S.*Z./(N);
dydt(4,:) = Kzr.*Z.*S./(N);

end

% Quarantine of healthy individuals
function dydt = zombies_quarantine(~,y,k)
% System of ODEs describing a zombie outbreak with quarantine

% S is susceptible population, D is naturally dead population, Z is zombie
% population, R is removed population, Q is quarantined population
% Inputs:
% y: a 1x5 [S, D, Z, R, Q]
% k: a 1x7 [Vb, Ksd, Ksz, Kdz, Kzr, Kqs, Ksq] 
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

N = S+Z+Q;

dydt(1,:) = Vb - Ksd.*S - Ksz.*S.*Z./N - Ksq.*S + Kqs.*Q;
dydt(2,:) = Ksd.*S - Kdz.*D;
dydt(3,:) = Ksz.*S.*Z./N + Kdz.*D - Kzr.*S.*Z./N;
dydt(4,:) = Kzr.*Z.*S./N;
dydt(5,:) = Ksq.*S - Kqs.*Q;

end 

%Increased Skill
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

N = S+Z+Q;

dydt(1,:) = Vb - Ksd.*S - Ksz_mod.*S.*Z./N - Ksq.*S + Kqs.*Q;
dydt(2,:) = Ksd.*S - Kdz.*D;
dydt(3,:) = Ksz_mod.*S.*Z./N + Kdz.*D - Kzr_mod.*S.*Z./N;
dydt(4,:) = Kzr_mod.*Z.*S./N;
dydt(5,:) = Ksq.*S - Kqs.*Q;
dydt(6,:) = -Alpha.*Ksz_mod;
dydt(7,:) = Beta.*Kzr_mod;

end