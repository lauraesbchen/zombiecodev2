%20.320 Project 2
%Modeling Train to Busan zombie outbreak
%2020 11 02

%% Nonquarantine
%Initial values
S_o = 51249999; %#people, population of South Korea
R_o = 0; %#no dead zombies initially
D_o = 0; %#no naturally dead people initially
Z_o = 1; %#one zombie at t=0
N = S_o + Z_o;%total population

%Rate constants
Vb  = 0; %assume birth rate zero on short timescale 
Ksd = 0; %assume death rate zero on short timescale
Ksz = 10/24; %hours-1
Kdz = 0; %set to zero, as no D population
Kzr = 0.00008333333333;%hours-1

y0 = [S_o, D_o, Z_o, R_o];%initial values
k = [Vb,Ksd, Ksz, Kdz, Kzr]; %rate constants

tspan = linspace(0,120); %hours

options = odeset('AbsTol', 1e-12, 'RelTol', 5e-12);

[t, p] = ode23s(@zombies, tspan, y0 , options, k);

%Linear plots
figure
plot(t, p(:,1), '-b') %Suceptible
hold on
plot(t, p(:,3), '-r') %Zombies
xlabel('Time (hours)')
ylabel('# people')
title('No quarantine')
legend(["Susceptible", "Zombies"])
hold off

figure
plot(t, p(:,4), '-k') %Removed
hold on
plot(t, p(:,2), '-c') %Dead
xlabel('Time (hours)')
ylabel('# people')
legend('Removed (killed zombies)', 'Naturally dead')
title('No quarantine')
hold off

%Fraction of population plots, semilogy
figure()
semilogy(t, p(:,1)./N, '-b') %Suceptible
hold on
semilogy(t, p(:,3)./N, '-r') %Zombies
semilogy(t, p(:,4)./N, '-k') %Removed
semilogy(t, p(:,2)./N, '-c') %Dead
xlabel('Time (hours)')
ylabel('Fraction of population')
title('No quarantine')
legend('Susceptible', 'Zombies','Removed (dead zombies)','Naturally dead')
hold off
%% Non quarantine changing rates
%Initial values
S_o = 51249999; %#people, population of South Korea
R_o = 0; %#no dead zombies initially
D_o = 0; %#no naturally dead people initially
Z_o = 1; %#one zombie at t=0
N = S_o + Z_o;%total population

%Rate constants
Vb  = 0; 
Ksd = 0; 
Ksz = 10/24; %hours-1
Kdz = 0; 
Kzr = 0.00008333333333;%hours-1
Alpha = 0.01;%accounting for increased skill at killing/avoiding zombies

tspan = 60*60; %minutes

options = odeset('AbsTol', 1e-6, 'RelTol', 5e-6);

[t, p_dr] = ode23s(@zombies_non_SIR, [0, tspan], [S_o, D_o, Z_o, R_o, Ksz], options, [Vb, ...
    Ksd, Kdz, Kzr, Alpha]);

figure()
loglog(t, p_dr(:,1), '-b')
hold on
loglog(t, p_dr(:,2), '-c')
hold on
legend(["Susceptible", "Dead"])
hold off
figure()
loglog(t, p_dr(:,3), '-r')
hold on
loglog(t, p_dr(:,4), '-k')
legend(["Zombies", "Removed"])


%% Quarantine 
%Initial values
S_o = 41249999;%#people, non-quarantined population of South Korea
R_o = 0;
D_o = 0;
Z_o = 1; %one zombie at t=0
Q_o = 10E6; %#people, estimated initial quarantined population
N = S_o + Q_o + Z_o; %total population

%Rate constants
Vb  = 0; %assume birth rate zero on short timescale 
Ksd = 0; %assume death rate zero on short timescale
Ksz = 10/24; %hours-1
Kdz = 0; %set to zero, as no D population
Kzr = 0.00008333333333;%hours-1
Kqs = 0; %assume no leaving quarantine on short timescale
Ksq = 1/120; %hours-1, one person enters quaratine every 5 days on average

y0 = [S_o, D_o, Z_o, R_o, Q_o]; %initial values
k = [Vb,Ksd, Ksz, Kdz, Kzr, Kqs, Ksq]; %rate constants

tspan = linspace(0,120); %hours

options = odeset('AbsTol', 1e-14, 'RelTol', 2.3e-14);

[t_q, p_q] = ode15s(@zombies_quarantine, tspan, y0 , options, k);

%Linear plots
figure()
plot(t_q, p_q(:,1), 'b') %suceptible
hold on
plot(t_q, p_q(:,3), 'r') %zombies
plot(t_q, p_q(:,5), '-g') %quarantined
xlabel('Time (hours)')
ylabel('# people')
title('Quarantine')
legend("Susceptible", "Zombies", "Quarantined")
hold off

figure
plot(t_q, p_q(:,4), 'k') %removed
hold on 
plot(t_q, p_q(:,2), 'c') %naturally dead
xlabel('Time (hours)')
ylabel('# people')
title('Quarantine')
legend("Removed (killed zombies)", "Naturally dead")
hold off

%Fraction of population plots, semilogy
figure
semilogy(t_q, p_q(:,1)./N, 'b') %suceptible
hold on
semilogy(t_q, p_q(:,3)./N, 'r') %zombies
semilogy(t_q, p_q(:,5)./N, '-g') %quarantined
semilogy(t_q, p_q(:,4)./N, 'k') %removed
semilogy(t_q, p_q(:,2)./N, 'c') %naturally dead
xlabel('Time (hours)')
ylabel('Fraction of population')
title('Quarantine')
legend("Susceptible", "Zombies", "Quarantined", "Removed (killed zombies)", "Naturally dead")
hold off
%% Functions
% Non-Quarantine of Zombies
function dydt = zombies(~,y,k)

% 
% System of ODEs describing a zombie outbreak
%%S is susceptible population, Z is zombie population, R is removed population
% Inputs:
% y: a 1x4 [S, D, Z, R]
% k: a 1x5 [Vb, Ksd, Ksz, Kdz, Kzr] 
%                            Vb  - Birth rate
%                            Ksd - background death rate
%                            Ksz - bite rate
%                            Kdz - zombification rate
%                            Kzr - "zombie destruction" rate
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

%assuming no birth rate
dydt(1,:) = Vb - Ksd.*S - Ksz.*S.*Z./(S+Z);
dydt(2,:) = Ksd.*S - Kdz.*D;
dydt(3,:) = Ksz.*S.*Z./(S+Z) + Kdz.*D - Kzr.*S.*Z./(S+Z);
dydt(4,:) = Kzr.*Z.*S./(S+Z);

end

% Quarantine People from Zombies

function dydt = zombies_quarantine(~,y,k)

% 
% System of ODEs describing a zombie outbreak
%%S is susceptible population, Z is zombie population, R is removed population
% Inputs:
% y: a 1x5 [S, D, Z, R, Q]
% k: a 1x7 [Kb, Ksd, Ksz, Kdz, Kzr, Kqs, Ksq] 
%                            Vb  - Birth rate
%                            Ksd - background death rate
%                            Ksz - bite rate
%                            Kdz - zombification rate
%                            Kzr - "zombie destruction" rate
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

%assuming no birth rate
dydt(1,:) = Vb - Ksd.*S - Ksz.*S.*Z./(S+Z+Q) - Ksq.*S + Kqs.*Q;
dydt(2,:) = Ksd.*S - Kdz.*D;
dydt(3,:) = Ksz.*S.*Z./(S+Z+Q) + Kdz.*D - Kzr.*S.*Z./(S+Z+Q);
dydt(4,:) = Kzr.*Z.*S./(S+Z+Q);
dydt(5,:) = Ksq.*S - Kqs.*Q;

end 

function dydt = zombies_changing_rates(~,y,k)

% 
% System of ODEs describing a zombie outbreak
%%S is susceptible population, Z is zombie population, R is removed population
% Inputs:
% y: a 1x4 [S, D, Z, R, Ksz]
% k: a 1x4 [Vb, Ksd, Kdz, Kzr, Alpha] 
%                            Vb  - Birth rate
%                            Ksd - background death rate
%                            Ksz - bite rate
%                            Kdz - zombification rate
%                            Kzr - "zombie destruction" rate
%                            
% Output (output is dydt):
% A 1X4  [dS/dt, dD/dt, dZ/dt, dR/dt, dKsz/dt]
 

dydt = zeros(5,1);
S = y(1);
D = y(2);
Z = y(3);
R = y(4);
Ksz_mod = y(5);

Vb = k(1);
Ksd = k(2);
Kdz = k(3);
Kzr = k(4);
Alpha = k(5);

%assuming no birth rate
dydt(1,:) = Vb - Ksd.*S - Ksz_mod.*S.*Z./(S+Z);
dydt(2,:) = Ksd.*S - Kdz.*D;
dydt(3,:) = Ksz_mod.*S.*Z./(S+Z) + Kdz.*D - Kzr.*S.*Z./(S+Z);
dydt(4,:) = Kzr.*Z.*S./(S+Z);
dydt(5,:) = -Alpha.*Ksz_mod;


end
