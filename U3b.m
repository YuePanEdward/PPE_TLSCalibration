% Daniel Willi 13.11.2017
%% Standardzeilen
clear
cla, clc

%% Daten laden
load U3Dat2

%% Parameter
y0 = 2.21; % y0
alpha = 0.05;
data(123,:) = []; % Ausreisser entfernen

%% Daten plotten
plot(data(:,1),data(:,2),'.')
xlabel('Zeit [Jahre]')
ylabel('Aktivitaet [kBq]')

%% Vorbereitungen 
% Funktion definieren
rad = @(t,k) y0*exp(-k*t);

% Naeherungswert rechnen
k = log(data(end,2)/y0)/-data(end,1);

% Stochastisches Modell
n = length(data); % Anzahl Beobachtungen
s0_2_priori = (0.012)^2; % A priori Varianz
Qll = eye(n);
P = inv(Qll);

dx = 100; % Initialisieren
it = 0; % Zaehler Iterationen initiieren
%% Ausgleichung
while abs(dx) > 1e-10 && it < 50
    A = zeros(n,1); % A initialisieren
    A(:,1) = -y0*data(:,1).*exp(-k*data(:,1)); % A belegen
    
    % Falls jemand numerisch ableiten mag
    A_num(:,1) = (rad(data(:,1),k+0.005) - rad(data(:,1),k-0.005))/0.01;
    
    % Verkuerzte Beobachtungen
    dl = data(:,2) - rad(data(:,1),k);

    % Ausgleichung
    Qxx = inv(A'*P*A); %#ok<*MINV>
    dx = Qxx*A'*P*dl;
    
    % Update
    it = it + 1; % Zaehler Iterationen
    k = k+dx; % Parameter aktualisieren mit neuem Wert
end


%% Residuen
t = data(:,1); % Vektor mit den Beobachtungszeiten

y_i = rad(t,k); % Ausgeglichene Beobachtungen
v = data(:,2) - y_i; % Residuen
s0_2_post = (v'*P*v)/(n-length(dx)); % Varianz a posteriori

%% Quotiententest
A_lim = finv(1-alpha/2,n-length(dx),inf);
B_lim = finv(alpha/2,n-length(dx),inf);
fprintf('Fisher Intervall mit alpha = %u%s:\n%.2f - %.2f\n'...
    ,100*alpha,'%',B_lim,A_lim)
fprintf('Test-Statistik: %.2f\n\n', s0_2_post/s0_2_priori)

%% Standardisierte Residuen
Qvv = Qll-A*Qxx*A';
v_norm = v./(sqrt(s0_2_priori.*diag(Qvv)));

%% Plots
hold on
plot(t,y_i) % Schaetzwerte ploten (Ausgeglichene Beobachtungen)

figure
plot(v_norm)
title('Standardisierte Residuen')

%% Varianzfortpflanzung zur Bestimmung der Halbwertszeit
T = -log(0.5)*k^(-1); % Schaetzwert

F = log(0.5)*k^(-2); % Design Matrix
Qtt = F*Qxx*F'; % Varianzfortpflanzung

Ktt = s0_2_post*Qtt;
s = sqrt(Ktt); % Stand. Abweichung

fprintf('T(1/2) = %.2f h +/- %.4f h (1 sigma)\n ',T,s)
