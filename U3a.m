% Daniel Willi 13.11.2017
clear % Zuerst wird der Workspace von alten Variablen befreit
cla, clc % Die Grafiken und die Konsole werden zur kgesetzt

%% Laden der Daten
load U3Dat1

%% Parameter
deg = 3; % Hier kann die Modellordnung variiert werden
s0_2_priori = 950^2; % Varianz a priori
alpha = 0.05;

%% Preprozessierung
t0 = data(1,1); % t0 wird bestimmt ...
data(:,1) = data(:,1) - t0; % ... und von den Daten abgezogen
n = length(data);

%% Uebersichtsplot
plot(data(:,1)+t0,data(:,2)/1000,'.k')
xlabel('Jahr')
ylabel('Gesundheitskosten [Milliarden CHF]')

%% Stochastisches Modell
Qll = eye(n);

%% Funktionales Modell
A = zeros(n,deg+1);
for k = 1:deg+1
    A(:,k) = data(:,1).^(k-1);
end

%% Ausgleichung
% Weil die Modelle linear sind, braucht es keine Iterationen
Qll_1 = inv(Qll);
Qxx = inv(A'*Qll_1*A); %#ok<*MINV>
dx = Qxx*A'*Qll_1*data(:,2);

%% Residuen
l_est = zeros(size(i)); % l_est sind die geschaetzen Messwerte
for k = 1:deg+1
    l_est = l_est + dx(k)*data(:,1).^(k-1);
end
% Die Schaetzung wird von den Messwerten abgezogen, es bleiben die Residuen
v = data(:,2) - l_est;
% Die Varianz a posteriori ist die gewichtete Quadratsumme der Residuen
%       geteilt durch die Freiheitsgrade
s0_2_post = (v'*Qll_1*v)/(n-length(dx));

%% Zweiter Teil der Uebung - Praediktion
% Erstellen einer Zeitreihe 50 Jahre zurueck und 50 Jahre in die Zukunft
i = (data(1,1):data(end,1)+50)';
% i = (data(1,1):data(end,1))';

F = zeros(size(i,1),deg+1); % F initialisieren
l_est = zeros(size(i)); % l_est initialisieren
for k = 1:deg+1
    l_est = l_est + dx(k)*i.^(k-1); % Berechnen der Modellwerte
    F(:,k) = i.^(k-1); % Berechnen der Ableitungen fuer dieselben Punkte
end

%% Varianzfortpflanzung
Kxx = s0_2_post*Qxx;
Kyy = F*Kxx*F';

% Die Standardabweichung der Messungen ist die Wurzel der entsprechenden
%       Diagonalelemente
s_yy = diag(Kyy).^(0.5);

%% Plots
hold on
col = get(groot,'DefaultAxesColorOrder');
p1 = plot(i+t0,l_est/1000,'Color',col(1,:));
p2 = plot(i+t0,(l_est+3*s_yy)/1000,'Color',col(3,:));
plot(i+t0,(l_est-3*s_yy)/1000,'Color',col(3,:));

legend([p1,p2],'Schaetzung','\pm 3 \sigma',...
    'Location','northwest')

%% Quotiententest
A_lim = finv(1-alpha/2,n-length(dx),inf);
B_lim = finv(alpha/2,n-length(dx),inf);
fprintf('Fisher Intervall mit alpha = %u%s:\n%.2f - %.2f\n'...
    ,100*alpha,'%',B_lim,A_lim)
fprintf('Test-Statistik: %.2f\n\n', s0_2_post/s0_2_priori)

%% Ausgabe der Parameter mit Genauigkeiten
for k = 1:length(dx)
    fprintf('a_%u = %.3f +/- %.4f MCHF/y-%u\n',k-1,dx(k),sqrt(Kxx(k,k)),k-1)
end


