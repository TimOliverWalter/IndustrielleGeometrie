clear all;
close all;

%Developer: Walter, Doerich, Mandic

% löse die Erhaltungsgleichung: u_t + u_u_x = 0
L = 230;            % Streckenlänge [Meter]
T = 60;             % Betrachtungszeit [Sekunden] 
CFL = 0.9;          % CFL Bedingung nach Literatur (Thomas Sonar)
Imax = 500;         % Anzahl Stützstellen (Wegachse)
VMax = 11.1;        % Maximale Geschwindigkeit (bei RhoMin) [Meter/Sekunde] 
RhoMax = 0.17;      % Maximale Dichte (bei Stau) [Auto/Meter]
RhoStart = 0.095;   % Mittlere Dichte (bei gleichmäßiger Verteilung der Autos auf die Strecke) [Auto/Meter]

x = linspace(0, L, Imax);   % x-Achse aufgeteilt in Imax Stellen
deltaX = L/(Imax - 1);      % (eine) Schrittweite auf x-Achse

a_Rho = ((-2)*RhoMax)*VMax/RhoMax + VMax;   % Strömungsgeschwindigkeit bzw. Wellengeschwindigkeit (am Anfang)
deltaT = CFL*deltaX/abs(a_Rho);             % (eine) Schrittweite auf y-Achse
Nmax = ceil(T/deltaT) + 1;                  % Anzahl Stützstellen (Zeitachse)

u0 = zeros(Imax, 1);        % Start-Verkehrsdichte-Array anlegen
u = zeros(Imax, Nmax);      % Verkehsdichte-Array anlegen
t = linspace(0, T, Nmax);   % t-Achse aufgeteilt in Nmax Stellen
v = zeros(Imax, Nmax);      % Geschwindigkeits-Array anlegen
wg = zeros(Imax, Nmax);     % Wellengeschwindigkeiten-Array anlegen
pos = zeros(Imax, Nmax);    % Positions-Array anlegen

%Setze Anfangsbedingungen in Anlehnung an Literatur (Thomas Sonar)
AnfangDichtestoerung = 80;
EndeDichtestoerung = 200;
freq = 2*pi/(EndeDichtestoerung-AnfangDichtestoerung);  %Frequenz Sinus

%Störung
for i = 1: Imax
    if(x(i)<AnfangDichtestoerung)
        u0(i,1) = RhoStart;
    elseif(x(i)<140)
        u0(i,1)= (1 + 0.3*sin(freq*(x(i) - 80)))* RhoStart; %Schwankung als Störung der Verkehrsdichte
    else
        u0(i,1)= RhoStart;
    end
end

%Anfangs- Verkehrsdichte, Fahrzeuggeschwindigkeit, Wellengeschwindigkeit
%und Position der Fahrzeuge
for i = 1: Imax
   u(i,1) = u0(i,1);
   v(i, 1) = -(VMax/RhoMax)*u0(i, 1) + VMax;
   wg(i, 1) = (((-2)*u0(i, 1)*VMax)/RhoMax + VMax);
   pos(i, 1) = i*deltaX;
end

for n = 2: Nmax
    %Periodische Randbedingungen (Vorlesung 12.11.20)
    %CIR-Verfahren
    a_Rho = (((-2)*u(Imax, n - 1)*VMax)/RhoMax + VMax);
    a_RhoPlus = max(a_Rho, 0);
    
    a_Rho = (((-2)*u(2, n - 1)*VMax)/RhoMax + VMax);
    a_RhoMinus = min(a_Rho, 0);
    
    %Euler-Verfahren mit Resultaten aus CIR-Verfahren
    u(1, n) = u(1, n - 1) - deltaT/deltaX*a_RhoPlus*(u(1, n - 1) - u(Imax, n - 1)) - deltaT/deltaX*a_RhoMinus*(u(2, n - 1) - u(1, n - 1));
    
    v(1, n) = -(VMax/RhoMax)*u(1, n) + VMax;    %zugehörige Fahrzeuggeschwindigkeit
    
    %Abspeichern der zugehörigen Wellengeschwindigkeit
    if (a_RhoPlus ~= 0)
            wg(1, n) = a_RhoPlus;
        elseif (a_RhoMinus ~= 0)
            wg(1, n) = a_RhoMinus;
        else
            wg(1, n) = 0;
    end
    
    %CIR-Verfahren
    a_Rho = (((-2)*u(Imax - 1, n - 1)*VMax)/RhoMax + VMax);   
    a_RhoPlus = max(a_Rho, 0);
        
    a_Rho = (((-2)*u(1, n - 1)*VMax)/RhoMax + VMax);
    a_RhoMinus = min(a_Rho, 0);
    
    %Euler-Verfahren mit Resultaten aus CIR-Verfahren
    u(Imax, n) = u(Imax, n - 1) - deltaT/deltaX*a_RhoPlus*(u(Imax, n - 1) - u(Imax - 1, n - 1)) - deltaT/deltaX*a_RhoMinus*(u(1, n - 1) - u(Imax, n - 1));
    
    v(Imax, n) = -(VMax/RhoMax)*u(Imax, n) + VMax;  %zugehörige Fahrzeuggeschwindigkeit
    
    %Abspeichern der zugehörigen Wellengeschwindigkeit
    if (a_RhoPlus ~= 0)
            wg(Imax, n) = a_RhoPlus;
        elseif (a_RhoMinus ~= 0)
            wg(Imax, n) = a_RhoMinus;
        else
            wg(Imax, n) = 0;
    end
    
    %Berechnung der numerischen Lösung von Rho
    for i = 2: Imax-1
        %CIR-Verfahren
        a_Rho = (((-2)*u(i - 1, n - 1)*VMax)/RhoMax + VMax);
        a_RhoPlus = max(a_Rho, 0);
        
        a_Rho = (((-2)*u(i + 1, n - 1)*VMax)/RhoMax + VMax);
        a_RhoMinus = min(a_Rho, 0);
        
        %Euler-Verfahren mit Resultaten aus CIR-Verfahren
        u(i, n) = u(i, n - 1) - deltaT/deltaX*a_RhoPlus*(u(i, n - 1) - u(i - 1, n - 1)) - deltaT/deltaX*a_RhoMinus*(u(i + 1, n - 1) - u(i, n - 1));
        
        %Deckelung von Rho (zum Verhindern numerischer Ungenauigkeiten)
        u(i, n) = min(u(i, n), RhoMax);
        
        v(i, n) = -(VMax/RhoMax)*u(i, n) + VMax;    %zugehörige Fahrzeuggeschwindigkeit
        
        %Abspeichern der zugehörigen Wellengeschwindigkeit
        if (a_RhoPlus ~= 0)
            wg(i, n) = a_RhoPlus;
        elseif (a_RhoMinus ~= 0)
            wg(i, n) = a_RhoMinus;
        else
            wg(i, n) = 0;
        end
        
    end
end

%Abspeichern der Fahrzeugpositionen
for i = 1 : Imax
    for n = 2 : Nmax
        pos(i, n) = mod((v(i, n)*deltaT + pos(i, n-1)), 230);
    end 
end


%%%%%%%
%PLOTS
%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name', 'Numerische Lösung der Verkehrsdichte','NumberTitle','off')
plot(x, u(1:Imax,1), 'g','LineWidth',2);
axis tight
set(gca,'nextplot','replacechildren');
xlabel('Streckenlänge_{m}','FontAngle','italic');
ylabel('Vekehrsdichte_{A/m}','FontAngle','italic');
for j = 1:Nmax
    plot(x, u(1:Imax,j), 'g','LineWidth',2);
    F(j) = getframe;
end
movie(F,0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name', 'Numerische Lösung der Verkehrsdichte in 6 Schritten','NumberTitle','off')
hold on;
xlabel('Streckenlänge_{m}','FontAngle','italic');
ylabel('Vekehrsdichte_{A/m}','FontAngle','italic');
plot(x, u(1:Imax,1), x, u(1:Imax,  ceil(Nmax/6)),'-.r', x, u(1:Imax, 2* ceil(Nmax/6)), '--m', x, u(1:Imax, 3* ceil(Nmax/6)), ':b', x, u(1:Imax, 4* ceil(Nmax/6)), '-.g', x, u(1:Imax, 5* ceil(Nmax/6)), ':r', x, u(1:Imax, Nmax), '-.m');
lgd = legend('0 Sekunden','8.5 Sekunden','17.1 Sekunden','25.7 Sekunden','34.3 Sekunden','42.9 Sekunden','51.4 Sekunden');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name', 'Strömungsgeschwindigkeit der Stauwelle in 6 Schritten','NumberTitle','off')
hold on;
xlabel('Streckenlänge_{m}','FontAngle','italic');
ylabel('Strömungsgeschwindigkeit_{m/s}','FontAngle','italic');
plot(x, wg(1:Imax, 1), x, wg(1:Imax,  ceil(Nmax/6)),'-.r', x, wg(1:Imax, 2* ceil(Nmax/6)), '--m', x, wg(1:Imax, 3* ceil(Nmax/6)), ':b', x, wg(1:Imax, 4* ceil(Nmax/6)), '-.g', x, wg(1:Imax, 5* ceil(Nmax/6)), ':r', x, wg(1:Imax, Nmax), '-.m');
lgd = legend('0 Sekunden','8.5 Sekunden','17.1 Sekunden','25.7 Sekunden','34.3 Sekunden','42.9 Sekunden','51.4 Sekunden');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name', 'Geschwindigkeiten der Fahrzeuge in 6 Schritten','NumberTitle','off')
hold on;
xlabel('Streckenlänge_{m}','FontAngle','italic');
ylabel('Geschwindigkeiten_{m/s}','FontAngle','italic');
plot(x, v(1:Imax,1), x, v(1:Imax,  ceil(Nmax/6)),'-.r', x, v(1:Imax, 2* ceil(Nmax/6)), '--m', x, v(1:Imax, 3* ceil(Nmax/6)), ':b', x, v(1:Imax, 4* ceil(Nmax/6)), '-.g', x, v(1:Imax, 5* ceil(Nmax/6)), ':r', x, v(1:Imax, Nmax), '-.m');
lgd = legend('0 Sekunden','8.5 Sekunden','17.1 Sekunden','25.7 Sekunden','34.3 Sekunden','42.9 Sekunden','51.4 Sekunden');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name', 'Position der Fahrzeuge auf der Strecke','NumberTitle','off')
hold on;
xlabel('Position_{m}','FontAngle','italic');
ylabel('Zeit_{s}','FontAngle','italic');
plot(pos(1:23:Imax, 1:Nmax), t);