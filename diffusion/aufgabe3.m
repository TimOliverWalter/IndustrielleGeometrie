close all;
clear;

%Developer: Walter, Doerich, Mandic

%veränderbare Variablen
T       = 1;        %obere Grenze für t (Betrachtungszeit bis 1 in a), bis 5 in c))
n       = 100;      %Anzahl der Summanden pro Zeitschritt
xMax    = 100;      %Anzahl Wegschritte
%abhängige Variablen (bzw. vorgegeben aus Aufgabenstellung)
a       = 2;                        %Vorgabe aus Aufgabenstellung
K       = 1;                        %Faktor vor Uxx
L       = 1;                        %obere Intervallgrenze für x
x       = linspace(0, L, xMax);     %Unterteilung Weg
deltaX  = L/xMax;                   %Abstand zwischen 2 Weg-Schritten
deltaT  = (0.5*deltaX^2)/K;         %Abstand zwischen 2 Zeit-Schritten
tMax    = ceil(T/deltaT);           %Anzahl Zeitschritte
t       = linspace(0, T, tMax);     %Unterteilung Zeit
d       = K*deltaT/deltaX^2;        %Faktor für explizite numerische Lösung 
B_n     = zeros(n, 1);              %Bausteine
u       = zeros(xMax, tMax);        %exakte Lösung
v       = zeros(xMax, tMax);        %numerische Lösung
%%%%%%%%%%%%%%%%%
%Funktionen für exakte Lösung
funAE   = @(t) a*exp(-t);           %Quellterm
funE    = @(n, t) exp(-(n^2)*(pi^2)*t);
funSin  = @(n, x) sin(n*pi*x);

%Randwerte
for i = 2:tMax
    v(1, i) = 1;
    v(xMax, i) = funAE(t(i));
end

%Numerische Lösung (explizit)
for k = 2:tMax
    for j = 2:xMax-1
        v(j,k)     = d*(v(j+1, k-1) - 2*v(j, k-1) + v(j-1, k-1)) + v(j, k-1) - funAE(t(k))*deltaT;
    end
end

%Bausteine berechnen
for i = 1:n
    B_n(i, 1)      = (-2*(a*(-1)^(i+1)+1))/(i*pi);
end

%Exakte Lösung
for k = 1:tMax
    for j = 1:xMax
        temp  = 0;
        for i = 1:n     %Summenwert (für n viele Summanden)
            temp        = temp + B_n(i, 1)*funE(i, t(k))*funSin(i, x(j));
        end
        u(j,k)     = (funAE(t(k))-1)*x(j) + 1 + temp;
    end
end

%Plot
figure('Name', 'exakt und explizite numerische Lösung','NumberTitle','off')

plot(x, u(1:xMax,  ceil(tMax/6)),'r')
xlabel('Betrachtetes Objekt','FontAngle','italic');
ylabel('Temperaturverteilung','FontAngle','italic');
hold on 
i1 = plot(x, v(1:xMax,  ceil(tMax/6)),'ro');

plot(x, u(1:xMax, 2* ceil(tMax/6)),'m')
i2 = plot(x, v(1:xMax, 2* ceil(tMax/6)),'m*');

plot(x, u(1:xMax, 3* ceil(tMax/6)),'k')
i3 = plot(x, v(1:xMax, 3* ceil(tMax/6)),'k.');

plot(x, u(1:xMax, 4* ceil(tMax/6)),'g')
i4 = plot(x, v(1:xMax, 4* ceil(tMax/6)),'gd');

plot(x, u(1:xMax, 5* ceil(tMax/6)),'y')
i5 = plot(x, v(1:xMax, 5* ceil(tMax/6)),'ys');
if (T == 1)
    lgd = legend([i1 i2 i3 i4 i5],'0.1666 Sekunden','0.3333 Sekunden','0.5 Sekunden','0.6666 Sekunden','0.8333 Sekunden','1 Sekunde');
else
    lgd = legend([i1 i2 i3 i4 i5],'0.8333 Sekunden','1.6666 Sekunden','2.5 Sekunden','3.3333 Sekunden','4.1666 Sekunden','5 Sekunden');
end
hold off
