close all;
clear;

%Developer: Walter, Doerich, Mandic

%veränderbare Variablen
T       = 0.1;                  %obere Grenze für t (Betrachtungszeit)
n       = 100;                  %Anzahl der Summanden pro Zeitschritt
eps     = 0;                    %Epsilon (aus Aufgabenstellung)
a       = 1;                    %variabel (aus Aufgabenstellung)
alpha   = -((9*pi^2)/4) - a^2;  %Konstante aus Separationsansatz (X * T = alpha)
xMax    = 100;                  %Anzahl Wegschritte
%abhängige Variablen (bzw. vorgegeben aus Aufgabenstellung)
L       = 1;                                                %obere Intervallgrenze für x
x       = linspace(0, L, xMax);                             %Unterteilung Weg
deltaX  = L/xMax;                                           %Abstand zwischen 2 Weg-Schritten
deltaT  = min(0.9*deltaX/(2*a), (0.5*0.5*deltaX^2)/eps);    %Abstand zwischen 2 Zeit-Schritten
tMax    = ceil(T/deltaT);                                   %Anzahl Zeitschritte
t       = linspace(0, T, tMax);                             %Unterteilung Zeit
u       = zeros(xMax, tMax);                                %exakte Lösung
v       = zeros(xMax, tMax);                                %numerische Lösung
%%%%%%%%%%%%%%%%%
%Funktionen für exakte Lösung
funSin  = @(x) sin((3*pi*x)/2); 
funEx   = @(x, t) exp(a*x + alpha*t);
%Funktionen für exakte Lösung
funE    = @(x) exp(a*x);

%Exakte Lösung
for k = 1:tMax
    for j = 1:xMax
        u(j, k) = funEx(x(j), t(k))*funSin(x(j));
    end
end

%Anfangswerte
for i = 1:xMax
   v(i, 1) = funE(x(i))*funSin(x(i));
end
%Randwerte
for i = 1:tMax
   v(1, i) = 0;
end

%Numerische Löung (explizit)
for k = 2:tMax-1
    if (a < 0)
        v(xMax, k) = 2*a*deltaX*v(xMax-1,k-1)+v(xMax-2,k-1);
        for j = 2:xMax-1   
            temp1   = (v(j+1, k-1) - 2*v(j, k-1) + v(j-1, k-1))/deltaX;
            temp2   = v(j+1, k-1) - v(j, k-1);
            v(j, k) = (eps*temp1 - 2*a*temp2)*(deltaT/deltaX) + v(j, k-1);
        end
    end
    if (a >= 0)
        v(xMax, k) = 2*a*deltaX*v(xMax-1,k-1)+v(xMax-2,k-1);
        for j = 2:xMax-1   
            temp1   = (v(j+1, k-1) - 2*v(j, k-1) + v(j-1, k-1))/deltaX;
            temp2   = v(j, k-1) - v(j-1, k-1);
            v(j, k) = (eps*temp1 - 2*a*temp2)*(deltaT/deltaX) + v(j, k-1);
        end
    end
end

%Plot
figure('Name', 'exakte Lösung und explizite numerische Lösung für a=1 und Epsilon=0','NumberTitle','off')
plot(x, u(1:xMax,1),'b')
xlabel('Betrachtetes Objekt','FontAngle','italic');
ylabel('Temperaturverteilung','FontAngle','italic');
hold on 
i1 = plot(x, v(1:xMax,1),'b+');

plot(x, u(1:xMax,  ceil(tMax/6)),'r')
i2 = plot(x, v(1:xMax,  ceil(tMax/6)),'ro');

plot(x, u(1:xMax, 2* ceil(tMax/6)),'m')
i3 = plot(x, v(1:xMax, 2* ceil(tMax/6)),'m*');

plot(x, u(1:xMax, 3* ceil(tMax/6)),'k')
i4 = plot(x, v(1:xMax, 3* ceil(tMax/6)),'k.');

plot(x, u(1:xMax, 4* ceil(tMax/6)),'g')
i5 = plot(x, v(1:xMax, 4* ceil(tMax/6)),'gd');

plot(x, u(1:xMax, 5* ceil(tMax/6)),'y')
i6 = plot(x, v(1:xMax, 5* ceil(tMax/6)),'ys');
lgd = legend([i1 i2 i3 i4 i5 i6],'0 Sekunden','0.02 Sekunden','0.04 Sekunden','0.06 Sekunden','0.08 Sekunden','0.1 Sekunden');
hold off