close all;
clear; 

%Developer: Walter, Doerich, Mandic

%veränderbare Variablen
T       = 0.005;    %obere Grenze für t (Betrachtungszeit (<= 0.005))
n       = 100;      %Anzahl der Summanden pro Zeitschritt
xMax    = 100;      %Anzahl Wegschritte
%abhängige Variablen (bzw. vorgegeben aus Aufgabenstellung)
K       = 1;                        %Faktor vor Uxx
L       = 1;                        %obere Intervallgrenze für x
x       = linspace(0, L, xMax);     %Unterteilung Weg
deltaX  = L/xMax;                   %Abstand zwischen 2 Weg-Schritten
deltaT  = (0.5*deltaX^2)/K;         %Abstand zwischen 2 Zeit-Schritten
tMax    = ceil(T/deltaT);           %Anzahl Zeitschritte
t       = linspace(0, T, tMax);     %Unterteilung Zeit
d       = K*deltaT/deltaX^2;        %Faktor für explizite numerische Lösung 
%Bausteine für Teilaufgaben i und ii 
B_nI    = zeros(n, 1);
B_nII   = zeros(n, 1);
%exakte Lösung für Teilaufgaben i und ii
uI      = zeros(xMax, tMax);
uII     = zeros(xMax, tMax);
%numerisch explizite Lösung für Teilaufgaben i und ii
vI      = zeros(xMax, tMax);
vII     = zeros(xMax, tMax);
%%%%%%%%%%%%%%%%%
%Funktionen für exakte Lösung
funE    = @(n, t) exp(-((n*pi/L)^2)*K*t);
funCos  = @(n, x) cos((n*pi*x)/L); 
%Funktionen für Anfangsbedingungen
funI    = @(x) 9+3*cos(pi*x) + 5*cos(4*pi*x);
funII   = @(x) x;

%Anfangswerte
for j = 1:xMax
    uI(j, 1)    = funI(x(j));
    uII(j, 1)   = funII(x(j));
    vI(j, 1)    = funI(x(j));
    vII(j, 1)   = funII(x(j));   
end

%Bausteine berechnen nach Neumann
for i = 1:n
    B_nI(i, 1)      = (2/L) * integral(@(x) funI(x).*funCos(i, x),0,L);
    B_nII(i, 1)     = (2/L) * integral(@(x) funII(x).*funCos(i, x),0,L);
end

%Exakte Lösung nach Neumann
for j = 1:xMax
    for k = 2:tMax-1
        tempI   = 0;
        tempII  = 0;
        for i = 1:n     %Summenwert (für n viele Summanden)
            tempI       = tempI + (B_nI(i, 1)*funE(i, t(k))*funCos(i, x(j)));
            tempII      = tempII + (B_nII(i, 1)*funE(i, t(k))*funCos(i, x(j)));
        end
        uI(j,k)     = tempI + 0.5*(2/L) * integral(@(x) funI(x),0,L);
        uII(j,k)    = tempII + 0.5*(2/L) * integral(@(x) funII(x),0,L);
    end
end

%Numerische Lösung (explizit)
for k = 2:tMax-1
    %Geisterzellen (Anfang)
    vI(1,k)     = d*(vI(2, k-1) - 2*vI(1, k-1) + vI(2, k-1)) + vI(1, k-1);
    vII(1,k)    = d*(vII(2, k-1) - 2*vII(1, k-1) + vII(2, k-1)) + vII(1, k-1);
    
    %Geisterzellen (Ende)
    vI(xMax,k)     = d*(vI(xMax-1, k-1) - 2*vI(xMax, k-1) + vI(xMax-1, k-1)) + vI(xMax, k-1);
    vII(xMax,k)    = d*(vII(xMax-1, k-1) - 2*vII(xMax, k-1) + vII(xMax-1, k-1)) + vII(xMax, k-1);
    for j = 2:xMax-1
        vI(j,k)     = d*(vI(j+1, k-1) - 2*vI(j, k-1) + vI(j-1, k-1)) + vI(j, k-1);
        vII(j,k)    = d*(vII(j+1, k-1) - 2*vII(j, k-1) + vII(j-1, k-1)) + vII(j, k-1);
    end
end

%Beim impliziten Verfahren führen wesentlich weniger Zeitschritte zu einem
%nahezu gleichwertigen Ergebnis.
deltaT_IM   = 10*(deltaX^2)/K;          %Abstand zwischen 2 Zeit-Schritten
tMax_IM     = ceil(T/deltaT_IM);        %Anzahl Zeitschritte
t_IM        = linspace(0, T, tMax_IM);  %Unterteilung Zeit
d_IM        = K*deltaT_IM/deltaX^2;     %Faktor für explizite numerische Lösung 

%numerisch implizite Lösung für Teilaufgaben i und ii
wI      = zeros(xMax, tMax_IM);
wII     = zeros(xMax, tMax_IM);

%Numerische Lösung (implizit) 
%Tridiagonalmatrix A erstellen
vec1    = zeros(xMax-2, 1);
vec2    = zeros(xMax-3, 1);

for j = 1:xMax
    vec1(j, 1) = 1+2*d_IM;
end
for j = 1:xMax-1
    vec2(j, 1) = -d_IM;
end
A1 = diag(vec1);
A2 = diag(vec2, 1);
A3 = diag(vec2, -1);
A = A1 + A2 + A3;

%R-Vektor erstellen und befüllen für Teilaufgaben i und ii
RI = zeros(xMax, 1);
RII = zeros(xMax, 1);
for i = 1:xMax-1
  RI(i, 1) = funI(x(i+1));
  RII(i, 1) = funII(x(i+1));
end

%Lösung berechnen für i-ten R-Vektor, diese Lösung wird neuer r-Vektor im
%i+1-ten Schritt
for i = 1:tMax_IM-1
    xVecI = A\RI(:,i); %x-Vector, wird neuer r-Vektor
    RI = [RI, xVecI];
    xVecII = A\RII(:,i); %x-Vector, wird neuer r-Vektor
    RII = [RII, xVecII];
end

%Speichern der impliziten numerischen Lösung
for i = 2:xMax-1
    wI(i,:) = RI(i-1,:);
    wII(i,:) = RII(i-1,:);
end

%%%Plots
%Plot i)
figure('Name', 'exakt und implizite numerische Lösung (i)','NumberTitle','off')
plot(x, uI(1:xMax,1),'b')
xlabel('Betrachtetes Objekt','FontAngle','italic');
ylabel('Temperaturverteilung','FontAngle','italic');
hold on 
i1 = plot(x, vI(1:xMax,1),'b+');

plot(x, uI(1:xMax,  ceil(tMax/6)),'r')
i2 = plot(x, wI(1:xMax,  ceil(tMax_IM/6)),'ro');

plot(x, uI(1:xMax, 2* ceil(tMax/6)),'m')
i3 = plot(x, wI(1:xMax, 2* ceil(tMax_IM/6)),'m*');

plot(x, uI(1:xMax, 3* ceil(tMax/6)),'k')
i4 = plot(x, wI(1:xMax, 3* ceil(tMax_IM/6)),'k.');

plot(x, uI(1:xMax, 4* ceil(tMax/6)),'g')
i5 = plot(x, wI(1:xMax, 4* ceil(tMax_IM/6)),'gd');

plot(x, uI(1:xMax, 5* ceil(tMax/6)),'y')
i6 = plot(x, wI(1:xMax, 5* ceil(tMax_IM/6)),'ys');
lgd = legend([i1 i2 i3 i4 i5 i6],'0 Sekunden','0.00083 Sekunden','0.00166 Sekunden','0.00249 Sekunden','0.00332 Sekunden','0.00415 Sekunden','0,005 Sekunden');
hold off

%Plot ii)
figure('Name', 'exakt und implizite numerische Lösung (ii)','NumberTitle','off')
plot(x, uII(1:xMax,1),'b')
xlabel('Betrachtetes Objekt','FontAngle','italic');
ylabel('Temperaturverteilung','FontAngle','italic');
hold on 
i1 = plot(x, vII(1:xMax,1),'b+');

plot(x, uII(1:xMax,  ceil(tMax/6)),'r')
i2 = plot(x, wII(1:xMax,  ceil(tMax_IM/6)),'ro');

plot(x, uII(1:xMax, 2* ceil(tMax/6)),'m')
i3 = plot(x, wII(1:xMax, 2* ceil(tMax_IM/6)),'m*');

plot(x, uII(1:xMax, 3* ceil(tMax/6)),'k')
i4 = plot(x, wII(1:xMax, 3* ceil(tMax_IM/6)),'k.');

plot(x, uII(1:xMax, 4* ceil(tMax/6)),'g')
i5 = plot(x, wII(1:xMax, 4* ceil(tMax_IM/6)),'gd');

plot(x, uII(1:xMax, 5* ceil(tMax/6)),'y')
i6 = plot(x, wII(1:xMax, 5* ceil(tMax_IM/6)),'ys');
lgd = legend([i1 i2 i3 i4 i5 i6],'0 Sekunden','0.00083 Sekunden','0.00166 Sekunden','0.00249 Sekunden','0.00332 Sekunden','0.00415 Sekunden','0,005 Sekunden');
hold off

%Plot i)
figure('Name', 'exakt und explizite numerische Lösung (i)','NumberTitle','off')
plot(x, uI(1:xMax,1),'b')
xlabel('Betrachtetes Objekt','FontAngle','italic');
ylabel('Temperaturverteilung','FontAngle','italic');
hold on 
i1 = plot(x, vI(1:xMax,1),'b+');

plot(x, uI(1:xMax,  ceil(tMax/6)),'r')
i2 = plot(x, vI(1:xMax,  ceil(tMax/6)),'ro');

plot(x, uI(1:xMax, 2* ceil(tMax/6)),'m')
i3 = plot(x, vI(1:xMax, 2* ceil(tMax/6)),'m*');

plot(x, uI(1:xMax, 3* ceil(tMax/6)),'k')
i4 = plot(x, vI(1:xMax, 3* ceil(tMax/6)),'k.');

plot(x, uI(1:xMax, 4* ceil(tMax/6)),'g')
i5 = plot(x, vI(1:xMax, 4* ceil(tMax/6)),'gd');

plot(x, uI(1:xMax, 5* ceil(tMax/6)),'y')
i6 = plot(x, vI(1:xMax, 5* ceil(tMax/6)),'ys');
lgd = legend([i1 i2 i3 i4 i5 i6],'0 Sekunden','0.00083 Sekunden','0.00166 Sekunden','0.00249 Sekunden','0.00332 Sekunden','0.00415 Sekunden','0,005 Sekunden');
hold off

%Plot ii)
figure('Name', 'exakt und explizite numerische Lösung (ii)','NumberTitle','off')
plot(x, uII(1:xMax,1),'b')
xlabel('Betrachtetes Objekt','FontAngle','italic');
ylabel('Temperaturverteilung','FontAngle','italic');
hold on 
i1 = plot(x, vII(1:xMax,1),'b+');

plot(x, uII(1:xMax,  ceil(tMax/6)),'r')
i2 = plot(x, vII(1:xMax,  ceil(tMax/6)),'ro');

plot(x, uII(1:xMax, 2* ceil(tMax/6)),'m')
i3 = plot(x, vII(1:xMax, 2* ceil(tMax/6)),'m*');

plot(x, uII(1:xMax, 3* ceil(tMax/6)),'k')
i4 = plot(x, vII(1:xMax, 3* ceil(tMax/6)),'k.');

plot(x, uII(1:xMax, 4* ceil(tMax/6)),'g')
i5 = plot(x, vII(1:xMax, 4* ceil(tMax/6)),'gd');

plot(x, uII(1:xMax, 5* ceil(tMax/6)),'y')
i6 = plot(x, vII(1:xMax, 5* ceil(tMax/6)),'ys');
lgd = legend([i1 i2 i3 i4 i5 i6],'0 Sekunden','0.00083 Sekunden','0.00166 Sekunden','0.00249 Sekunden','0.00332 Sekunden','0.00415 Sekunden','0,005 Sekunden');
hold off