%% Simulation eines Antriebsstrangmodells in MATLAB/Simulink
% V2: Zur thermischen Analyse einzelner Betriebspunkte n_EM/M_EM.

clear all;                                                                  %#ok<CLALL>
close all;
clc;

%% I. Setup des Antriebsstrangmodells

%% Allgemeine Paramter

%Pfad der Simulationsparameter und Funktionen hinzufügen
addpath(strcat(pwd,'\Simulationsparameter_Komponenten'));
addpath(strcat(pwd,'\Funktionen_Plot'));

%--------------------------------------------------------------------------
% Paramter Modell
%--------------------------------------------------------------------------

Modell='Antriebsstrangmodell_V2';                          % Name des Modells
open_system(Modell);                                                        % Oeffnen des Modells

%--------------------------------------------------------------------------
% Parameter Subsystem "Input"
%--------------------------------------------------------------------------
%%% Hier Eingaben tätigen für Umgebungsbedingungen, EM-Betriebspunkte (n,M), Belastungsdauer
%%% (t_Simulation), Schrittweite der Drehzahl und des Drehmoments und Temperaturgrenzen von Rotor und Stator.
T_Umgebung=273.15 + 40;                    % Umgebungstemperatur in K 
T_Start_EMaschine = 273.15 + 40;           % Starttemperatur E-Maschine in K
T_Start_Coolant = 273.15 + 40;             % Starttemperatur Kühlmittel in K
%n_EM = 5000;                               % min^-1
%M_EM = 45;                                 % Nm
t_Simulation = 240*60;                      % s
T_Stator_Wicklung_Limit = 273.15 + 120;    % K

n_min = 500;                               % min^-1
n_max = 6000;                              % min^-1
n_Stuetzstellen = 12;                      % -
M_Diskretisierung = 1;                    % Nm

% n_EM = [0,n_EM;t_Simulation,n_EM];
% M_EM = [0,M_EM;t_Simulation,M_EM];
n_Sim = linspace(n_min,n_max,n_Stuetzstellen);
M_Sim = zeros(1,length(n_Sim));

%--------------------------------------------------------------------------
% Parameter Subsystem "Output"
%--------------------------------------------------------------------------

Datenfrequenz_Output_max=100;                                               % Maximale Datenfrequenz der Outputs in 1/s -> Maximale Anzahl an Datenpunkten der Simulationsergebnisse, die mittels "To Workspace"-Bloecken nach Ende der Simulation von Simulink an MATLAB uebergeben werden, pro Sekunde
% -> Datenfrequenz_Output_max legt die maximale Anzahl an Datenpunkten der Simulationsergebnisse, die mittels "To Workspace"-Bloecken nach Ende der Simulation von Simulink an MATLAB uebergeben werden, pro Sekunde global fuer jeden Output fest
% -> Es handelt sich nicht um einen genauen Wert sondern um einen Maximalwert, weil die tatsaechliche Anzahl an Datenpunkten der Simulationsergebnisse, die mittels "To Workspace"-Bloecken nach Ende der Simulation von Simulink an MATLAB uebergeben werden, pro Sekunde zudem vom Zeitschritt der Simulation abhaengt
% -> Es ist zu beachten, dass sich Datenfrequenz_Output_max nur auf die Simulationsergebnisse, die mittels "To Workspace"-Bloecken nach Ende der Simulation von Simulink an MATLAB uebergeben werden, bezieht
% -> Das Antriebsstrangmodell in Simulink wird unabh. von Datenfrequenz_Output_max zu jedem einzelnen Auswertungszeitpunkt ausgewertet

%--------------------------------------------------------------------------
% Parameter Zeit
%--------------------------------------------------------------------------
dt_max_Anwender=1;                                                       % Durch den Anwender festgelegter maximaler Zeitschritt der Simulation in s -> Je kleiner der Zeitschritt der Simulation ist, desto besser ist die Ergebnisqualitaet der Simulation, desto hoeher ist jedoch der Simulationsaufwand
% -> Fuer eine numerisch stabile Simulation gibt es noch andere Bedingungen, die den maximalen Zeitschritt der Simulation festlegen
% -> Fuer jede dieser Bedingungen wird ein maximaler Zeitschritt der Simulation bestimmt
% -> Aus Gruenden des Simulationsaufwandes ist ein moeglichst großer Zeitschritt der Simualtion von Vorteil 
% -> Deswegen wird fuer den Zeitschritt der Simualtion zuerst der minimale aller vorhandenen maximalen Zeitschritte der Simulation ausgewaehlt
% -> Fuer den wirklichen Zeitschritt der Simulation wird nur 99,9% des minimalen aller vorhandenen maximalen Zeitschritte der Simulation verwendet, um Rundungsfehler zu umgehen und somit sicher fuer numerische Stabilitaet der diskreten Simulation zu sorgen

%% Parameter Subsystem "E-Maschine Verlustleistung"

% -> Im verwendeten Antriebstrangmodell wird genau mit 1 E-Maschine gerechnet

%--------------------------------------------------------------------------
% Definition allgemeiner Parameter
%--------------------------------------------------------------------------
Verlustleistung_Fit = 1; % Fit-Parameter um Temperaturverläufe aus Messung anzunähern, mit ihm wird Gesamtverlustleistung multipliziert

load('KennfeldMaschine_NEmo_V5.mat'); %NEmo                                           % Laden des Wirkungsgradkennfeldes der E-Maschine

% Wikrungsgradkennfeld Antrieb
eta_EMaschine_Break_M_EMaschine=flipud(KennfeldMaschine.M(:,1))';                     % Der Wirkungsgrad der E-Maschine ist abh. vom Drehmoment der E-Maschine in Nm
eta_EMaschine_Break_n_EMaschine=KennfeldMaschine.n(end-1,:);                          % Der Wirkungsgrad der E-Maschine ist abh. von der Drehzahl der E-Maschine 1/min
eta_EMaschine_Table=flip(KennfeldMaschine.etages,1);                                  % Der Wirkungsgrad der E-Maschine in Abh. vom Drehmoment der E-Maschine und in Abh. von der Drehzahl der E-Maschine

% Wirkungsgradkennfeld Rekuperation
eta_EMaschine_Break_M_EMaschine_reku=-KennfeldMaschine.M(:,1)';
eta_EMaschine_Break_n_EMaschine_reku=KennfeldMaschine.n(end-1,:);
eta_EMaschine_Table_reku=KennfeldMaschine.etages;

clearvars KennfeldMaschine;

eta_EMaschine_Faktor_Volllast=0.96;                                         % Faktor der Volllast fuer das Wirkungsgradkennfeld der E-Maschine
% -> Außerhalb der Volllastkennlinie ist das Wirkungsgradkennfeld mit 'NaN' bedatet
% -> Sollten in der Simualtion Betriebspunkte außerhalb des Kennfelds liegen, muss die Simulation noch immer funktionieren
% -> Es wird fuer diese Betriebspunkte dann der Wirkungsgrad am Rand des Kennfelds benutzt bzw. genauer bei eta_EMaschine_Faktor_Volllast * Volllast
% -> Dies ist notwendig, da sowohl Wirkungsgradkennfeld wie auch Lookup Table nicht stetig sind und somit zwischen Stuetzstellen interpolieren muessen
% -> Dies kann dazu fuehren, dass ein Wert, der auf der Volllastlinie sitzt, in 'NaN' im Wirkungsgradkennfeld resultiert
% -> Um dies zu vermeiden, wird mit eta_EMaschine_Faktor_Volllast * Volllast gerechnet, wenn Betriebspunkte außerhalb des Kennfelds liegen, um 'NaN' Wirkungsgrade zu vermeiden
% -> eta_EMaschine_Faktor_Volllast ist stark abh. von der Guete (Abstand der Stuetzstellen) des Wirkungsgradkennfeldes
% -> Fuer das aktuell verwendete Wirkungsgradkennfeld ***.mat darf maximal mit 98% der Volllast gerechnet werden, um 'NaN' Wirkungsgrade zu vermeiden
% -> Deswegen ist eta_EMaschine_Faktor_Volllast zu 0.99 festgelegt
%--------------------------------------------------------------------------
% Bestimmung allgemeiner Paramater
%--------------------------------------------------------------------------
% Maximales Drehmoment beim Antreiben
M_EMaschine_max=zeros(2,numel(eta_EMaschine_Break_n_EMaschine));            % Maximales Drehmoment der E-Maschine nach Schoetz "Thermische Modellierung und Optimierung elektrischer Antriebsstraenge"
i=1;
while i<=numel(eta_EMaschine_Break_n_EMaschine)
    j=2;
    while  j<=numel(eta_EMaschine_Break_M_EMaschine) && isnan(eta_EMaschine_Table(j,i))==0
        j=j+1;        
    end
    M_EMaschine_max(1,i)=eta_EMaschine_Break_n_EMaschine(i);
    M_EMaschine_max(2,i)=eta_EMaschine_Break_M_EMaschine(j-1);
    i=i+1;
end 
% Maximales Rekuperationsmoment ist im Modell als feste Kurve hinterlegt

% P_EMaschine_max=[M_EMaschine_max(1,:);(M_EMaschine_max(1,:).*M_EMaschine_max(2,:)*pi/30)]; % Maximale Leistung der E-Maschine = 2 * pi * maximales Drehmoment der E-Maschine * Drehzahl der E-Maschine / 60 mit Drehzal in 1/min nach Schoetz "Thermische Modellierung und Optimierung elektrischer Antriebsstraenge"

%% Parameter Subsystem "E-Maschine Verlustteilung"
% Gleiche Breakpoints wie vom Wirkungsgradkennfeld genutzt
%eta_EMaschine_Break_M_EMaschine=flipud(KennfeldMaschine.M(:,1))';
%eta_EMaschine_Break_n_EMaschine=KennfeldMaschine.n(end-1,:);

load('Ser_Verlustteilung_KF_V5.mat') % NEmo                                 % Laden der Verlustteilungskennfelder der E-Maschine

% Antrieb
share_EMaschine_Wicklung_Table = ser.mot.share_Wicklungsverluste;
share_EMaschine_Ummagnet_Table = ser.mot.share_Ummagnetisierungsverluste;
share_EMaschine_WirkZus_Table = ser.mot.share_WirkZusVerl;
share_EMaschine_WirkMech_Table = ser.mot.share_WirkmechVerluste;

% Rekuperation
share_EMaschine_Wicklung_reku_Table = ser.gen.share_Wicklungsverluste;
share_EMaschine_Ummagnet_reku_Table = ser.gen.share_Ummagnetisierungsverluste;
share_EMaschine_WirkZus_reku_Table = ser.gen.share_WirkZusVerl;
share_EMaschine_WirkMech_reku_Table = ser.gen.share_WirkmechVerluste;

%% Parameter Subsystem "E-Maschine Thermik"
%--------------------------------------------------------------------------
% Bestimmung Einzelkomponenten des thermischen Netzwerkes
%--------------------------------------------------------------------------
% Massen der EM-Komponenten und zugehörige spezifische Wärmekapazitäten,
% sowie Flächen für Konvektion in E-Maschine. Zusätzliche Annahmen zur
% Maschinengeometrie sowie Materialkennwerte können im Skript vorgegeben
% werden.
create_Emotor_parameters; % Skript, welches Parameter aus Entwurfsparameter der E-Maschinenauslegung von Svenja berechnet
% -> dieses Skript greift auf die Datei "Entwurf_ASM_20190704_164434.mat" zurueck, welche im Ordner "Simulationsparameter_Komponenten" abgelegt sein muss. Sie enthält die Entwurfsparameter der EM-Auslegung.

%--------------------------------------------------------------------------
% Definition Parameter für Luft-Konvektion in der E-Maschine
%--------------------------------------------------------------------------
% Physikal. Eigenschaften Luft bei T=60°C
lam_Luft = 0.028592; % W/(m*K)
cp_Luft = 1008.8; % J/(kg*K)
rho_Luft = 1.0468; % kg/m^3
nue_Luft = 19.230034390524e-6; % m^2/s

%--------------------------------------------------------------------------
% Definition Parameter für temperaturabhängige ohmsche Verluste in Rotor und Stator
%--------------------------------------------------------------------------
% Eingabe durch den Nutzer:
alpha_Cu = 3.9e-3; % K^-1; Temperaturkoeffizient des elektrischen Widerstands für Kupfer; Quelle: https://www.chemie.de/lexikon/Temperaturkoeffizient.html
alpha_Al = 4.0e-3; % K^-1; Temperaturkoeffizient des elektrischen Widerstands für Aluminium; Quelle: https://www.chemie.de/lexikon/Temperaturkoeffizient.html

% Bestimmung der Parameter aus der Eingabe und den Entwurfsparametern der E-Maschine
load('Entwurf_ASM_20190704_164434.mat') % Entwurfsparameter der E-Maschine aus Svenjas Berechnungstool
T_Entwurf_Stator = 273.15+Entwurf.Optionen.theta_1; % K; Temperatur Leitermaterial Stator
T_Entwurf_Rotor = 273.15+Entwurf.Optionen.theta_2; % K; Temperatur Leitermaterial Rotor
if strncmp(Entwurf.Optionen.Stator_Leitermaterial.Bezeichnung,'Kupfer',3) == 1
    alpha_Stator = alpha_Cu;
elseif strncmp(Entwurf.Optionen.Stator_Leitermaterial.Bezeichnung,'Alu',3) == 1
    alpha_Stator = alpha_Al;
else
    disp('Das Stator-Leitermaterial ist nicht bekannt. Es kann kein Temperaturkoeffizient zugeordnet werden.')
    keyboard;
end
if strncmp(Entwurf.Optionen.Rotor_Leitermaterial.Bezeichnung,'Kupfer',3) == 1
    alpha_Rotor = alpha_Cu;
elseif strncmp(Entwurf.Optionen.Rotor_Leitermaterial.Bezeichnung,'Alu',3) == 1
    alpha_Rotor = alpha_Al;
else
    disp('Das Rotor-Leitermaterial ist nicht bekannt. Es kann kein Temperaturkoeffizient zugeordnet werden.')
    keyboard;
end
clearvars Entwurf

%--------------------------------------------------------------------------
% Bestimmung Parameter der Wasserkühlung
%--------------------------------------------------------------------------
% Diese können offline (d.h. vor der Simulation) berechnet werden, da keine
% Abhängigkeit von Betriebsparametern vorliegt.
% Geometrie des Kühlkanals an der E-Maschine wurde bereits in Skript
% "create_motor_parameters" definiert. (siehe vorheriger Absatz)

% --> Kennwerte des Kühlmittels: hier Wasser/Glykol 50%/50%
lam_WC = 0.416; % W*(m*K)^-1; spezifische Wärmeleitfähigkeit
nue_WC = 0.00000397; % m^2*s^-1; kinematische Viskosität
cp_WC = 3322.5; % J*(kg*K)^-1; spezifische Wärmekapazität
rho_WC = 1065.8; % kg*m^-3; Dichte Kühlmittel

% --> Kennwerte des Kühlkreislaufes
vpkt_WC = 10; % l*min^-1; Kühlmittelvolumenstrom, wenn Kühlmittelpumpe angeschaltet
m_WC_pumpon = 9; % kg; Kühlmittelmenge im Kreislauf
radiator_power = 50; % W/K; Kühlleistung des Kreislaufes bzgl. der E-Maschine
% Starttemperatur Kühlmittelpumpe
pump_start_limit = 273.15 + 33; % K; 33 °C aus NEmo-Daten abgelesen; Wickelkopftemperatur, ab der die Pumpe beginnt zu arbeiten

% --> Berechnung der Kennwerte
% Fall 1: Wärmeübergangskoeffizient, wenn Kühlmittelpumpe angeschaltet
dh_WC = 4*A_WC/p_WC; % m; hydraulischer Durchmesser
v_WC = vpkt_WC/(1e-3*60)/A_WC; % m*s^-1; Fluidgeschwindigkeit
Re_WC = dh_WC*v_WC/nue_WC; % -; Reynold-Zahl
Pr_WC = rho_WC*cp_WC*nue_WC/lam_WC; % -; Prandtl-Zahl
if Re_WC < 2300 % laminare Strömung im Kühlmantel
    Nu_WC = (3.66^3+1.615^3*Re_WC*Pr_WC*dh_WC/l_WC)^(1/3); % -; Nusselt-Zahl
elseif Re_WC >= 2300 && Pr_WC < 1.5 % turbulente Strömun (Fall 1)
    Nu_WC = 0.0214*(Re_WC^0.8-100)*Pr_WC^0.4*(1+((dh_WC/l_WC)^2)^(1/3)); % -; Nusselt-Zahl
elseif Re_WC >= 2300 && Pr_WC >= 1.5 % turbulente Strömun (Fall 2)
    Nu_WC = 0.0120*(Re_WC^0.87-280)*Pr_WC^0.4*(1+((dh_WC/l_WC)^2)^(1/3)); % -; Nusselt-Zahl
end
h_WC_pumpon = Nu_WC*lam_WC/dh_WC; % W*(m^2*K)^-1; Wärmeübertragungskoeffizient Kühlwasser im Kühlmantel
% Fall 2: Wärmeübergangskoeffizient, wenn Kühlmittelpumpe ausgeschaltet
v_WC = 0;
Re_WC = 0;
Nu_WC = 3.66; % -; Nusselt-Zahl
h_WC_pumpoff = Nu_WC*lam_WC/dh_WC; % W*(m^2*K)^-1; Wärmeübertragungskoeffizient Kühlwasser im Kühlmantel
% Masse des Wasser im Kühlmantel (benötigt, wenn Kühlmittelpumpe ausgeschaltet)
m_WC = A_WC*l_WC*rho_WC; % kg

clearvars dh_WC A_WC p_WC lam_EC nue_WC rho_WC vpkt_WC v_WC Re_WC Pr_WC l_WC Nu_WC

%% III. Bestimmung relevanter Parameter des Modells
%% Allgemeine Paramter

%--------------------------------------------------------------------------
% Bestimmung des Zeitschrittes der Simulation
%--------------------------------------------------------------------------
dt=0.999*min(dt_max_Anwender); % Verwendeter Zeitschritt der Simulation in s

%--------------------------------------------------------------------------
% Bestimmung des Parameters "Decimation" in den "To Workspace"-Bloecken im Subsystem "Output"
%--------------------------------------------------------------------------

% -> Der Parameter "Decimation" in einem "To Workspace"-Block gibt an, zu jedem wievielten Auswertungszeitpunkt die Simulationsergebnisse waehrend der Simulation in Simulink zwischengespeichert und nach Ende der Simulation an MATLAB uebergeben werden (-> Zeitschritte der Simulation pro Datenpunkt der Simulationsergebnisse)
% -> Der Parameter "Decimation" in einem "To Workspace"-Block muss eine natuerliche Zahl groeßer 0 sein
% -> Dieser Parameter wird ueber Decimation_Output fuer jeden "To Workspace"-Block im Subsystem "Output" global bestimmt
% -> Decimation_Output ist abh. von Datenfrequenz_Output_max und von dt und wird ueber folgenden mathematischen Zusammenhang bestimmt:
% -> Decimation_Output=ceil(1/(Datenfrequenz_Output_max*dt))
% -> Da der Parameter "Decimation" in einem "To Workspace"-Block eine natuerliche Zahl groeßer 0 sein muss, wird bei der Bestimmung von Decimation_Output aufgerundet, um immer eine natuerliche Zahl und immer mindestens einen Wert = 1 zu erhalten
% -> Wenn bei der Bestimmung von Decimation_Output ein Aufrundvorgang stattfindet (Ergebnis der Division ist keine natuerliche Zahl), unterscheidet sich die tatsaechliche Datenfrequenz der Outputs von Datenfrequenz_Output_max, die vom Anwender festgelegt wird
% -> Durch Aufrunden wird die tatsaechliche Datenfrequenz der Outputs kleiner als Datenfrequenz_Output_max
% -> Deswegen handelt es sich bei der Angabe von Datenfrequenz_Output_max korrekterweise um einen Maximalwert
Decimation_Output=ceil(1/(Datenfrequenz_Output_max*dt));

%% V. Allgemeine Einstellungen zum Vermeiden von Warnungen bei der Simulation des Thermomanagemensystemmodells

set_param(Modell,'UnconnectedInputMsg','none');                             % Unverbundene Eingaenge sollen keine Warnung ausgeben -> Diese koennen wegen dem allgemein modellierten Thermomanagementsystemmodell vorhanden sein, fuehren aber zu keinem Fehler, weil die entsprechenden unverbundenen Eingaenge fuer die Simulation nicht benoetigt werden
set_param(Modell,'UnconnectedOutputMsg','none');                            % Unverbundene Ausgaenge sollen keine Warnung ausgeben -> Diese koennen wegen dem allgemein modellierten Thermomanagementsystemmodell vorhanden sein, fuehren aber zu keinem Fehler, weil die entsprechenden unverbundenen Ausgaenge fuer die Simulation nicht benoetigt werden

%% VI. Simulation des Thermomanagementsystemmodells

keyboard;

for c=1:length(n_Sim)
    n_EM = [0,n_Sim(c);t_Simulation,n_Sim(c)];
    
    M_EM = interp1(M_EMaschine_max(1,:), M_EMaschine_max(2,:), n_Sim(c))-1;
    %test
    %M_EM = M_EM-M_Diskretisierung;
    %ende
    M_EM = [0,M_EM;t_Simulation,M_EM];
    sim(Modell);
    while max(EMaschine_Thermik_Output.T_EMaschine_Wickelkopf__K_.Data) > T_Stator_Wicklung_Limit
        disp('nochmal')
        M_EM = M_EM-[0,M_Diskretisierung;0,M_Diskretisierung];
        sim(Modell);
    end
    M_Sim(c) = M_EM(1,2);
    disp('Erfolg!')
    Darstellung_der_Simulationsergebnisse_V2;  
end

figure
plot(n_Sim,M_Sim)