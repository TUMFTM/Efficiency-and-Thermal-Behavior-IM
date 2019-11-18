% Aktuelle und finale Version V_1_4
% Skript zur iterativen Bestimmung des maximalen Drehmoments im S1-Betrieb
% Für S1-Betrieb muss Schalter in Modell im Getriebe-Block umgelegt werden!

%% Simulation eines Antriebsstrangmodells in MATLAB/Simulink
clear all;                                                                  %#ok<CLALL>
close all;
clc;

%% I. Setup des Antriebsstrangmodells

%% Allgemeine Parameter
%Pfad der Simulationsparameter und Funktionen hinzufügen
addpath(strcat(pwd,'\Simulationsparameter_Komponenten'));
addpath(strcat(pwd,'\Funktionen_Plot'));

%--------------------------------------------------------------------------
% Wahl des Simulationsmodells (Simulink)
%--------------------------------------------------------------------------
Modell='Antriebsstrangmodell_V1_4_n5_WHasym'; % Name des Modells
open_system(Modell); % Oeffnen des Modells

%--------------------------------------------------------------------------
% Parameter Subsystem "Input"
%--------------------------------------------------------------------------
T_Umgebung=273.15 + 23; % Umgebungstemperatur in K 
T_Start_EMaschine = 273.15 + 23; % Wenn fahrzyklus case=1, dann erster Wert aus Messung der Statortemperatur

%--------------------------------------------------------------------------
% Wahl des Fahrzyklus (durch nicht-auskommentieren einer beider Optionen)
%--------------------------------------------------------------------------
fahrzyklus = 0; % Wahl des Zyklus aus Standardfahrzyklen
%fahrzyklus = 1; % Laden von Geschwindigkeitsverlauf aus Messfahrt

switch fahrzyklus
    case 0
        % -> ArtMw130: Artemis Motorway mit maximaler Fahrzeuggeschwindigkeit = 130 km/h
        % -> ArtMw150: Artemis Motorway mit maximaler Fahrzeuggeschwindigkeit = 150 km/h
        % -> ArtRoad: Artemis Rural Road
        % -> ArtUrban: Artemis Urban
        % -> ECE_R15: Urban Driving Cycle (UDC)
        % -> EUDC: Extra Urban Driving Cycle
        % -> FTP72: Federal Test Procedure 72
        % -> FTP75: Federal Test Procedure 75
        % -> HWFET: Highway Fuel Economy Test
        % -> LA92: California Unified Cycle (UC/UCDS)
        % -> LA92_short: Kuerzere Version des LA92
        % -> NEDC: Neuer europaeischer Fahrzyklus (NEFZ)
        % -> NYCC: New York City Cycle
        % -> SC03: SFTP SC03
        % -> US06: SFTP US06
        % -> WLTP_Class2: Worldwide Harmonized Light Vehicles Test Procedure Klasse 2
        % -> WLTP_Class2_Low: Worldwide Harmonized Light Vehicles Test Procedure Klasse 2 Teil 1
        % -> WLTP_Class2_Middle: Worldwide Harmonized Light Vehicles Test Procedure Klasse 2 Teil 2
        % -> WLTP_Class2_High: Worldwide Harmonized Light Vehicles Test Procedure Klasse 2 Teil 3
        % -> WLTP_Class3: Worldwide Harmonized Light Vehicles Test Procedure Klasse 3
        % -> WLTP_Class3_modNEmo: Worldwide Harmonized Light Vehicles Test Procedure Klasse 3, herunterskaliert für vmax=110 km/h
        % -> WLTP_Class3_modNEmo_Depleting: Worldwide Harmonized Light Vehicles Test Procedure Klasse 3, herunterskaliert für vmax=110 km/h, Testzyklus zur Bestimmung des Energieverbrauchs des Fahrzeugs, wenn Reichweite größer 3 WLTP Vollzyklen (3*22,25 km)
                
        Fahrzyklus='WLTP_Class3_modNEmo';                                                   % Name des Fahrzyklus
        v_Fahrzeug=getfield(load('drivingcycles',Fahrzyklus),Fahrzyklus);           % Geschwindigkeitsverlauf des Fahrzyklus in km/h -> Fahrzeuggeschwindigkeit in km/h
        % dreimal hintereinander Zyklus fahren:
        % v_Fahrzeug=timeseries(vertcat(v_Fahrzeug.Data,v_Fahrzeug.Data,v_Fahrzeug.Data),vertcat(v_Fahrzeug.Time,v_Fahrzeug.Time+v_Fahrzeug.Time(end)+1,v_Fahrzeug.Time+v_Fahrzeug.Time(end)*2+2));
        %

    case 1
        create_drivingcycles_NEmo;
end
v_max=110; % km*h^-1; 110 km/h für NEmo; Maximalgeschwindigkeit des Fahrzeugs in km/h (Die Geschwindigkeit des gewählten Zyklus wird auf die angegebene Maximalgeschwindigkeit limitiert)
alpha_Strecke=[0,0;v_Fahrzeug.Time(end),0]; % Steigung der Strecke in %
M_EM=[0,0;v_Fahrzeug.Time(end),0]; % Dummy-Variable; sie wird für S1-Betriebssimulation benötigt (Momentenvorgabe)

%--------------------------------------------------------------------------
% Parameter Zeit - Simulationsdauer nach Fahrzyklusende + Zeitschritt
%--------------------------------------------------------------------------
t_Simulation=v_Fahrzeug.Time(end)+1; % Simulationszeit in s -> Die Simulationszeit wird abh. von der Dauer des Fahrzyklus bzw. von der Dauer des Geschwindigkeitsverlaufs des Fahrzyklus festgelegt
dt_max_Anwender=1; % Durch den Anwender festgelegter maximaler Zeitschritt der Simulation in s -> Je kleiner der Zeitschritt der Simulation ist, desto besser ist die Ergebnisqualitaet der Simulation, desto hoeher ist jedoch der Simulationsaufwand
dt=0.999*dt_max_Anwender; % Verwendeter Zeitschritt der Simulation in s

%--------------------------------------------------------------------------
% Parameter Subsystem "Output" - maximale Datenfrequenz im Output
%--------------------------------------------------------------------------
Datenfrequenz_Output_max=100; % Maximale Datenfrequenz der Outputs in 1/s -> Maximale Anzahl an Datenpunkten der Simulationsergebnisse, die mittels "To Workspace"-Bloecken nach Ende der Simulation von Simulink an MATLAB uebergeben werden, pro Sekunde
Decimation_Output=ceil(1/(Datenfrequenz_Output_max*dt));

%% Parameter Subsystem "Rad"
%--------------------------------------------------------------------------
% Definition allgemeiner Fahrzeugparameter
%--------------------------------------------------------------------------
r_dynamisch=0.2774; % m; NEmo, Dynamischer Radius des Reifens
rho_Luft=1.2; % kg/m^3; Dichte von Luft (-> verwendeter Wert bei 20°C)
A_Stirn=2; %m^2; NEmo ,Stirnflaeche des Fahrzeugs
c_W=0.37; % -; NEmo, Luftwiderstandsbeiwert des Fahrzeugs
m_Fahrzeug=825; % kg; Masse des Fahrzeugs (Leergewicht des leichtesten Smarts fortwo 451: 825 kg)
m_Zuladung=80; % kg; Masse der Zuladung (Fahrer, Messtechnik, ...)
f_R=0.013; % -; (ca. 0.013 nach PISCHINGER) Rollwiderstandsbeiwert der Reifen
e_Fahrzeug=1.0478; % -; Drehmassenzuschlagsfaktor des Fahrzeugs, e_NEmo=1+(0.5*(2*14.7+2*15.6)*r_dynamisch^2+i_Getriebe^2*0.5*17.12*0.054^2)/(r_dynamisch^2*m_Fahrzeug); mit Radgewicht vorne 14.7kg, hinten 15.6 kg; m_Rotor=17.12kg, r_Rotor=0.054m

%--------------------------------------------------------------------------
% Definition Parameter für Simulation eines Allradfahrzeugs
%--------------------------------------------------------------------------
COG_x = 1.4267; % m; Schwerpunktskoordinate in Fahrrichtung (gemessen von Vorderachse); 1.4267 für Tesla Model S
COG_z = 0.4572; % m; Schwerpunktshöhe des Fahrzeugs über der Fahrbahn; 0.4572 für Tesla Model S
WB = 2.96; % m; Radstand des Fahrzeugs; 2.96 für Tesla Model S
% Für Allrad-Simulation müssen Schalter im Simulink-Modell im Rad-Block umgelegt werden

%% Parameter Subsystem "Getriebe"
%--------------------------------------------------------------------------
% Definition Übersetzung + Wirkungsgrad
%--------------------------------------------------------------------------
i_Getriebe=5.697; % -; NEmo=5.697; Uebersetzung des Getriebes, n_EM/n_Rad
eta_Getriebe=0.96; % -; Wirkungsgrad des Getriebes

%% Parameter Subsystem "E-Maschine Verlustleistung"
%--------------------------------------------------------------------------
% Laden des Wirkungsgradkennfelds der E-Maschine aus Datei
%--------------------------------------------------------------------------
load('KennfeldMaschine_NEmo_V5.mat'); % Laden des Wirkungsgradkennfeldes der E-Maschine

%--------------------------------------------------------------------------
% krb-Korrektur - Korrektur der mechanischen Verlustleistung
%--------------------------------------------------------------------------
% Korrektur des experimentellen mechanischen Reibkoeffizienten
% Dieser wurde in der Maschinenberechnung fälschlicher Weise mit 10 Ws^2/m^4 
% angenommen, beträgt bei der wassergekühlten ASM des NEmo jedoch
% ca. 0.53 Ws^2/m^4 wie durch analytische Nachrechnung von Lager- und
% Lüftungsverluste herausgefunden wurde. Daher werden die mechanischen
% Verluste nach der Verlusttteilung um einen Korrekturfaktor verkleinert.
krb_Korrektur = 0.53/10; % 0.53/10 für den NEmo

% Wikrungsgradkennfeld Antrieb
eta_EMaschine_Break_M_EMaschine=flipud(KennfeldMaschine.M(:,1))';                     % Der Wirkungsgrad der E-Maschine ist abh. vom Drehmoment der E-Maschine in Nm
eta_EMaschine_Break_n_EMaschine=KennfeldMaschine.n(end-1,:);                          % Der Wirkungsgrad der E-Maschine ist abh. von der Drehzahl der E-Maschine 1/min
eta_EMaschine_Table=flip(KennfeldMaschine.etages,1);                                  % Der Wirkungsgrad der E-Maschine in Abh. vom Drehmoment der E-Maschine und in Abh. von der Drehzahl der E-Maschine

% Wirkungsgradkennfeld Rekuperation
eta_EMaschine_Break_M_EMaschine_reku=-KennfeldMaschine.M(:,1)';
eta_EMaschine_Break_n_EMaschine_reku=KennfeldMaschine.n(end-1,:);
eta_EMaschine_Table_reku=KennfeldMaschine.etages;

clearvars KennfeldMaschine;

eta_EMaschine_Faktor_Volllast=0.98;                                         % Faktor der Volllast fuer das Wirkungsgradkennfeld der E-Maschine
% -> Außerhalb der Volllastkennlinie ist das Wirkungsgradkennfeld mit 'NaN' bedatet
% -> Sollten in der Simualtion Betriebspunkte außerhalb des Kennfelds liegen, muss die Simulation noch immer funktionieren
% -> Es wird fuer diese Betriebspunkte dann der Wirkungsgrad am Rand des Kennfelds benutzt bzw. genauer bei eta_EMaschine_Faktor_Volllast * Volllast
% -> Dies ist notwendig, da sowohl Wirkungsgradkennfeld wie auch Lookup Table nicht stetig sind und somit zwischen Stuetzstellen interpolieren muessen
% -> Dies kann dazu fuehren, dass ein Wert, der auf der Volllastlinie sitzt, in 'NaN' im Wirkungsgradkennfeld resultiert
% -> Um dies zu vermeiden, wird mit eta_EMaschine_Faktor_Volllast * Volllast gerechnet, wenn Betriebspunkte außerhalb des Kennfelds liegen, um 'NaN' Wirkungsgrade zu vermeiden
% -> eta_EMaschine_Faktor_Volllast ist stark abh. von der Guete (Abstand der Stuetzstellen) des Wirkungsgradkennfeldes
% -> Fuer das aktuell verwendete Wirkungsgradkennfeld ***.mat darf maximal mit 100% der Volllast gerechnet werden, um 'NaN' Wirkungsgrade zu vermeiden
% -> Deswegen ist eta_EMaschine_Faktor_Volllast zu 1 festgelegt

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
%--------------------------------------------------------------------------
% Laden der Verlustteilung der E-Maschine aus Datei
%--------------------------------------------------------------------------
load('Ser_Verlustteilung_KF_V5_korr.mat') % Laden der Verlustteilungskennfelder der E-Maschine

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

% Gleiche Breakpoints für Lookup-Table wie vom Wirkungsgradkennfeld genutzt
% eta_EMaschine_Break_M_EMaschine=flipud(KennfeldMaschine.M(:,1))';
% eta_EMaschine_Break_n_EMaschine=KennfeldMaschine.n(end-1,:);

%% Parameter Subsystem "E-Maschine Thermik"
%--------------------------------------------------------------------------
% Definition Parameter für Luft-Konvektion in der E-Maschine
%--------------------------------------------------------------------------
% Physikal. Eigenschaften Luft bei T=60°C
lam_Luft = 0.028592; %28.08e-3;%0.028592; % W*(m*K)^-1
cp_Luft = 1008.8; %1008;%1008.8; % J*(kg*K)^-1
rho_Luft_EM = 1.0468; %1.078;%1.0468; % kg*m^-3
nue_Luft = 19.230034390524e-6; %18.22;%19.230034390524e-6; % m^2*s^-1
% Parameter für Innenluftkonvektion
eta_IA = 0.5; % -; Lüfterwirkungsgrad der Lüfterkontour auf dem Kurzschlussring (Verhältnis Innenluftgeschwindigkeit zur Umfangsgeschwindigkeit des Kurzschlussringes); 50 % als Startwert vorgeschlagen von MELLOR 1991

%--------------------------------------------------------------------------
% Bestimmung Einzelkomponenten des thermischen Netzwerkes,
% Massen der EM-Komponenten und zugehörige spezifische Wärmekapazitäten und Flächen für Konvektion in E-Maschine,
% Parameter der Maschinengeometrie und Materialkennwerte können im Skript vorgegeben werden!
%--------------------------------------------------------------------------
create_Emotor_parameters_3_n5_WHasym; % Skript, welches Parameter aus Entwurfsparameter der E-Maschinenauslegung von Svenja berechnet; Wert lam_Luft muss übergeben werden

%--------------------------------------------------------------------------
% Aufteilung der mechanischen Verlustleistung auf Lager- und Lüftungsverluste,
% im Skript können Designparameter: mue_Lager & k_AGP abgeändert werden
%--------------------------------------------------------------------------
create_mech_loss_dis; % Skript, welches die Verlustteilung der mechanischen Verluste berechnet

%--------------------------------------------------------------------------
% Definition Parameter für temperaturabhängige Lüftungsverluste
%--------------------------------------------------------------------------
alpha_Luft = -4e-3; % K^-1; Temperaturkoeffizient der Lüftungsverluste (also von Luft): -4 bis -5% pro 10 K Temperaturzunahme nach AUINGER 1999

%--------------------------------------------------------------------------
% Definition Parameter für temperaturabhängige Ummagnetisierungsverluste
%--------------------------------------------------------------------------
alpha_Ummag = -0.75e-3; % K^-1; Temperaturkoeffizient der Ummagnetisierungsverluste: -4 bis -8% pro 80 K Temperaturzunahme nach AUINGER 1999

%--------------------------------------------------------------------------
% Definition Parameter für temperaturabhängige ohmsche Verluste in Rotor und Stator
%--------------------------------------------------------------------------
T_Entwurf_Stator = 273.15+20; % K; Temperatur Leitermaterial Stator in Maschinenberechung
T_Entwurf_Rotor = 273.15+20; % K; Temperatur Leitermaterial Rotor in Maschinenberechung
alpha_Cu = 3.9e-3; % K^-1; Temperaturkoeffizient des elektrischen Widerstands für Kupfer; Quelle: https://www.chemie.de/lexikon/Temperaturkoeffizient.html
alpha_Al = 4.0e-3; % K^-1; Temperaturkoeffizient des elektrischen Widerstands für Aluminium; Quelle: https://www.chemie.de/lexikon/Temperaturkoeffizient.html

% Bestimmung der Parameter aus der Eingabe und den Entwurfsparametern der E-Maschine
load('Entwurf_ASM_20190704_164434.mat') % Entwurfsparameter der E-Maschine aus Svenjas Berechnungstool
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
m_WC_pumpon = 4.2*1e-3*rho_WC; % kg; Kühlmittelmenge im Kreislauf; 4.2 l KM im NEmo (Wert vom Verbrenner-Smart)
radiator_power = 50; % W/K; Kühlleistung des Kreislaufes bzgl. der E-Maschine
% Starttemperatur Kühlmittelpumpe
pump_start_limit = 273.15 + 33; % K; 33 °C aus NEmo-Daten abgelesen

% Berechnung der Kennwerte
% Fall 1: Wärmeübergangskoeffizient, wenn Kühlmittelpumpe angeschaltet
dh_WC = 4*A_WC/p_WC; % m; hydraulischer Durchmesser
v_WC = vpkt_WC*1e-3/60/A_WC; % m*s^-1; Fluidgeschwindigkeit
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

%% II. Simulation des Antriebsstrangmodells
keyboard;

%--------------------------------------------------------------------------
% Eingaben für S1-Simulation
%--------------------------------------------------------------------------
set_param(Modell,'FastRestart','on')       % erlaubt schnelle Iterationen, da nicht jedes mal kompiliert werden muss
t_sim = 3*60*60;                           % s; Belastungsdauer im S1-Betrieb
n_min = 500;                               % min^-1; niedrigste betrachtete Drehzahl
n_max = 6000;                              % min^-1; höchste betrachtete Drehzahl
n_Stuetzstellen = 12;                      % -; Anzahl an Drehzahl-Stützstellen
M_Diskretisierung = 1;                     % Nm; Schrittweite in Iteration des Drehmoments; es wird vom maximalen Drehmoment laut Kennfeld absteigend iteriert, bis die Grenztemperaturen nicht mehr überschritten werden
T_Stator_Wicklung_Limit = 273.15 + 120;    % K; Grenztemperatur der Statorwicklung; 120 °C entspricht Isolationsklasse B

n_Sim = linspace(n_min,n_max,n_Stuetzstellen);
M_Sim = zeros(1,length(n_Sim));
% Hier kann eine erwartete Volllastlinie eingegeben werden, damit die Lösung mit weniger Iterationen gefunden wird.
M_EMaschine_max_mod=[500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000;...
                     109, 109, 109, 109, 109, 102,  90,  80,  70,  60,  53,  50];
for c=1:length(n_Sim)
    v_Fahrzeug=timeseries([0,(2*pi/60*n_Sim(c)/i_Getriebe*r_dynamisch*3.6),(2*pi/60*n_Sim(c)/i_Getriebe*r_dynamisch*3.6)],[0,(n_Sim(c)/200),t_sim]);
    t_Simulation=v_Fahrzeug.Time(end);                                      % Simulationszeit in s -> Die Simulationszeit wird abh. von der Dauer des Fahrzyklus bzw. von der Dauer des Geschwindigkeitsverlaufs des Fahrzyklus festgelegt
    alpha_Strecke=[0,0;v_Fahrzeug.Time(end),0];                                 % Steigung der Strecke in %
    
    % Folgende Zeile auskommentieren, wenn EINE erwartete Volllastlinie eingegeben wurde
    M_EM = interp1(M_EMaschine_max(1,:), M_EMaschine_max(2,:), n_Sim(c))-1;
    % Folgende Zeile auskommentieren, wenn KEINE erwartete Volllastlinie eingegeben wurde
    %M_EM=interp1(M_EMaschine_max_mod(1,:), M_EMaschine_max_mod(2,:), n_Sim(c));
    M_EM = [0,0;(n_Sim(c)/200),M_EM;t_Simulation,M_EM];
    
    Sim=sim(Modell);
    while max(Sim.EMaschine_Thermik_Output.T_EMaschine_WH__K_.T_EMaschine_WH_1__K_.Data) > T_Stator_Wicklung_Limit
        disp('nochmal')
        M_EM = M_EM-[0,0;0,M_Diskretisierung;0,M_Diskretisierung];
        Sim=sim(Modell);
    end
    M_Sim(c) = M_EM(end,end);
    save(fullfile('Ergebnisse',strcat('S1_n',num2str(n_Sim(c)))),'Sim');
    fprintf('n=%.0f min^_{-1}, M_{max}=%.0f Nm\n',n_Sim(c),M_Sim(c))
  end

set_param(Modell,'FastRestart','off')

%% III. Darstellung der Ergebnisse der S1-Simulation
figure
yyaxis left
plot(n_Sim,M_Sim,'x-')
hold on
% Darstellung Kennlinie des maximalen Maschinenmoments (Kennfeldgrenze)
M_EM_max = interp1(M_EMaschine_max(1,:), M_EMaschine_max(2,:), n_Sim)-1;
plot(n_Sim,M_EM_max,'--')
% Darstellung Überlastfaktor
yyaxis right
plot(n_Sim,(M_EM_max./M_Sim),'x-')
ylabel('Überlastfaktor:  M_{EM,max}/M_{EM,S1}  / -')
yyaxis left
xlabel('Drehzahl E-Maschine  n_{EM}  / min^{-1}')
ylabel('max. Drehmoment S1-Betrieb  M_{EM}  / Nm')
xlim([0 6000])
title('max. Drehmoment E-Maschine in S1-Betrieb')
