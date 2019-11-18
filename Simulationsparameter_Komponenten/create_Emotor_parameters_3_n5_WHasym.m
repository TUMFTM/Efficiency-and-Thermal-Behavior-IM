% V2: anderer Wärmeleitwert lambda für Wicklung in Statornut
% V3: für Innenluftkonvektion nach Mellor(1991), es werden andere geometrische Angaben dafür benötigt
% n5: Für den Einsatz im thermischen Netzwerk mit n=5 Diskretisierung auf aktiver Länge
% neu ist auch, dass in m_SSL das Gewicht des Vergussmaterials berücksichtigt wird

% Skript zur Berechnung der E-Maschinenparameter,
% welche für die Modellbildung des thermischen Netzwerkes
% benötigt werden.
% Berechnung basierend auf Output der Maschinenauslegung von Svenja.
% Im Bereich 1 können Werte durch Nutzer angepasst werden.
% Im Bereich 2 findet die Berechnung statt.

load('Entwurf_ASM_20190704_164434.mat') % Entwurfsparameter der E-Maschine aus Svenjas Berechnungstool

%% Bereich 1: ==> Hier können Werte abgeändert werden!
%% Externe Annahmen für die Berechnung der thermischen Parameter
% Asymmetrie der Wickelköpfe (bzgl. der Verteilung der Wicklungsmasse)
WHasym = 0.2; % 0 entspricht keiner Asymmetrie; 1 wäre volle Asymmetrie mit vollem Gewicht auf WH1 

% Diskretisierung auf aktiver Länge
n_dis = 5;

% Ergänzende Abmaße zur Berechnung der Maschinengeometrie
LF = 1.22; % -; Längenfaktor: E-Maschinen-Baulänge/aktive Länge
b_HO = 0.006; % m; radiale Gehäusedicke
b_EC = 0.006; % m; Stirnseitige Gehäusedicke
b_SH = 0.004; % m; Wandstärke Rotorwelle
b_BR = 0.01; % m; radiale Lagerbreite
b_WC = 0.005; % m; radiale Breite des Kühlkanals um das Gehäuse der E-Maschine
d_Lack = 0.1; % mm; Dicke des Isolierlackes um die Leiterdrähte der Statorwicklung
A_WH_Faktor = 1.5; % -; Multiplikator für Oberfläche des Wickelkopfes zur Innenluft, dadurch Berücksichtigung einer unregelmäßigen Oberfläche (Mellor 1991)
l_BR_eq = 0.3; % mm; äquivalenter Luftspalt zur Modellierung des Kontaktwiderstandes im Kugellager; empfohlener Startwert nach BOGLIETTI 2007: 0.3 mm

% Dichten zur Berechnung von Massen
% --> Komponenten
rho_HO = 7860; % kg*m^-3; Dichte Gehäuse; hier Stahl
rho_EC = 7860; % kg*m^-3; Dichte Gehäuse-Stirndeckel; hier Stahl
rho_SH = 7860; % kg*m^-3; Dichte Rotorwelle; hier Stahl
rho_Verguss = 1250; % kg*m^-3; Dichte Vergussharz um Statorwicklung in den Statornuten
rho_Lack = rho_Verguss; % kg*m^-3; Dichte Lackisolation um Leiter der Statorwicklung

% Wärmeleitfähigkeiten zur Berechnung thermischer Kenngrößen
% --> Werkstoffe
lam_Cu = 380; % W*(m*K)^-1; Wärmeleitfähigkeit Kupfer
lam_Al = 215; % W*(m*K)^-1; Wärmeleitfähigkeit Aluminium
% --> Komponenten
lam_iso = 0.145; % W*(m*K)^-1; Wärmeleitfähigkeit Nutisolierung
lam_ES = 30; % W*(m*K)^-1; Wärmeleitfähigkeit Elektroblech; hier VACOFLUX 50
lam_HO = 42; % W*(m*K)^-1; Wärmeleitfähigkeit Gehäuse; hier Stahl
lam_EC = 42; % W*(m*K)^-1; Wärmeleitfähigkeit Gehäuse-Stirndeckel; hier Stahl
lam_SH = 42; % W*(m*K)^-1; Wärmeleitfähigkeit Rotorwelle; hier Stahl
% Spezial: Leitfähigkeit Statorwicklung quer zum Leiter
lam_SSL_quer = 2; % W*(m*K)^-1
% Spezial: Leitfähigkeit Wickelkopf in axialer Richtung
lam_WH_axial = 33; % W*(m*K)^-1
lam_WH_radial = 12; % W*(m*K)^-1

% Wärmeübergangskoeffizient zwischen Statorrücken und Gehäuse
h_contactHO = 1500; % W*(m^2*K)^-1; Wärmeübergang am Kontakt Statorrücken zu Gehäuse: nach Mellor 1991: 300...3000, je nach Flächenpressung in Kontaktfläche

% Spezifische Wärmekapazitäten zur Berechnung thermischer Kenngrößen
% --> Werkstoffe
cp_Cu = 385; % J*(kg*K)^-1; spez. Wärmekapazität Kupfer
cp_Al = 896; % J*(kg*K)^-1; spez. Wärmekapazität Aluminium
cp_Verguss = 1389; % J*(kg*K)^-1; spez. Wärmekapazität Vergussmasse
cp_Lack = cp_Verguss; % J*(kg*K)^-1; spez. Wärmekapazität Lackisolation um Leiter der Statorwicklung
% --> Komponenten
cp_ES = 441; % J*(kg*K)^-1; spez. Wärmekapazität Elektroblech; hier VACOFLUX50
cp_HO = 490; % J*(kg*K)^-1; spez. Wärmekapazität Gehäuse; hier Stahl
cp_EC = 490; % J*(kg*K)^-1; spez. Wärmekapazität Gehäuse-Stirndeckel; hier Stahl
cp_SH = 490; % J*(kg*K)^-1; spez. Wärmekapazität Rotorwelle; hier Stahl

%% Bereich 2:
%% Spezifische Wärmekapazitäten und -leitfähigkeiten der Komponenten
% Zuordnung der Leitermaterialen in Rotor & Stator
if strncmp(Entwurf.Optionen.Rotor_Leitermaterial.Bezeichnung,'Kupfer',3) == 1
    cp_RSL = cp_Cu; % spez. Wärmekapazität Leiter in Rotornuten
    cp_ER = cp_Cu; % spez. Wärmekapazität Leiter im Kurzschlussring
    lam_RSL = lam_Cu; % spez. Wärmeleitfähigkeit Leiter in Rotornuten
    lam_ER = lam_Cu; % spez. Wärmeleitfähigkeit Leiter im Kurzschlussring
elseif strncmp(Entwurf.Optionen.Rotor_Leitermaterial.Bezeichnung,'Alu',3) == 1
    cp_RSL = cp_Al; % spez. Wärmekapazität Leiter in Rotornuten
    cp_ER = cp_Al; % spez. Wärmekapazität Leiter im Kurzschlussring
    lam_RSL = lam_Al; % spez. Wärmeleitfähigkeit Leiter in Rotornuten
    lam_ER = lam_Al; % spez. Wärmeleitfähigkeit Leiter im Kurzschlussring
else
    disp('Das Rotor-Leitermaterial ist nicht bekannt. Ermittlung der geometrischen und thermischen Maschinenparameter nicht möglich.')
    keyboard;
end

if strncmp(Entwurf.Optionen.Stator_Leitermaterial.Bezeichnung,'Kupfer',3) == 1
    cp_1L = cp_Cu; % spez. Wärmekapazität Wicklung in Statornuten
    lam_SSL_axial = lam_Cu; % spez. Wärmeleitfähigkeit Wicklung in Statornuten
elseif strncmp(Entwurf.Optionen.Stator_Leitermaterial.Bezeichnung,'Alu',3) == 1
    cp_1L = cp_Al; % spez. Wärmekapazität Wicklung in Statornuten
    lam_SSL_axial = lam_Al; % spez. Wärmeleitfähigkeit Wicklung in Statornuten
else
    disp('Das Stator-Leitermaterial ist nicht bekannt. Ermittlung der geometrischen und thermischen Maschinenparameter nicht möglich.')
    keyboard;
end

% Komponenten aus Elektroblech
cp_SYO = cp_ES;
cp_STO = cp_ES;
cp_RTO = cp_ES;
cp_RYO = cp_ES;
lam_SYO = lam_ES;
lam_STO = lam_ES;
lam_RTO = lam_ES;
lam_RYO = lam_ES;

%% Berechnung der Massen der einzelnen Komponenten
m_HO_a = 1/n_dis*rho_HO*Entwurf.Geometrie.l*pi/4*((Entwurf.Geometrie.D_1a+2*b_HO)^2-Entwurf.Geometrie.D_1a^2); % kg; Masse Gehäuse auf aktiver Länge (ohne Stirnseitige Deckel!)
m_HO_p = 1/2*rho_HO*(LF-1)*Entwurf.Geometrie.l*pi/4*((Entwurf.Geometrie.D_1a+2*b_HO)^2-Entwurf.Geometrie.D_1a^2); % kg; Masse Gehäuse passiv (Überhang links & rechts über aktive Länge) (ohne Stirnseitige Deckel!)
m_EC = rho_HO*b_EC*pi/4*((Entwurf.Geometrie.D_1a+2*b_HO)^2-(Entwurf.Geometrie.D_2i+2*b_BR)^2); % kg; Masse Stirndeckel
m_SYO = 1/n_dis*Entwurf.Optionen.Stator_Eisenmaterial.rho_Fe*Entwurf.Geometrie.l*pi/4*(Entwurf.Geometrie.D_1a^2-(Entwurf.Geometrie.D_1i+2*Entwurf.Geometrie.Nut_1.h_1n*1e-3)^2); % kg; Masse Statorjoch
m_STO = 1/n_dis*Entwurf.Optionen.Stator_Eisenmaterial.rho_Fe*Entwurf.Geometrie.l*Entwurf.Geometrie.A_1z; % kg; Masse Statorzähne
    m_SSL_L = 1/n_dis*Entwurf.Geometrie.l*Entwurf.Wicklung.N_1*Entwurf.Optionen.Stator_Leitermaterial.rho_Le*Entwurf.Wicklung.A_1L*1e-6*Entwurf.Wicklung.z_1n; % kg; Masse des Leiters in den Statornuten
    m_SSL_V = 1/n_dis*Entwurf.Geometrie.l*Entwurf.Wicklung.N_1*(1-Entwurf.Richtwerte.phi_1n)*Entwurf.Geometrie.Nut_1.A_1n*1e-6*rho_Verguss; % kg; Masse Vergussmaterial + Leiterisolierung + Nutpapierisolierung
m_SSL = m_SSL_L + m_SSL_V; % kg; Masse in Statornuten gesamt
    m_WH1_L = (1+WHasym)*1/2*Entwurf.Optionen.Stator_Leitermaterial.rho_Le*Entwurf.Wicklung.l_1w*Entwurf.Wicklung.A_1L*1e-6*Entwurf.Wicklung.z_1n*Entwurf.Wicklung.N_1; % kg; Masse Leiter in Wicklungskopf
    m_WH2_L = (1-WHasym)*1/2*Entwurf.Optionen.Stator_Leitermaterial.rho_Le*Entwurf.Wicklung.l_1w*Entwurf.Wicklung.A_1L*1e-6*Entwurf.Wicklung.z_1n*Entwurf.Wicklung.N_1; % kg; Masse Leiter in Wicklungskopf
    A_Lack = pi/4*((Entwurf.Wicklung.d_1L*1e-3+2*d_Lack*1e-3)^2-(Entwurf.Wicklung.d_1L*1e-3)^2); % m^2; Hilfsgröße: Querschnittsfläche des Lackes im Axialschnitt des Leiters
    m_WH1_Lack = (1+WHasym)*1/2*(Entwurf.Wicklung.l_1w*Entwurf.Wicklung.z_1n*Entwurf.Wicklung.N_1)*A_Lack*rho_Lack; % kg; Masse des Isolierlackes im Wickelkopf
    m_WH2_Lack = (1-WHasym)*1/2*(Entwurf.Wicklung.l_1w*Entwurf.Wicklung.z_1n*Entwurf.Wicklung.N_1)*A_Lack*rho_Lack; % kg; Masse des Isolierlackes im Wickelkopf
m_WH1 = m_WH1_L + m_WH1_Lack; % kg; Gesamtmasse des Wickelkopfes
m_WH2 = m_WH2_L + m_WH2_Lack; % kg; Gesamtmasse des Wickelkopfes
m_RTO = 1/n_dis*Entwurf.Optionen.Rotor_Eisenmaterial.rho_Fe*Entwurf.Geometrie.l*(Entwurf.Geometrie.A_2r-Entwurf.Wicklung.N_2*Entwurf.Geometrie.Nut_2.A_2n_tat*1e-6-pi/4*((Entwurf.Geometrie.D_2i+2*Entwurf.Geometrie.Nut_2.h_2r*1e-3)^2-Entwurf.Geometrie.D_2i^2)); % kg; Masse Rotorzähne
m_RSL = 1/n_dis*Entwurf.Optionen.Rotor_Leitermaterial.rho_Le*Entwurf.Geometrie.l*Entwurf.Wicklung.N_2*Entwurf.Wicklung.A_2s*1e-6; % kg; Masse Rotorstäbe in Rotornuten
m_ER = Entwurf.Optionen.Rotor_Leitermaterial.rho_Le*pi*Entwurf.Geometrie.D_2r*Entwurf.Wicklung.A_2r*1e-6; % kg; Masse Kurzschlussring Motor
m_RYO = 1/n_dis*Entwurf.Optionen.Rotor_Eisenmaterial.rho_Fe*Entwurf.Geometrie.l*pi/4*((Entwurf.Geometrie.D_2i+2*Entwurf.Geometrie.Nut_2.h_2r*1e-3)^2-Entwurf.Geometrie.D_2i^2); % kg; Masse Rotorrücken (Rotorjoch)
m_SH_a = 1/n_dis*rho_SH*Entwurf.Geometrie.l*pi/4*(Entwurf.Geometrie.D_2i^2-(Entwurf.Geometrie.D_2i-2*b_SH)^2); % kg; Masse Rotorwelle auf aktiver Länge
m_SH_p = 1/2*rho_SH*((LF-1)*Entwurf.Geometrie.l+2*b_EC)*pi/4*(Entwurf.Geometrie.D_2i^2-(Entwurf.Geometrie.D_2i-2*b_SH)^2); % kg; Masse Rotorwelle passiv (Überhang links & rechts über aktive Länge)

%% Berechnung spezifische Wärmekapazität
cp_SSL = (cp_1L*m_SSL_L+cp_Verguss*m_SSL_V)/m_SSL; % J*(kg*K)^-1
cp_WH1 = (cp_1L*m_WH1_L+cp_Lack*m_WH1_Lack)/m_WH1; % J*(kg*K)^-1
cp_WH2 = (cp_1L*m_WH2_L+cp_Lack*m_WH2_Lack)/m_WH2; % J*(kg*K)^-1

%% Berechnung Gesamtmasse und spezifische Wärmekapazität der gesamten Maschine
m_EM_ges = n_dis*m_HO_a+2*m_HO_p+2*m_EC+n_dis*m_SYO+n_dis*m_STO+n_dis*m_SSL+m_WH1+m_WH2+n_dis*m_RTO+n_dis*m_RSL+2*m_ER+n_dis*m_RYO+n_dis*m_SH_a+2*m_SH_p; % kg; Gesamtmasse E-Maschine
cp_EM_ges = ((n_dis*m_HO_a+2*m_HO_p)*cp_HO+2*m_EC*cp_EC+n_dis*m_SYO*cp_SYO+n_dis*m_STO*cp_STO+n_dis*m_SSL*cp_SSL+m_WH1*cp_WH1+m_WH2*cp_WH2+n_dis*m_RTO*cp_RTO+n_dis*m_RSL*cp_RSL+2*m_ER*cp_ER+n_dis*m_RYO*cp_RYO+(n_dis*m_SH_a+2*m_SH_p)*cp_SH)/m_EM_ges; % J*(kg*K)^-1; mittlere spez. Wärmekapazität E-Maschine
% Ausgabe Gesamtmasse als Information an den Nutzer
fprintf('\nDie Gesamtmasse der E-Maschine beträgt %.2f kg.\n\n',m_EM_ges);

%% Berechnung der Wärmeleitwiderstände zwischen den einzelnen Komponenten
R_toWH = Entwurf.Geometrie.l/(2*n_dis)/lam_SSL_axial/(Entwurf.Wicklung.N_1*Entwurf.Wicklung.z_1n*Entwurf.Wicklung.A_1L*1e-6); % K*W^-1; Wärmeleitwiderstand Statorwicklung axial zum Wickelkopf
R_toER = Entwurf.Geometrie.l/(2*n_dis)/lam_RSL/(Entwurf.Wicklung.N_2*Entwurf.Wicklung.A_2s*1e-6); % K*W^-1; Wärmeleitwiderstand Rotorstäbe axial zum Kurzschlussring
R_toSTO = 1/(2*Entwurf.Wicklung.N_1)*(0.5*Entwurf.Geometrie.Nut_1.b_1n_m*1e-3/lam_SSL_quer/(Entwurf.Geometrie.l/n_dis*Entwurf.Geometrie.Nut_1.h_1l*1e-3)+Entwurf.Geometrie.Nut_1.d_1iso*1e-3/lam_iso/(Entwurf.Geometrie.l/n_dis*Entwurf.Geometrie.Nut_1.h_1l*1e-3)); % K*W^-1; Wärmeleitwiderstand Statorwicklung tangential (Umfangsrichtung) zu Statorzähnen
R_toRTO = 1/(2*Entwurf.Wicklung.N_2)*(0.5*Entwurf.Geometrie.Nut_2.b_2n_m*1e-3/lam_RSL/(Entwurf.Geometrie.l/n_dis*Entwurf.Geometrie.Nut_2.h_2l*1e-3)+0*Entwurf.Geometrie.Nut_2.d_2iso*1e-3/lam_iso/(Entwurf.Geometrie.l/n_dis*Entwurf.Geometrie.Nut_2.h_2l*1e-3)); % K*W^-1; Wärmeleitwiderstand tangential (Umfangsrichtung) Rotorstab zum Rotorzahn
R_toSAGP = 1/Entwurf.Wicklung.N_1*(0.5*Entwurf.Geometrie.Nut_1.h_1l*1e-3/lam_SSL_quer/(Entwurf.Geometrie.l/n_dis*0.5*(Entwurf.Geometrie.Nut_1.b_1n_u+Entwurf.Geometrie.Nut_1.b_1n_m)*1e-3)); % K*W^-1; Wärmeleitwiderstand radial Statorwicklung zum Luftspalt
R_toRAGP = 1/Entwurf.Wicklung.N_2*(0.5*Entwurf.Geometrie.Nut_2.h_2l*1e-3/lam_RSL/(Entwurf.Geometrie.l/n_dis*0.5*(Entwurf.Geometrie.Nut_2.b_2n_u+Entwurf.Geometrie.Nut_2.b_2n_m)*1e-3)); % K*W^-1; Wärmeleitwiderstand radial Rotorstäbe zum Luftspalt
R_toSYO = 1/Entwurf.Wicklung.N_1*(0.5*Entwurf.Geometrie.Nut_1.h_1l*1e-3/lam_SSL_quer/(Entwurf.Geometrie.l/n_dis*0.5*(Entwurf.Geometrie.Nut_1.b_1n_o+Entwurf.Geometrie.Nut_1.b_1n_m)*1e-3)+Entwurf.Geometrie.Nut_1.d_1iso*1e-3/lam_iso/(Entwurf.Geometrie.l/n_dis*Entwurf.Geometrie.Nut_1.b_1n_o*1e-3)); % K*W^-1; Wärmeleitwiderstand radial Statorwicklung zum Statorjoch (Statorrücken)
R_toRYO = 1/Entwurf.Wicklung.N_2*(0.5*Entwurf.Geometrie.Nut_2.h_2l*1e-3/lam_RSL/(Entwurf.Geometrie.l/n_dis*0.5*(Entwurf.Geometrie.Nut_2.b_2n_o+Entwurf.Geometrie.Nut_2.b_2n_m)*1e-3)+0*Entwurf.Geometrie.Nut_2.d_2iso*1e-3/lam_iso/(Entwurf.Geometrie.l/n_dis*Entwurf.Geometrie.Nut_2.b_2n_o*1e-3)); % K*W^-1; Wärmeleitwiderstand radial Rotorstäbe zum Rotorrücken
R_toSSL2 = 1/(2*Entwurf.Wicklung.N_1)*0.5*Entwurf.Geometrie.Nut_1.b_1z_m*1e-3/lam_STO/(Entwurf.Geometrie.l/n_dis*Entwurf.Geometrie.Nut_1.h_1n*1e-3); % K*W^-1; Wärmeleitwiderstand Statorzahn tangential (Umfangsrichtung) zur Statorwicklung
R_toRSL2 = 1/(2*Entwurf.Wicklung.N_2)*0.5*Entwurf.Geometrie.Nut_2.b_2z_m*1e-3/lam_RTO/(Entwurf.Geometrie.l/n_dis*Entwurf.Geometrie.Nut_2.h_2n*1e-3); % K*W^-1; Wärmeleitwiderstand Rotorzahn tangential (Umfangsrichtung) zu Rotorstäben
R_toSYO2 = 1/Entwurf.Wicklung.N_1*0.5*Entwurf.Geometrie.Nut_1.h_1n*1e-3/lam_STO/(Entwurf.Geometrie.l/n_dis*0.5*(Entwurf.Geometrie.Nut_1.b_1z_o+Entwurf.Geometrie.Nut_1.b_1z_m)*1e-3); % K*W^-1; Wärmeleitwiderstand radial Statorzahn zum Statorjoch (Statorrücken)
R_toRYO2 = 1/Entwurf.Wicklung.N_2*0.5*Entwurf.Geometrie.Nut_2.h_2n*1e-3/lam_RTO/(Entwurf.Geometrie.l/n_dis*0.5*(Entwurf.Geometrie.Nut_2.b_2z_o+Entwurf.Geometrie.Nut_2.b_2z_m)*1e-3); % K*W^-1; Wärmeleitwiderstand radial Rotorzahn zum Rotorrücken
R_toSAGP2 = 1/Entwurf.Wicklung.N_1*0.5*Entwurf.Geometrie.Nut_1.h_1n*1e-3/lam_STO/(Entwurf.Geometrie.l/n_dis*0.5*(Entwurf.Geometrie.Nut_1.b_1z_u+Entwurf.Geometrie.Nut_1.b_1z_m)*1e-3); % K*W^-1; Wärmeleitwiderstand radial Statorzahn zum Luftspalt
R_toRAGP2 = 1/Entwurf.Wicklung.N_2*0.5*Entwurf.Geometrie.Nut_2.h_2n*1e-3/lam_RTO/(Entwurf.Geometrie.l/n_dis*0.5*(Entwurf.Geometrie.Nut_2.b_2z_u+Entwurf.Geometrie.Nut_2.b_2z_m)*1e-3); % K*W^-1; Wärmeleitwiderstand radial Rotorzahn zum Luftspalt
    D_1r_i = Entwurf.Geometrie.D_1i+2*Entwurf.Geometrie.Nut_1.h_1n*1e-3; % m; Hilfsgröße: Innendurchmesser Statorjoch (Statorrücken)
    D_1r_m = (1/2*(Entwurf.Geometrie.D_1a^2+D_1r_i^2))^0.5; % m; Hilfsgröße: Durchmesser zum radialen Massenschwerpunkt des Statorjochs
R_toSSLTO = 0.5*(D_1r_m-D_1r_i)/lam_SYO/(Entwurf.Geometrie.l/n_dis*pi*0.5*(D_1r_i+D_1r_m)); % K*W^-1; Wärmeleitwiderstand vom Statorjoch zu Statorzähnen und -wicklung
R_toHO = 0.5*(Entwurf.Geometrie.D_1a-D_1r_m)/lam_SYO/(Entwurf.Geometrie.l/n_dis*pi*0.5*(Entwurf.Geometrie.D_1a+D_1r_m)); % K*W^-1; Wärmeleitwiderstand vom Statorjoch radial zum Gehäuse
    D_2r_a = Entwurf.Geometrie.D_2i+2*Entwurf.Geometrie.Nut_2.h_2r*1e-3; % m; Hilfsgröße: Außendurchmesser Rotorrücken
    D_2r_m = (1/2*(Entwurf.Geometrie.D_2i^2+D_2r_a^2))^0.5; % m; Hilfsgröße: Durchmesser zum radialen Massenschwerpunkt des Rotorrückens
R_toRSLTO = 0.5*(D_2r_a-D_2r_m)/lam_RYO/(Entwurf.Geometrie.l/n_dis*pi*0.5*(D_2r_a+D_2r_m)); % K*W^-1; Wärmeleitwiderstand vom Rotorrücken radial zu Rotorzähnen und -stäben
R_toSH = 0.5*(D_2r_m-Entwurf.Geometrie.D_2i)/lam_RYO/(Entwurf.Geometrie.l/n_dis*pi*0.5*(D_2r_m+Entwurf.Geometrie.D_2i)); % K*W^-1; Wärmeleitwiderstand vom Rotorrücken radial zur Rotorwelle
    D_EC_m = (1/2*((Entwurf.Geometrie.D_1a+2*b_HO)^2+(Entwurf.Geometrie.D_2i+2*b_BR)^2))^0.5; % m; Hilfsgröße: Radius zum radialen Massenschwerpunkt des Gehäuse-Stirndeckels
R_toHO2 = (0.5*(Entwurf.Geometrie.D_1a+2*b_HO-D_EC_m))/lam_EC/(b_EC*pi*0.5*(Entwurf.Geometrie.D_1a+2*b_HO+D_EC_m)); % K*W^-1; Wärmeleitwiderstand Gehäuse-Stirndeckel radial nach außen zum Gehäuse
R_toBR2 = (D_EC_m-0.5*Entwurf.Geometrie.D_2i-b_BR)/lam_EC/(b_EC*pi*(0.5*Entwurf.Geometrie.D_2i+b_BR+D_EC_m)); % K*W^-1; Wärmeleitwiderstand Gehäuse-Stirndeckel radial nach innen zum Lager
R_toBR_a = (0.5*Entwurf.Geometrie.l/n_dis)/lam_SH/(pi/4*(Entwurf.Geometrie.D_2i^2-(Entwurf.Geometrie.D_2i-2*b_SH)^2)); % K*W^-1; Wärmeleitwiderstand Rotorwelle axial zum Lager auf aktiver Länge
R_toBR_p = 0.5*(0.5*(LF-1)*Entwurf.Geometrie.l+b_EC)/lam_SH/(pi/4*(Entwurf.Geometrie.D_2i^2-(Entwurf.Geometrie.D_2i-2*b_SH)^2)); % K*W^-1; Wärmeleitwiderstand Rotorwelle axial zum Lager auf passivem Überhang links und rechts über aktive Länge
    D_HO_m = (1/2*(Entwurf.Geometrie.D_1a^2+(Entwurf.Geometrie.D_1a+2*b_HO)^2))^0.5; % m; Hilfsgröße; Durchmesser zum radialen Massenschwerpunkt des Gehäuses
R_toSYO3 = 0.5*(D_HO_m-Entwurf.Geometrie.D_1a)/lam_HO/(Entwurf.Geometrie.l/n_dis*pi*0.5*(D_HO_m+Entwurf.Geometrie.D_1a)); % K*W^-1; Wärmeleitwiderstand des Gehäuses radial nach innen zum Statorjoch (auf aktiver Länge)
R_toWC_a = 0.5*(Entwurf.Geometrie.D_1a+2*b_HO-D_HO_m)/lam_HO/(Entwurf.Geometrie.l/n_dis*pi*0.5*(D_HO_m+Entwurf.Geometrie.D_1a+2*b_HO)); % K*W^-1; Wärmeleitwiderstand des Gehäuses radial nach außen zur Wasserkühlung (auf aktiver Länge)
R_toWC_p = 0.5*(Entwurf.Geometrie.D_1a+2*b_HO-D_HO_m)/lam_HO/(0.5*Entwurf.Geometrie.l*(LF-1)*pi*0.5*(D_HO_m+Entwurf.Geometrie.D_1a+2*b_HO)); % K*W^-1; Wärmeleitwiderstand des Gehäuses radial nach außen zur Wasserkühlung (auf passiver Länge)
R_toEC_a = (0.5*Entwurf.Geometrie.l/n_dis)/lam_HO/(pi/4*((Entwurf.Geometrie.D_1a+2*b_HO)^2-Entwurf.Geometrie.D_1a^2)); % K*W^-1; Wärmeleitwiderstand des Gehäuses axial zu den Gehäuse-Stirndeckeln auf aktiver Länge
R_toEC_p = (0.5*(LF-1)*Entwurf.Geometrie.l*0.5)/lam_HO/(pi/4*((Entwurf.Geometrie.D_1a+2*b_HO)^2-Entwurf.Geometrie.D_1a^2)); % K*W^-1; Wärmeleitwiderstand des Gehäuses axial zu den Gehäuse-Stirndeckeln auf passivem Überhang links und rechts über aktive Länge
    b_WH1 = (m_WH1_L/Entwurf.Optionen.Stator_Leitermaterial.rho_Le+m_WH1_Lack/rho_Lack)/(pi*(Entwurf.Geometrie.D_1i+Entwurf.Geometrie.Nut_1.h_1l*1e-3))/(Entwurf.Geometrie.Nut_1.h_1l*1e-3); % m; Hilfsgröße; Breite (axial) des Wickelkopfes
    b_WH2 = (m_WH2_L/Entwurf.Optionen.Stator_Leitermaterial.rho_Le+m_WH2_Lack/rho_Lack)/(pi*(Entwurf.Geometrie.D_1i+Entwurf.Geometrie.Nut_1.h_1l*1e-3))/(Entwurf.Geometrie.Nut_1.h_1l*1e-3); % m; Hilfsgröße; Breite (axial) des Wickelkopfes
R_toSSL_1 = 0.5*b_WH1/lam_WH_axial/(pi/4*((Entwurf.Geometrie.D_1i+2*Entwurf.Geometrie.Nut_1.h_1l*1e-3)^2-Entwurf.Geometrie.D_1i^2)); % K*W^-1; Wärmeleitwiderstand axial im Wickelkopf in beide Richtungen (Nut und Innenluft)
R_toSSL_2 = 0.5*b_WH2/lam_WH_axial/(pi/4*((Entwurf.Geometrie.D_1i+2*Entwurf.Geometrie.Nut_1.h_1l*1e-3)^2-Entwurf.Geometrie.D_1i^2)); % K*W^-1; Wärmeleitwiderstand axial im Wickelkopf in beide Richtungen (Nut und Innenluft)
    b_ER = (Entwurf.Wicklung.A_2r*1e-6)/(Entwurf.Geometrie.Nut_2.h_2l*1e-3); % m; Hilfsgröße; Breite (axiale Richtung) des Kurzschlussringes
R_toRSL = 0.5*b_ER/lam_ER/(pi/4*((Entwurf.Geometrie.D_2r+Entwurf.Geometrie.Nut_2.h_2l*1e-3)^2-(Entwurf.Geometrie.D_2r-Entwurf.Geometrie.Nut_2.h_2l*1e-3)^2)); % K*W^-1; Wärmeleitwiderstand axial im Kurzschlussring in beide Richtungen (Nut und Innenluft)
R_contactHO = 1/(h_contactHO*Entwurf.Geometrie.l/n_dis*pi*Entwurf.Geometrie.D_1a); % K*W^-1; Kontaktwiderstand zwischen Gehäuse und Statorrücken aufgrund der nicht-stoffschlüssigen Verbindung (Pressverband=kraftschlüssig)
R_toIA_1 = (0.5*Entwurf.Geometrie.Nut_1.h_1l*1e-3)/lam_WH_radial/(b_WH1*pi*(Entwurf.Geometrie.D_1i+1.5*Entwurf.Geometrie.Nut_1.h_1l*1e-3)); % K*W^-1; Wärmeleitwiderstand Wickelkopf radial zur Außenfläche
R_toIA_2 = (0.5*Entwurf.Geometrie.Nut_1.h_1l*1e-3)/lam_WH_radial/(b_WH2*pi*(Entwurf.Geometrie.D_1i+1.5*Entwurf.Geometrie.Nut_1.h_1l*1e-3)); % K*W^-1; Wärmeleitwiderstand Wickelkopf radial zur Außenfläche
R_toIA2_1 = (0.5*Entwurf.Geometrie.Nut_1.h_1l*1e-3)/lam_WH_radial/(b_WH1*pi*(Entwurf.Geometrie.D_1i+0.5*Entwurf.Geometrie.Nut_1.h_1l*1e-3)); % K*W^-1; Wärmeleitwiderstand Wickelkopf radial zur Innenfläche
R_toIA2_2 = (0.5*Entwurf.Geometrie.Nut_1.h_1l*1e-3)/lam_WH_radial/(b_WH2*pi*(Entwurf.Geometrie.D_1i+0.5*Entwurf.Geometrie.Nut_1.h_1l*1e-3)); % K*W^-1; Wärmeleitwiderstand Wickelkopf radial zur Innenfläche
R_toIA3 = 0.5*Entwurf.Geometrie.Nut_2.h_2l*1e-3/lam_ER/(b_ER*pi*(Entwurf.Geometrie.D_2r-0.5*Entwurf.Geometrie.Nut_2.h_2l*1e-3)); % K*W^-1; Wärmeleitwiderstand Kurzschlussring radial zur Innenfläche
R_toIA4 = 0.5*Entwurf.Geometrie.Nut_2.h_2l*1e-3/lam_ER/(b_ER*pi*(Entwurf.Geometrie.D_2r+0.5*Entwurf.Geometrie.Nut_2.h_2l*1e-3)); % K*W^-1; Wärmeleitwiderstand Kurzschlussring radial zur Außenfläche
R_toIA5 = (D_HO_m-0.5*Entwurf.Geometrie.D_1a)/lam_HO/(0.5*(LF-1)*Entwurf.Geometrie.l*pi*(D_HO_m+0.5*Entwurf.Geometrie.D_1a)); % K*W^-1; Wärmeleitwiderstand des Gehäuses radial nach innen zur Innenluft (auf passiver Länge)
    D_SH_m = (0.5*(Entwurf.Geometrie.D_2i^2+(Entwurf.Geometrie.D_2i-2*b_SH)^2))^0.5; % m; Hilfsgröße: Durchmesser zum Massenschwerpunkt der Rotorwelle im Radialschnitt
R_toRYO3 = n_dis/lam_SH*0.5*(Entwurf.Geometrie.D_2i-D_SH_m)/(Entwurf.Geometrie.l*pi*0.5*(Entwurf.Geometrie.D_2i+D_SH_m)); % K*W^-1; Wärmeleitwiderstand Rotorwelle radial nach außen zum Rotorrücken
R_toIA6 = 0.5*(Entwurf.Geometrie.D_2i-D_SH_m)/lam_SH/(0.5*(LF-1)*Entwurf.Geometrie.l*pi*0.5*(Entwurf.Geometrie.D_2i+D_SH_m)); % K*W^-1; Wärmeleitwiderstand Rotorwelle radial nach außen zur Innenluft (auf passiver Länge)
R_contactBR = (l_BR_eq*1e-3)/lam_Luft/(pi*(Entwurf.Geometrie.D_2i+b_BR)*b_EC); % K*W^-1; Kontaktwiderstand im Kugellager zwischen Innen- und Außenschale
R_toIA7 = 0.5*b_EC/lam_EC/(pi/4*(Entwurf.Geometrie.D_1a^2-(Entwurf.Geometrie.D_2i+2*b_BR)^2)); % K*W^-1; Wärmeleitwiderstand Gehäusedeckel axial in Richtung Innenluft

%% Berechnung der relevanten Größen für den konvektiven Wärmeübergang in der E-Maschine
% Gehäuse Außenfläche -- Konvektion zum Kühlmittel
A_HO_a = Entwurf.Geometrie.l/n_dis*pi*(Entwurf.Geometrie.D_1a+2*b_HO); % m^2
A_HO_p = (LF-1)*Entwurf.Geometrie.l/2*pi*(Entwurf.Geometrie.D_1a+2*b_HO); % m^2
% Konvektion im Luftspalt zwischen Stator und Rotor
A_SSL = pi*Entwurf.Geometrie.D_1i*Entwurf.Geometrie.l/n_dis*Entwurf.Geometrie.Nut_1.b_1n_u/(Entwurf.Geometrie.Nut_1.b_1n_u+Entwurf.Geometrie.Nut_1.b_1z_u); % m^2; Statorwicklung zum Luftspalt
A_STO = pi*Entwurf.Geometrie.D_1i*Entwurf.Geometrie.l/n_dis*Entwurf.Geometrie.Nut_1.b_1z_u/(Entwurf.Geometrie.Nut_1.b_1n_u+Entwurf.Geometrie.Nut_1.b_1z_u); % m^2; Statorzähne zum Luftspalt
A_RSL = pi*Entwurf.Geometrie.D_2a*Entwurf.Geometrie.l/n_dis*Entwurf.Geometrie.Nut_2.b_2n_u/(Entwurf.Geometrie.Nut_2.b_2n_u+Entwurf.Geometrie.Nut_2.b_2z_u); % m^2; Rotostäbe zum Luftspalt
A_RTO = pi*Entwurf.Geometrie.D_2a*Entwurf.Geometrie.l/n_dis*Entwurf.Geometrie.Nut_2.b_2z_u/(Entwurf.Geometrie.Nut_2.b_2n_u+Entwurf.Geometrie.Nut_2.b_2z_u); % m^2; Rotorzähne zum Luftspalt
A_AGP = n_dis*(A_SSL+A_STO+A_RSL+A_RTO); % m^2; Gesamtoberfläche Luftspalt: benötigt für Aufteilung der Lüfterverluste
% Konvektion in Innenluft
A_SH = 0.5*(LF-1)*Entwurf.Geometrie.l*pi*Entwurf.Geometrie.D_2i; % m^2; Oberfläche der Rotorwelle im Kontakt mit Innenluft (auf einer Seite)
    A_ER_axial = pi/4*((Entwurf.Geometrie.D_2r+Entwurf.Geometrie.Nut_2.h_2l*1e-3)^2-(Entwurf.Geometrie.D_2r-Entwurf.Geometrie.Nut_2.h_2l*1e-3)^2); % m^2; Stirnfläche Kurzschlussring zu Endcap-Luft
    A_ER_tang_a = b_ER*pi*(Entwurf.Geometrie.D_2r+Entwurf.Geometrie.Nut_2.h_2l*1e-3); % m^2; Kurzschlussring Mantelfläche (tangential) außen
    A_ER_tang_i = b_ER*pi*(Entwurf.Geometrie.D_2r-Entwurf.Geometrie.Nut_2.h_2l*1e-3); % m^2; Kurzschlussring Mantelfläche (tangential) innen
A_ER = A_ER_axial+A_ER_tang_a+A_ER_tang_i; % m^2; Oberfläche des Kurzschlussrings im Kontakt mit Innenluft(auf einer Seite)
    A_WH_axial = A_WH_Faktor*pi/4*((Entwurf.Geometrie.D_1i+2*Entwurf.Geometrie.Nut_1.h_1l*1e-3)^2-Entwurf.Geometrie.D_1i^2); % m^2; Wickelkopf Stirnfläche zur Innenluft
    A_WH1_tang_a = A_WH_Faktor*b_WH1*pi*(Entwurf.Geometrie.D_1i+2*Entwurf.Geometrie.Nut_1.h_1l*1e-3); % m^2; Wickelkopf Mantelfläche (tangential) außen
    A_WH1_tang_i = A_WH_Faktor*b_WH1*pi*Entwurf.Geometrie.D_1i; % m^2; Wickelkopf Mantelfläche (tangential) innen
    A_WH2_tang_a = A_WH_Faktor*b_WH2*pi*(Entwurf.Geometrie.D_1i+2*Entwurf.Geometrie.Nut_1.h_1l*1e-3); % m^2; Wickelkopf Mantelfläche (tangential) außen
    A_WH2_tang_i = A_WH_Faktor*b_WH2*pi*Entwurf.Geometrie.D_1i; % m^2; Wickelkopf Mantelfläche (tangential) innen
A_WH1 = A_WH_axial+A_WH1_tang_a+A_WH1_tang_i; % m^2; Oberfläche des Wickelkopfes im Kontakt mit Innenluft (auf einer Seite)
A_WH2 = A_WH_axial+A_WH2_tang_a+A_WH2_tang_i; % m^2; Oberfläche des Wickelkopfes im Kontakt mit Innenluft (auf einer Seite)
A_EC = pi/4*(Entwurf.Geometrie.D_1a^2-(Entwurf.Geometrie.D_2i+2*b_BR)^2); % m^2; Fläche Gehäusedeckel-Innenseite zu Innenluft
A_HO2 = 0.5*(LF-1)*Entwurf.Geometrie.l*pi*Entwurf.Geometrie.D_1a; % m^2; Oberfläche Gehäuse im Kontakt mit Innenluft (auf einer Seite)
A_IA = A_WH1+A_WH2+2*(A_ER+A_SH+A_EC+A_HO2); % m^2; Gesamtoberfläche Innenluft: benötigt für Aufteilung der Lüfterverluste

% geometr. Größen, die zur Online-Berechnung der Innenluftkonvektion benötigt werden
D_2r = Entwurf.Geometrie.D_2r; % m; Mittlerer Durchmesser Kurzschlussring
% geometr. Größen, die zur Online-Berechnung der Luftspaltkonvektion benötigt werden
D_2a = Entwurf.Geometrie.D_2a; % m; Außendurchmesser des Rotors
delta = Entwurf.Geometrie.delta*1e-3; % m; Geometrische Luftspaltlänge

%% Berechnung der Kühlkanal-Geometrie
A_WC = LF*Entwurf.Geometrie.l*b_WC; % m^2; Querschnitt des Kühlkanals
l_WC = pi*(Entwurf.Geometrie.D_1a+2*b_HO+b_WC); % m; Länge des Kühlkanals
p_WC = 2*LF*Entwurf.Geometrie.l+2*b_WC; % m; Umfang des Kühlkanals

%% Aufteilung der Wicklungsverluste auf Rotor und Stator
P_Verlust_Rotor_Ringe = 2*Entwurf.Wicklung.N_2*Entwurf.EMAG.R_2r*Entwurf.EMAG.I_2r^2; % W
P_Verlust_Rotor_Staebe = Entwurf.Wicklung.N_2*Entwurf.EMAG.R_2s*Entwurf.EMAG.I_2s^2; % W
P_Verlust_Stator_Wicklung = Entwurf.Bemessungswerte.m*Entwurf.EMAG.R_1*Entwurf.EMAG.I_1Str^2; % W
Anteil_Stator = P_Verlust_Stator_Wicklung/(P_Verlust_Stator_Wicklung+P_Verlust_Rotor_Ringe+P_Verlust_Rotor_Staebe); % zur Aufteilung der Stromwärmeverluste auf Stator und Rotor
Anteil_Rotor_Ringe = P_Verlust_Rotor_Ringe/(P_Verlust_Rotor_Ringe+P_Verlust_Rotor_Staebe); % zur Aufteilung der Stromwärmeverluste im Rotor auf Ringe und Stäbe

%% Überprüfung Plausibilität der geometrischen Abmaße
% if b_WH > ((LF-1)*Entwurf.Geometrie.l/2) || b_ER > ((LF-1)*Entwurf.Geometrie.l/2)
%     fprintf('Kollision von Wickelkopf und/oder Kurzschlussring mit dem Gehäusedeckel.\nDer Längenfaktor LF sollte vergrößert werden.\n\n')
% end

clearvars lam_Al lam_Cu lam_E* lam_HO lam_iso lam_R* lam_S* lam_W*
clearvars cp_Al cp_Cu cp_ES cp_1L cp_Lack cp_Verguss
clearvars rho_HO rho_EC rho_SH rho_Verguss rho_Lack
clearvars b_* LF Entwurf D_1J* D_2J* P_Verlust* D_SH_m D_2r_* D_EC_m D_HO_m D_1r_*
clearvars A_Lack m_SSL_V m_WH_Lack
clearvars A_WH_Faktor