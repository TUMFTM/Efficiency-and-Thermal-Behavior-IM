% Skript, welches die Verlustteilung der mechanischen Verluste berechnet
% Ausgangsgröße: chi_Lager (relativer Verlustanteil des Lagers am gesamter
% mechanischer Verlustleistung, verwendet in 1D lookup table)
load('Entwurf_ASM_20190704_164434.mat') % Entwurfsparameter der E-Maschine aus Svenjas Berechnungstool

%% Bereich 1: ==> Hier können Werte abgeändert werden!
% Designparameter:
mue_Lager = 0.005; % Reibungskoeffizient der Lagerreibung; 0.005 für NEmo
k_AGP = 2.5; % Rauhigkeitskoeffizient der Oberflächen im Luftspalt; 2.5 für "slottet rotor" (nach NERG)

%% Bereich 2: Hier keine Werte ändern! (Berechnung)
% Berechnung Lagerverlustleistung
% Ansatz nach SKF einfach bzw. PYRHONEN
n_Lager=linspace(0,(v_max/3.6/r_dynamisch*i_Getriebe/2/pi*60)+10); % 1/min
omega_Lager=2*pi/60.*n_Lager; % 1/s
F_Lager = 0.5*9.81*(2*(m_ER+m_SH_p)+n_dis*(m_SH_a+m_RSL+m_RTO+m_RYO));
P_V_Lager = 2.*(0.5*mue_Lager*F_Lager*Entwurf.Geometrie.D_2i).*omega_Lager; % W; Verlust BEIDER Lager

% Berechnung Lüftungsverlust im Luftspalt
% Ansatz nach PYRHONEN
n_AGP=linspace(0,6100);
omega_AGP=2*pi/60.*n_AGP;
Re_delta_AGP = (Entwurf.Geometrie.D_2a*Entwurf.Geometrie.delta*1e-3/(2*nue_Luft)).*omega_AGP;
C_M_AGP=zeros(1,length(Re_delta_AGP));
for c=1:length(Re_delta_AGP)
    if Re_delta_AGP(c) ~= 0
        if Re_delta_AGP(c) < 64
            C_M_AGP(c) = 10*(2*Entwurf.Geometrie.delta*1e-3/Entwurf.Geometrie.D_2a)^0.3/Re_delta_AGP(c);
        elseif Re_delta_AGP(c) >= 64 && Re_delta_AGP(c) < 500
            C_M_AGP(c) = 2*(2*Entwurf.Geometrie.delta*1e-3/Entwurf.Geometrie.D_2a)^0.3/Re_delta_AGP(c)^0.6;
        elseif Re_delta_AGP(c) >= 500 && Re_delta_AGP(c) < 1e4
            C_M_AGP(c) = 1.03*(2*Entwurf.Geometrie.delta*1e-3/Entwurf.Geometrie.D_2a)^0.3/Re_delta_AGP(c)^0.5;
        else
            C_M_AGP(c) = 0.065*(2*Entwurf.Geometrie.delta*1e-3/Entwurf.Geometrie.D_2a)^0.3/Re_delta_AGP(c)^0.2;
        end
    end
end
P_V_AGP = (1/32*k_AGP*pi*rho_Luft_EM*Entwurf.Geometrie.D_2a^4*Entwurf.Geometrie.l).*C_M_AGP.*omega_AGP.^3;

% Berechnung Lüftungsverlust in der Innenluft
% dazu zunächst Bestimmung h_fconv (HTC der erzwungenen Konvektion) im AGP
v_AGP=2*pi/60*Entwurf.Geometrie.D_2a/2.*n_AGP;
Pr_AGP = cp_Luft*nue_Luft*rho_Luft_EM/lam_Luft;
Re_AGP = Entwurf.Geometrie.delta*1e-3/nue_Luft.*v_AGP;
Ta_AGP = ((Entwurf.Geometrie.delta*1e-3/(Entwurf.Geometrie.D_2a/2))^0.5).*Re_AGP;
Nu_AGP = zeros(1,length(Ta_AGP));
for c=1:length(Ta_AGP)
    if Ta_AGP(c) <= 100
        Nu_AGP(c) = Pr_AGP^0.27*0.202.*Ta_AGP(c)^0.63;
    else
        Nu_AGP(c) = Pr_AGP^0.27*0.386.*Ta_AGP(c)^0.5;
    end    
end
h_fconv_AGP = lam_Luft/(Entwurf.Geometrie.delta*1e-3).*Nu_AGP;
% dazu Bestimmung h_fconv in der Innenluft
v_IA = eta_IA*Entwurf.Geometrie.D_2r/2*2*pi/60.*n_AGP;
h_fconv_IA = 15.5*0.29.*v_IA;
% dazu Berechnung Anteil Innenluft am Gesamt-Lüftungsverlust
chi_Lueft_IA = h_fconv_IA.*A_IA./(h_fconv_AGP.*A_AGP+h_fconv_IA.*A_IA);
% Berechnung Innenluft-Lüftungsverlust
P_V_IA = P_V_AGP./(1./chi_Lueft_IA-1);

% Berechnung Endergebnis des Skriptes: Anteil Lager am gesamten mech. Verlust
chi_Lager = P_V_Lager./(P_V_Lager+P_V_AGP+P_V_IA);
chi_Lager(1) = 1;
chi_Lager = [n_Lager;chi_Lager];

clearvars *_Lager Entwurf *_AGP *_IA -except chi_Lager A_AGP A_IA eta_IA