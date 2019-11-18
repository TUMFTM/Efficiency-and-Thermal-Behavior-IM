% für "fast restart" Simulation mit anderem Output File

%--------------------------------------------------------------------------
% Darstellung der Umgebungstemperatur
%--------------------------------------------------------------------------

% figure
% plot([0;t_Simulation],[T_Umgebung;T_Umgebung]-273.15);
% grid on;
% axis tight;
% title('Umgebungstemperatur');
% xlabel('Zeit in s');
% ylabel('Temperatur in °C');


%--------------------------------------------------------------------------
% Darstellung der Fahrzeuggeschwindigkeit
%--------------------------------------------------------------------------

figure
plot([v_Fahrzeug.Time;t_Simulation],[v_Fahrzeug.Data(:);0]);
grid on;
axis tight;
hold on
plot(Rad_Output.v_Fahrzeug_ist__m_s__1_.Time,(Rad_Output.v_Fahrzeug_ist__m_s__1_.Data(:).*3.6))
title('Fahrzeuggeschwindigkeit');
xlabel('Zeit in s');
ylabel('Geschwindigkeit in km/h');
legend('Sollgeschw.','Istgeschw.')

%--------------------------------------------------------------------------
% Darstellung der Steigung der Strecke
%--------------------------------------------------------------------------

% figure
% plot([alpha_Strecke(:,1);t_Simulation],[alpha_Strecke(:,2);0]);
% grid on;
% axis tight;
% title('Steigung der Strecke');
% xlabel('Zeit in s');
% ylabel('Steigung in %');

%--------------------------------------------------------------------------
% Darstellung der Leistungsanforderungen an die Raeder
%--------------------------------------------------------------------------
% 
% if exist('Rad_Output')==1
%    figure
%    plot(Rad_Output.P_Rad__W_.Time,Rad_Output.P_Rad__W_.Data(:)/1000);
%    grid on;
%    axis tight;
%    title('Leistungsanforderung an die Raeder');
%    xlabel('Zeit in s');
%    ylabel('Leistung in kW');
% end
% 

%--------------------------------------------------------------------------
% Darstellung der Leistungsanforderungen an die 1. E-Maschine
%--------------------------------------------------------------------------

figure
plot(EMaschine_Verlustleistung_Output.M_EMaschine_ist__Nm_.Time,(2*pi/60/1000.*Getriebe_Output.n_EMaschine__min__1_.Data.*EMaschine_Verlustleistung_Output.M_EMaschine_ist__Nm_.Data));
%Vergleich mit gemessener EM Leistung
if exist('pmech_meas')==1
    hold on
    plot(pmech_meas)
    legend('Simulation','Messung')
end
%
grid on;
axis tight;
title('Leistungsanforderung an die 1. E-Maschine');
xlabel('Zeit in s');
ylabel('Leistung in kW');

%--------------------------------------------------------------------------
% Darstellung der fehlenden Leistung bei Überschreiten des maximalen Moments
%--------------------------------------------------------------------------

figure
plot(Getriebe_Output.n_EMaschine__min__1_.Time,(2*pi/60/1000.*Getriebe_Output.n_EMaschine__min__1_.Data(:).*EMaschine_Verlustleistung_Output.Regelabweichung__Nm_.Data(:)));
grid on;
axis tight;
title('Leistungsüberschreitung');
xlabel('Zeit in s');
ylabel('Mechanische Leistung E-Maschine kW');

%--------------------------------------------------------------------------
% Darstellung Drehzahlverlauf der E-Maschine
%--------------------------------------------------------------------------

%    figure
%    plot(Getriebe_Output.n_EMaschine__min__1_.Time,Getriebe_Output.n_EMaschine__min__1_.Data(:));
%    grid on;
%    axis tight;
%    title('Drehzahl der 1. E-Maschine');
%    xlabel('Zeit in s');
%    ylabel('Drehzahl in 1/min');

%--------------------------------------------------------------------------
% Darstellung Drehzahlverlauf der E-Maschine
%--------------------------------------------------------------------------

%    figure
%    plot(Getriebe_Output.M_EMaschine__Nm_.Time,Getriebe_Output.M_EMaschine__Nm_.Data(:));
%    grid on;
%    axis tight;
%    title('Drehmoment der 1. E-Maschine');
%    xlabel('Zeit in s');
%    ylabel('Moment in Nm');

%--------------------------------------------------------------------------
% Darstellung Verlustleistung der E-Maschine vor & nach Korrektur durch temperaturabhängige Verlustleistung
%--------------------------------------------------------------------------
 
figure
plot(EMaschine_Verlustleistung_Output.P_EMaschine_Verlust_korr__W_.Time,EMaschine_Verlustleistung_Output.P_EMaschine_Verlust_korr__W_.Data(:)/1000);
hold on
plot(EMaschine_Thermik_Output.P_EMaschine_Verlust_ist__W_.Time,EMaschine_Thermik_Output.P_EMaschine_Verlust_ist__W_.Data./1000);
grid on;
axis tight;
title('Verlustleistung E-Maschine ohne & mit temperaturabh. Verlustleistung');
legend('P_{V,EM,ohne}','P_{V,EM,mit}')
xlabel('Zeit in s');
ylabel('Leistung in kW');

%--------------------------------------------------------------------------
% Darstellung Wickelkopftemperatur der E-Maschine
%--------------------------------------------------------------------------

figure
plot(EMaschine_Thermik_Output.T_EMaschine_WH__K_.T_EMaschine_WH_1__K_.Time, EMaschine_Thermik_Output.T_EMaschine_WH__K_.T_EMaschine_WH_1__K_.Data-273.15)
% Vergleich mit Messdaten, falls Messfahrt simuliert wurde
if exist('T_Verlauf_EMaschine')==1
    hold on
    plot(T_Verlauf_EMaschine(:,1),T_Verlauf_EMaschine(:,2))
    legend('Simulation','Messung')
end
grid on
axis tight
title('Wickelkopftemperatur Stator der E-Maschine')
xlabel('Zeit in s');
ylabel('Temperatur in °C');

%--------------------------------------------------------------------------
% Darstellung Kühlmitteltemperatur der E-Maschine
%--------------------------------------------------------------------------

figure
plot(EMaschine_Thermik_Output.T_EMaschine_Kuehlmittel__K_.Time, EMaschine_Thermik_Output.T_EMaschine_Kuehlmittel__K_.Data-273.15)
% Vergleich mit Messdaten, falls Messfahrt simuliert wurde
if exist('T_Verlauf_Coolant')==1
    hold on
    plot(T_Verlauf_Coolant(:,1),T_Verlauf_Coolant(:,2)+8)
    legend('Simulation','Messung')
end
grid on
axis tight
title('Kühlmitteltemperatur')
xlabel('Zeit in s');
ylabel('Temperatur in °C');
