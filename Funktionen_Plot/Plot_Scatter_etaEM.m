% Skript zur Datstellung der aufgetretenen Wirkungsgrade in der Simulation
% im E-Maschinenkennfeld

% Berechnung motorischer und generatorischer Wirkungsgrad an den
% aufgetretenden Betriebspunkten der Simulation
eta_mot=(1./(1+EMaschine_Thermik_Output.P_EMaschine_Verlust_ist__W_.Data./Energiebilanz_Output.P_EMaschine_mech_mot__W_.Data));
eta_gen=(1+EMaschine_Thermik_Output.P_EMaschine_Verlust_ist__W_.Data./Energiebilanz_Output.P_EMaschine_mech_reku__W_.Data);

% Grafische Darstellung
figure
scatter3(Getriebe_Output.n_EMaschine__min__1_.Data,EMaschine_Verlustleistung_Output.M_EMaschine_ist__Nm_.Data,Getriebe_Output.n_EMaschine__min__1_.Time,10,(eta_mot).*100);
ylim([0 90])
colorbar
grid off
box on
xlabel('Drehzahl E-Maschine  n_{EM}  / min^{-1}')
ylabel('Drehmoment E-Maschine  M_{EM}  / Nm')
h=colorbar;
ylabel(h,'Wirkungsgrad E-Maschine  / %')
legend('\eta_{EM,mot}')
title('Wirkungsgrade der E-Maschine im Zyklus (motorisch)')

figure
scatter3(Getriebe_Output.n_EMaschine__min__1_.Data,EMaschine_Verlustleistung_Output.M_EMaschine_ist__Nm_.Data,Getriebe_Output.n_EMaschine__min__1_.Time,10,(eta_gen).*100);
ylim([-55 0])
colorbar
grid off
box on
xlabel('Drehzahl E-Maschine  n_{EM}  / min^{-1}')
ylabel('Drehmoment E-Maschine  M_{EM}  / Nm')
h=colorbar;
ylabel(h,'Wirkungsgrad E-Maschine  / %')
legend('\eta_{EM,gen}')
title('Wirkungsgrade der E-Maschine im Zyklus (generatorisch)')
