%% Hier Titel der Grafik angeben!
Titel = '6000 min^{-1}, ohne d_2iso, Artemis Urban, T_{max}';

%% Abbildung des thermischen Netzwerkes
strBR1 = sprintf('  BR1\n  ');
strEC1 = sprintf('  EC1\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_EC__K_.T_EMaschine_EC_1__K_.Data)-273.15);
strSHp1 = sprintf('  SHp1\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_SH__K_.T_EMaschine_SH_p_1__K_.Data)-273.15);
strIA1 = sprintf('  IA1\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_IA__K_.T_EMaschine_IA_1__K_.Data)-273.15);
strHOp1 = sprintf('  HOp1\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_HO__K_.T_EMaschine_HO_p_1__K_.Data)-273.15);
strER1 = sprintf('  ER1\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_ER__K_.T_EMaschine_ER_1__K_.Data)-273.15);
strWH1 = sprintf('  WH1\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_WH__K_.T_EMaschine_WH_1__K_.Data)-273.15);
strSHa1 = sprintf('  SHa1\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_SH__K_.T_EMaschine_SH_a_1__K_.Data)-273.15);
strRYO1 = sprintf('  RYO1\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_RYO__K_.T_EMaschine_RYO_1__K_.Data)-273.15);
strRTO1 = sprintf('  RTO1\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_RTO__K_.T_EMaschine_RTO_1__K_.Data)-273.15);
strAGP1 = sprintf('  AGP1\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_AGP__K_.T_EMaschine_AGP_1__K_.Data)-273.15);
strSTO1 = sprintf('  STO1\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_STO__K_.T_EMaschine_STO_1__K_.Data)-273.15);
strSYO1 = sprintf('  SYO1\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_SYO__K_.T_EMaschine_SYO_1__K_.Data)-273.15);
strHOa1 = sprintf('  HOa1\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_HO__K_.T_EMaschine_HO_a_1__K_.Data)-273.15);
strRSL1 = sprintf('    RSL1\n    %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_RSL__K_.T_EMaschine_RSL_1__K_.Data)-273.15);
strSSL1 = sprintf('    SSL1\n    %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_SSL__K_.T_EMaschine_SSL_1__K_.Data)-273.15);
strSHa2 = sprintf('  SHa2\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_SH__K_.T_EMaschine_SH_a_2__K_.Data)-273.15);
strRYO2 = sprintf('  RYO2\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_RYO__K_.T_EMaschine_RYO_2__K_.Data)-273.15);
strRTO2 = sprintf('  RTO2\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_RTO__K_.T_EMaschine_RTO_2__K_.Data)-273.15);
strAGP2 = sprintf('  AGP2\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_AGP__K_.T_EMaschine_AGP_2__K_.Data)-273.15);
strSTO2 = sprintf('  STO2\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_STO__K_.T_EMaschine_STO_2__K_.Data)-273.15);
strSYO2 = sprintf('  SYO2\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_SYO__K_.T_EMaschine_SYO_2__K_.Data)-273.15);
strHOa2 = sprintf('  HOa2\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_HO__K_.T_EMaschine_HO_a_2__K_.Data)-273.15);
strRSL2 = sprintf('    RSL2\n    %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_RSL__K_.T_EMaschine_RSL_2__K_.Data)-273.15);
strSSL2 = sprintf('    SSL2\n    %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_SSL__K_.T_EMaschine_SSL_2__K_.Data)-273.15);
strSHa3 = sprintf('  SHa3\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_SH__K_.T_EMaschine_SH_a_3__K_.Data)-273.15);
strRYO3 = sprintf('  RYO3\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_RYO__K_.T_EMaschine_RYO_3__K_.Data)-273.15);
strRTO3 = sprintf('  RTO3\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_RTO__K_.T_EMaschine_RTO_3__K_.Data)-273.15);
strAGP3 = sprintf('  AGP3\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_AGP__K_.T_EMaschine_AGP_3__K_.Data)-273.15);
strSTO3 = sprintf('  STO3\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_STO__K_.T_EMaschine_STO_3__K_.Data)-273.15);
strSYO3 = sprintf('  SYO3\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_SYO__K_.T_EMaschine_SYO_3__K_.Data)-273.15);
strHOa3 = sprintf('  HOa3\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_HO__K_.T_EMaschine_HO_a_3__K_.Data)-273.15);
strWC = sprintf('  WC\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_Kuehlmittel__K_.Data)-273.15);
strRSL3 = sprintf('    RSL3\n    %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_RSL__K_.T_EMaschine_RSL_3__K_.Data)-273.15);
strSSL3 = sprintf('    SSL3\n    %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_SSL__K_.T_EMaschine_SSL_3__K_.Data)-273.15);
strSHa4 = sprintf('  SHa4\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_SH__K_.T_EMaschine_SH_a_4__K_.Data)-273.15);
strRYO4 = sprintf('  RYO4\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_RYO__K_.T_EMaschine_RYO_4__K_.Data)-273.15);
strRTO4 = sprintf('  RTO4\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_RTO__K_.T_EMaschine_RTO_4__K_.Data)-273.15);
strAGP4 = sprintf('  AGP4\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_AGP__K_.T_EMaschine_AGP_4__K_.Data)-273.15);
strSTO4 = sprintf('  STO4\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_STO__K_.T_EMaschine_STO_4__K_.Data)-273.15);
strSYO4 = sprintf('  SYO4\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_SYO__K_.T_EMaschine_SYO_4__K_.Data)-273.15);
strHOa4 = sprintf('  HOa4\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_HO__K_.T_EMaschine_HO_a_4__K_.Data)-273.15);
strRSL4 = sprintf('    RSL4\n    %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_RSL__K_.T_EMaschine_RSL_4__K_.Data)-273.15);
strSSL4 = sprintf('    SSL4\n    %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_SSL__K_.T_EMaschine_SSL_4__K_.Data)-273.15);
strSHa5 = sprintf('  SHa5\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_SH__K_.T_EMaschine_SH_a_5__K_.Data)-273.15);
strRYO5 = sprintf('  RYO5\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_RYO__K_.T_EMaschine_RYO_5__K_.Data)-273.15);
strRTO5 = sprintf('  RTO5\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_RTO__K_.T_EMaschine_RTO_5__K_.Data)-273.15);
strAGP5 = sprintf('  AGP5\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_AGP__K_.T_EMaschine_AGP_5__K_.Data)-273.15);
strSTO5 = sprintf('  STO5\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_STO__K_.T_EMaschine_STO_5__K_.Data)-273.15);
strSYO5 = sprintf('  SYO5\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_SYO__K_.T_EMaschine_SYO_5__K_.Data)-273.15);
strHOa5 = sprintf('  HOa5\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_HO__K_.T_EMaschine_HO_a_5__K_.Data)-273.15);
strRSL5 = sprintf('    RSL5\n    %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_RSL__K_.T_EMaschine_RSL_5__K_.Data)-273.15);
strSSL5 = sprintf('    SSL5\n    %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_SSL__K_.T_EMaschine_SSL_5__K_.Data)-273.15);
strER2 = sprintf('  ER2\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_ER__K_.T_EMaschine_ER_2__K_.Data)-273.15);
strWH2 = sprintf('  WH2\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_WH__K_.T_EMaschine_WH_2__K_.Data)-273.15);
strSHp2 = sprintf('  SHp2\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_SH__K_.T_EMaschine_SH_p_2__K_.Data)-273.15);
strIA2 = sprintf('  IA2\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_IA__K_.T_EMaschine_IA_2__K_.Data)-273.15);
strHOp2 = sprintf('  HOp2\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_HO__K_.T_EMaschine_HO_p_2__K_.Data)-273.15);
strBR2 = sprintf('  BR2\n  ');
strEC2 = sprintf('  EC2\n  %.1f °C',max(EMaschine_Thermik_Output.T_EMaschine_EC__K_.T_EMaschine_EC_2__K_.Data)-273.15);

xC=[1,1,2,2,2,3,3,4,4,5,5,5,5,5,5,5,6,6,7,7,7,7,7,7,7,8,8,9,9,9,9,9,9,9,9,10,10,11,11,11,11,11,11,11,12,12,13,13,13,13,13,13,13,14,14,15,15,15,16,16];
yC=[2,5,1,5,9,4,6,4,6,1,2,3,5,7,8,9,4,6,1,2,3,5,7,8,9,4,6,1,2,3,5,7,8,9,10,4,6,1,2,3,5,7,8,9,4,6,1,2,3,5,7,8,9,4,6,1,5,9,2,5];
labels={strBR1,strEC1,strSHp1,strIA1,strHOp1,strER1,strWH1,strRSL1,strSSL1,strSHa1,strRYO1,strRTO1,strAGP1,strSTO1,strSYO1,strHOa1,strRSL2,strSSL2,strSHa2,strRYO2,strRTO2,strAGP2,strSTO2,strSYO2,strHOa2,strRSL3,strSSL3,strSHa3,strRYO3,strRTO3,strAGP3,strSTO3,strSYO3,strHOa3,strWC,strRSL4,strSSL4,strSHa4,strRYO4,strRTO4,strAGP4,strSTO4,strSYO4,strHOa4,strRSL5,strSSL5,strSHa5,strRYO5,strRTO5,strAGP5,strSTO5,strSYO5,strHOa5,strER2,strWH2,strSHp2,strIA2,strHOp2,strBR2,strEC2};

figure('Name','Thermisches Netzwerk','Color','white');
p=plot(xC,yC,'o');
p.MarkerFaceColor='b';
p.Color='b';
axis off
text(xC,yC,labels,'VerticalAlignment','middle','HorizontalAlignment','left')
title(Titel,'Position',[2,10],'HorizontalAlignment','left')
xlim([1,16])

% Linien
hold on
aussen_x=[1,1,16,16,1];
aussen_y=[1,9,9,1,1];
aussen=plot(aussen_x,aussen_y,'b');

radial1_x=[5,5];
radial1_y=[1,10];
radial1=plot(radial1_x,radial1_y,'b');
radial2_x=[7,7];
radial2_y=[1,10];
radial2=plot(radial2_x,radial2_y,'b');
radial3_x=[9,9];
radial3_y=[1,10];
radial3=plot(radial3_x,radial3_y,'b');
radial4_x=[11,11];
radial4_y=[1,10];
radial4=plot(radial4_x,radial4_y,'b');
radial5_x=[13,13];
radial5_y=[1,10];
radial5=plot(radial5_x,radial5_y,'b');

ssl1_x=[4,4,5,5,4];
ssl1_y=[6,7,8,7,6];
ssl1=plot(ssl1_x,ssl1_y,'b');
ssl2_x=[6,6,7,7,6];
ssl2_y=[6,7,8,7,6];
ssl2=plot(ssl2_x,ssl2_y,'b');
ssl3_x=[8,8,9,9,8];
ssl3_y=[6,7,8,7,6];
ssl3=plot(ssl3_x,ssl3_y,'b');
ssl4_x=[10,10,11,11,10];
ssl4_y=[6,7,8,7,6];
ssl4=plot(ssl4_x,ssl4_y,'b');
ssl5_x=[12,12,13,13,12];
ssl5_y=[6,7,8,7,6];
ssl5=plot(ssl5_x,ssl5_y,'b');

wh_x=[2,15];
wh_y=[6,6];
wh=plot(wh_x,wh_y,'b');

rsl1_x=[4,4,5,5,4];
rsl1_y=[4,3,2,3,4];
rsl1=plot(rsl1_x,rsl1_y,'b');
rsl2_x=[6,6,7,7,6];
rsl2_y=[4,3,2,3,4];
rsl2=plot(rsl2_x,rsl2_y,'b');
rsl3_x=[8,8,9,9,8];
rsl3_y=[4,3,2,3,4];
rsl3=plot(rsl3_x,rsl3_y,'b');
rsl4_x=[10,10,11,11,10];
rsl4_y=[4,3,2,3,4];
rsl4=plot(rsl4_x,rsl4_y,'b');
rsl5_x=[12,12,13,13,12];
rsl5_y=[4,3,2,3,4];
rsl5=plot(rsl5_x,rsl5_y,'b');

er_x=[2,15];
er_y=[4,4];
er=plot(er_x,er_y,'b');

agp1_x=[5,5];
agp1_y=[3,7];
agp1=plot(agp1_x,agp1_y,'--w');
agp2_x=[7,7];
agp2_y=[3,7];
agp2=plot(agp2_x,agp2_y,'--w');
agp3_x=[9,9];
agp3_y=[3,7];
agp3=plot(agp3_x,agp3_y,'--w');
agp4_x=[11,11];
agp4_y=[3,7];
agp4=plot(agp4_x,agp4_y,'--w');
agp5_x=[13,13];
agp5_y=[3,7];
agp5=plot(agp5_x,agp5_y,'--w');

agp1b_x=[4,5,4];
agp1b_y=[4,5,6];
agp1b=plot(agp1b_x,agp1b_y,'b');
agp1b=plot(agp1b_x,agp1b_y,'--w');
agp2b_x=[6,7,6];
agp2b_y=[4,5,6];
agp2b=plot(agp2b_x,agp2b_y,'b');
agp2b=plot(agp2b_x,agp2b_y,'--w');
agp3b_x=[8,9,8];
agp3b_y=[4,5,6];
agp3b=plot(agp3b_x,agp3b_y,'b');
agp3b=plot(agp3b_x,agp3b_y,'--w');
agp4b_x=[10,11,10];
agp4b_y=[4,5,6];
agp4b=plot(agp4b_x,agp4b_y,'b');
agp4b=plot(agp4b_x,agp4b_y,'--w');
agp5b_x=[12,13,12];
agp5b_y=[4,5,6];
agp5b=plot(agp5b_x,agp5b_y,'b');
agp5b=plot(agp5b_x,agp5b_y,'--w');

ia1a_x=[2,2,1];
ia1a_y=[1,5,5];
ia1a=plot(ia1a_x,ia1a_y,'b');
ia1a=plot(ia1a_x,ia1a_y,'--w');
ia1b_x=[2,2];
ia1b_y=[5,9];
ia1b=plot(ia1b_x,ia1b_y,'b');
ia1b=plot(ia1b_x,ia1b_y,'--w');
ia1c_x=[2,3];
ia1c_y=[6,6];
ia1c=plot(ia1c_x,ia1c_y,'--w');
ia1d_x=[2,3];
ia1d_y=[4,4];
ia1d=plot(ia1d_x,ia1d_y,'--w');

ia2a_x=[15,15,16];
ia2a_y=[1,5,5];
ia2a=plot(ia2a_x,ia2a_y,'b');
ia2a=plot(ia2a_x,ia2a_y,'--w');
ia2b_x=[15,15];
ia2b_y=[5,9];
ia2b=plot(ia2b_x,ia2b_y,'b');
ia2b=plot(ia2b_x,ia2b_y,'--w');
ia2c_x=[15,14];
ia2c_y=[6,6];
ia2c=plot(ia2c_x,ia2c_y,'--w');
ia2d_x=[2,3];
ia2d_y=[4,4];
ia2d=plot(ia2d_x,ia2d_y,'--w');

wc_x=[2,2,15,15];
wc_y=[9,10,10,9];
wc=plot(wc_x,wc_y,'b');

clearvars aussen* radial* ssl* wh* rsl* er* agp* ia* wc* str*