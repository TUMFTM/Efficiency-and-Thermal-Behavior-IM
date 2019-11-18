% Skript um Geschwindikeitsverlauf aus Messfahrten in Simulationsinput
% ("double time series data type") umzuwandeln

%% Lade Messfahrt:
[file,path] = uigetfile({'*.mat'},...
    'Messfahrt auswählen');
addpath(path);
load(file);
fullf=fullfile(path,file);

%% Daten Beginn finden, auf einheitliches t interpolieren und auf einheitliche Länge kürzen
% vehicle.accelerator
ind=0;
c=1; % Laufvariable
while ind == 0 && c < size(vehicle.accelerator,1)
    if vehicle.accelerator(c,2) > 0.5 && vehicle.accelerator(c,2) < 101
        ind = 1;
        Start = c;
    else
        c = c+1;
    end
end
vehicle.accelerator(1:Start-1,:) = [];
Start_T = vehicle.accelerator(1,1);
vehicle.accelerator(:,1)=vehicle.accelerator(:,1)-Start_T; % + letzter Wert aus ...
samplingtime = vehicle.accelerator(end,1)/(size(vehicle.accelerator,1)-1);

% em.r.speed
while em.r.speed(2,1) <= Start_T
    em.r.speed(1,:) = [];
end
em.r.speed(:,1) = em.r.speed(:,1)-Start_T;
ip1=interp1(em.r.speed(:,1), em.r.speed(:,2), vehicle.accelerator(:,1));
em.r.speed = horzcat(vehicle.accelerator(:,1),ip1);
if isnan(em.r.speed(end,2))
    while isnan(em.r.speed(end,2))
        em.r.speed(end,:) = [];
    end
end

% em.r.torque
    while em.r.torque(2,1) <= Start_T
        em.r.torque(1,:)=[];
    end
    em.r.torque(:,1) = em.r.torque(:,1)-Start_T;
    ip1=interp1(em.r.torque(:,1), em.r.torque(:,2), vehicle.accelerator(:,1));
    em.r.torque = horzcat(vehicle.accelerator(:,1),ip1);
    if isnan(em.r.torque(end,2))
        while isnan(em.r.torque(end,2))
            em.r.torque(end,:) = [];
        end
    end
    
% em.r.temp.stator
while em.r.temp.stator(2,1) <= Start_T
    em.r.temp.stator(1,:) = [];
end
em.r.temp.stator(:,1) = em.r.temp.stator(:,1)-Start_T;
ip1=interp1(em.r.temp.stator(:,1), em.r.temp.stator(:,2), vehicle.accelerator(:,1));
em.r.temp.stator = horzcat(vehicle.accelerator(:,1),ip1);
if isnan(em.r.temp.stator(end,2))
    while isnan(em.r.temp.stator(end,2))
        em.r.temp.stator(end,:) = [];
    end
end

% em.r.temp.coolant
while em.r.temp.coolant(2,1) <= Start_T
    em.r.temp.coolant(1,:) = [];
end
em.r.temp.coolant(:,1) = em.r.temp.coolant(:,1)-Start_T;
ip1=interp1(em.r.temp.coolant(:,1), em.r.temp.coolant(:,2), vehicle.accelerator(:,1));
em.r.temp.coolant = horzcat(vehicle.accelerator(:,1),ip1);
if isnan(em.r.temp.coolant(end,2))
    while isnan(em.r.temp.coolant(end,2))
        em.r.temp.coolant(end,:) = [];
    end
end

% Ablängen aller Variablen auf Länge der Kürzesten
ende=min(length(vehicle.accelerator),min(length(em.r.speed),min(length(em.r.torque),min(length(em.r.temp.stator),length(em.r.temp.coolant)))));
vehicle.accelerator(ende+1:end,:) = [];
em.r.speed(ende+1:end,:) = [];
em.r.torque(ende+1:end,:) = [];
em.r.temp.stator(ende+1:end,:) = [];
em.r.temp.coolant(ende+1:end,:) = [];

%% Glättung
windowspeed = 30;
em.r.speedsmooth = em.r.speed;
em.r.speedsmooth(:,2) = smoothdata(em.r.speed(:,2),'gaussian',windowspeed);

windowtorque = 5;
em.r.torquesmooth = em.r.torque;
em.r.torquesmooth(:,2) = smoothdata(em.r.torque(:,2),'gaussian',windowtorque);

windowtempstator = 300;
em.r.temp.statorsmooth = em.r.temp.stator;
em.r.temp.statorsmooth(:,2) = smoothdata(em.r.temp.stator(:,2),'gaussian', windowtempstator);

%% Berechnung mechanische Leitung EM
% zum Vergleich von Simulation und Messung
em.r.pmechsmooth = em.r.speedsmooth(:,1);
em.r.pmechsmooth(:,2) = em.r.torquesmooth(:,2).*em.r.speedsmooth(:,2)*(2*pi/60/1000); % in kW
pmech_meas=timeseries(em.r.pmechsmooth(:,2),em.r.pmechsmooth(:,1));

%% Umrechnung in Fahrzeuggeschwindigkeit in km/h
v_Fahrzeug = timeseries((em.r.speedsmooth(:,2).*(2*pi*wheel.rl.radius/gears.r.ratio*3.6/60)),em.r.speedsmooth(:,1));

% Starttemperatur der EM (Annahme, dass alle Komponenten in EM gleichmäßig
% durchgewärmt sind --> gilt nur nach ausreichend langer Standzeit)
T_Start_EMaschine = 273.15 + em.r.temp.stator(1,2);
% Ausgaben zum Verlgleich Simulation & Messung
T_Verlauf_EMaschine = em.r.temp.statorsmooth;
T_Verlauf_Coolant = em.r.temp.coolant;

clearvars aux bat c em ende file fullf gears ind inv ip1 path samplingtime Start Start_T vehicle wheel window*