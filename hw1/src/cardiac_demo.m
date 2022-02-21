% Cardiac functions with the demo. measure

%% (1). "Hard-coded" values from demo video's measurement
L_ed = 8.40;
L_es = 7.14;
D_ed = 5.40;
D_es = 3.61;
VTI = 23.8;
D_ao = 2.0;
HR = 60;


%% (2). Calculation Volumes
EDV = pi/6 * D_ed^2 * L_ed;
ESV = pi/6 * D_es^2 * L_es;

SV = EDV - ESV;
EF = SV / EDV;
CO = SV * HR;
CO_doppler = VTI * HR * (pi*D_ao^2/4);


%% (3). Save to output
of = "../results/volumes_demo.txt";
Types = {'EDV'; 'ESV'; 'SV'; 'EF'; 'CO'; 'CO_doppler'};
Volumes = {EDV; ESV; SV; EF; CO; CO_doppler};

t = table(Types, Volumes);
writetable(t, of);



