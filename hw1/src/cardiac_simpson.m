% Caldiac functions with modified Simpson's rule

%% (1). Load sample video, get cm/pixel unit-transformation
f = "../data/Apical4chVolunteer.mp4";
xrange = [100, 200];
yrange = [550, 600];
[height, width, scale] = calibrate(f, xrange, yrange);
clear f xrange yrange
%}

%% (2). Measure ventricular diameters (ED, ES) in each videos
files = "../data/" + [
         "Apical4chVolunteer.mp4"   % Apical view
         "SA_ApexVolunteer.mp4"     % Apex
         "SA_MitralVolunteer.mp4"   % MV
         "SA_PapilaryVolunteer.mp4" % PV
        ];

ed_diams = struct();
es_diams = struct();

for i = 1:length(files)
    fig = figure(i);
    f = files(i);
    [edd, esd] = measureDiam(f, height, width, scale);
    
    if contains(f, "Apical")
        ed_diams.long_axis = edd;
        es_diams.long_axis = esd;
    elseif contains(f, "Apex")
        ed_diams.apex = edd;
        es_diams.apex = esd;
    elseif contains(f, "Mitral")
        ed_diams.mv = edd;
        es_diams.mv = esd;
    else  % PV view
        ed_diams.pm = edd;
        es_diams.pm = esd;
    end
    
    close(fig);
end


%% (3). Calculate Volumes
[ED_vol, ES_vol, stroke_vol] = calcStrokeVol(ed_diams, es_diams);
EF = calcEjectFrac(stroke_vol, ED_vol);
CO = calcCardiacOutput(stroke_vol);


%% (4). Save to output
of1 = "../results/diameters.txt";
Views = {'Mitral'; 'Papillary'; 'Apex'; 'Apical'};
Diastolic = {ed_diams.mv; ed_diams.pm; ed_diams.apex; ed_diams.long_axis};
Systolic = {es_diams.mv; es_diams.pm; es_diams.apex; es_diams.long_axis};

t1 = table(Views, Diastolic, Systolic);
writetable(t1, of1);

of2 = "../results/volumes.txt";
Types = {'EDV'; 'ESV'; 'SV'; 'EF'; 'CO'};
Volumes = {ED_vol; ES_vol; stroke_vol; EF; CO};

t2 = table(Types, Volumes);
writetable(t2, of2);


%% util functions
function A = calcArea(diam)
    A = pi * diam^2 / 4;
end


function V = calcSimpVol(diams)
    A_mv = calcArea(diams.mv);
    A_pm = calcArea(diams.pm);
    A_ap = calcArea(diams.apex);
    L = diams.long_axis;
    
    % Ventricular volume by modified Simpson's eqn.
    V = (A_mv+A_pm)*L/3 + A_ap*L/6 + pi*L^3/162;
    
end


function [EDV, ESV, SV] = calcStrokeVol(ed_diams, es_diams)
    EDV = calcSimpVol(ed_diams);
    ESV = calcSimpVol(es_diams);
    
    assert (EDV >= ESV, "diastolic volumn should >= systolic vol.");
    
    SV = EDV - ESV; % stroke vol.
end


function EF = calcEjectFrac(SV, EDV)
    EF = SV / EDV; % ejection fraction
end


function CO = calcCardiacOutput(SV)
    HR = 60; % heart rate 
    CO = SV * HR; % cardiac output
end
