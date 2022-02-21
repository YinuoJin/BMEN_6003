% Calculate balloon volumes

%% (1). Measure with the larger machine (IE33)
L_1 = [9.30, 10.20];
D1_1 = [8.27, 8.64];
D2_1 = [5.79, 5.19];

se_L = std(L_1) / sqrt(length(L_1));
se_D1 = std(D1_1) / sqrt(length(D1_1));
se_D2 = std(D2_1) / sqrt(length(D2_1));

V1 = calcVol(mean(L_1), mean(D1_1), mean(D2_1));
se_1 = calcVol(se_L, se_D1, se_D2);


%% (2). Measure with the smaller machine (t3000)
L_2 = 9.40;
D1_2 = 8.24;
D2_2 = 5.34;

V2 = calcVol(L_2, D1_2, D2_2);

%% (3). Save to output
of = "../results/balloon.txt";
Types = {'Machine_1'; 'Machine_2'; 'M1_SE'};
Volumes = {V1; V2; se_1};
t = table(Types, Volumes);
writetable(t, of);


function V = calcVol(l, d1, d2)
    V = 4/3*pi * l/2 * d1/2 * d2/2;
end

