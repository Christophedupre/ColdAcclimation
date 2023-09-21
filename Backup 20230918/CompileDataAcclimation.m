% This function compiles into 1 file all the experimental data.
function CompileDataAcclimation();
disp('Compiling all the experimental data into "data.mat"');
%%% FULL EXPERIMENTAL DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTROLS
PosCtrl = [readtable('Acclimation28D.csv').AcclimA1 readtable('Acclimation28D.csv').AcclimA2 readtable('Acclimation28D.csv').AcclimA3 readtable('Acclimation28D.csv').AcclimA4 readtable('Acclimation28D.csv').AcclimA5 readtable('Acclimation28D.csv').AcclimA6 readtable('Acclimation28D.csv').AcclimA7 readtable('Acclimation28D.csv').AcclimA8 readtable('Acclimation28D.csv').AcclimA9 readtable('Acclimation28D.csv').AcclimA10 readtable('Acclimation28D.csv').AcclimA11 readtable('Acclimation28D.csv').AcclimA12 readtable('Acclimation28D.csv').AcclimA13 readtable('Acclimation28D.csv').AcclimA14 readtable('Acclimation28D.csv').AcclimA15]./5;
NoPCR = [readtable('Acclimation28D.csv').NoPCRA1 readtable('Acclimation28D.csv').NoPCRA2 readtable('Acclimation28D.csv').NoPCRA3 readtable('Acclimation28D.csv').NoPCRA4 readtable('Acclimation28D.csv').NoPCRA5 readtable('Acclimation28D.csv').NoPCRA6 readtable('Acclimation28D.csv').NoPCRA7 readtable('Acclimation28D.csv').NoPCRA8 readtable('Acclimation28D.csv').NoPCRA9 readtable('Acclimation28D.csv').NoPCRA10 readtable('Acclimation28D.csv').NoPCRA11 readtable('Acclimation28D.csv').NoPCRA12 readtable('Acclimation28D.csv').NoPCRA13 readtable('Acclimation28D.csv').NoPCRA14 readtable('Acclimation28D.csv').NoPCRA15]./5;
JustPCR = [readtable('Acclimation28D.csv').JustPCRA1 readtable('Acclimation28D.csv').JustPCRA2 readtable('Acclimation28D.csv').JustPCRA3 readtable('Acclimation28D.csv').JustPCRA4 readtable('Acclimation28D.csv').JustPCRA5 readtable('Acclimation28D.csv').JustPCRA6 readtable('Acclimation28D.csv').JustPCRA7 readtable('Acclimation28D.csv').JustPCRA8 readtable('Acclimation28D.csv').JustPCRA9 readtable('Acclimation28D.csv').JustPCRA10 readtable('Acclimation28D.csv').JustPCRA11 readtable('Acclimation28D.csv').JustPCRA12 readtable('Acclimation28D.csv').JustPCRA13 readtable('Acclimation28D.csv').JustPCRA14 readtable('Acclimation28D.csv').JustPCRA15]./5;
Pulses0_2 = [readtable('Acclimation28D.csv').P02A1 readtable('Acclimation28D.csv').P02A2 readtable('Acclimation28D.csv').P02A3 readtable('Acclimation28D.csv').P02A4 readtable('Acclimation28D.csv').P02A5 readtable('Acclimation28D.csv').P02A6 readtable('Acclimation28D.csv').P02A7 readtable('Acclimation28D.csv').P02A8 readtable('Acclimation28D.csv').P02A9 readtable('Acclimation28D.csv').P02A10 readtable('Acclimation28D.csv').P02A11 readtable('Acclimation28D.csv').P02A12 readtable('Acclimation28D.csv').P02A13 readtable('Acclimation28D.csv').P02A14 readtable('Acclimation28D.csv').P02A15]./5;;
Pulses1_1 = [readtable('Acclimation28D.csv').P11A1 readtable('Acclimation28D.csv').P11A2 readtable('Acclimation28D.csv').P11A3 readtable('Acclimation28D.csv').P11A4 readtable('Acclimation28D.csv').P11A5 readtable('Acclimation28D.csv').P11A6 readtable('Acclimation28D.csv').P11A7 readtable('Acclimation28D.csv').P11A8 readtable('Acclimation28D.csv').P11A9 readtable('Acclimation28D.csv').P11A10 readtable('Acclimation28D.csv').P11A11 readtable('Acclimation28D.csv').P11A12 readtable('Acclimation28D.csv').P11A13 readtable('Acclimation28D.csv').P11A14 readtable('Acclimation28D.csv').P11A15]./5;;
% RAMP
Ramp = [readtable('Acclimation28D.csv').RampA1 readtable('Acclimation28D.csv').RampA2 readtable('Acclimation28D.csv').RampA3 readtable('Acclimation28D.csv').RampA4 readtable('Acclimation28D.csv').RampA5 readtable('Acclimation28D.csv').RampA6 readtable('Acclimation28D.csv').RampA7 readtable('Acclimation28D.csv').RampA8 readtable('Acclimation28D.csv').RampA9 readtable('Acclimation28D.csv').RampA10 readtable('Acclimation28D.csv').RampA11 readtable('Acclimation28D.csv').RampA12 readtable('Acclimation28D.csv').RampA13 readtable('Acclimation28D.csv').RampA14 readtable('Acclimation28D.csv').RampA15]./5;;
% CONSTANT
Cst4C = [readtable('Acclimation28D.csv').Cst4CA1 readtable('Acclimation28D.csv').Cst4CA2 readtable('Acclimation28D.csv').Cst4CA3 readtable('Acclimation28D.csv').Cst4CA4 readtable('Acclimation28D.csv').Cst4CA5 readtable('Acclimation28D.csv').Cst4CA6 readtable('Acclimation28D.csv').Cst4CA7 readtable('Acclimation28D.csv').Cst4CA8 readtable('Acclimation28D.csv').Cst4CA9 readtable('Acclimation28D.csv').Cst4CA10 readtable('Acclimation28D.csv').Cst4CA11 readtable('Acclimation28D.csv').Cst4CA12 readtable('Acclimation28D.csv').Cst4CA13 readtable('Acclimation28D.csv').Cst4CA14 readtable('Acclimation28D.csv').Cst4CA15];
Cst6C = [readtable('Acclimation28D.csv').Cst6CA1 readtable('Acclimation28D.csv').Cst6CA2 readtable('Acclimation28D.csv').Cst6CA3 readtable('Acclimation28D.csv').Cst6CA4 readtable('Acclimation28D.csv').Cst6CA5 readtable('Acclimation28D.csv').Cst6CA6 readtable('Acclimation28D.csv').Cst6CA7 readtable('Acclimation28D.csv').Cst6CA8 readtable('Acclimation28D.csv').Cst6CA9 readtable('Acclimation28D.csv').Cst6CA10 readtable('Acclimation28D.csv').Cst6CA11 readtable('Acclimation28D.csv').Cst6CA12 readtable('Acclimation28D.csv').Cst6CA13 readtable('Acclimation28D.csv').Cst6CA14 readtable('Acclimation28D.csv').Cst6CA15];
Cst8C = [readtable('Acclimation28D.csv').Cst8CA1 readtable('Acclimation28D.csv').Cst8CA2 readtable('Acclimation28D.csv').Cst8CA3 readtable('Acclimation28D.csv').Cst8CA4 readtable('Acclimation28D.csv').Cst8CA5 readtable('Acclimation28D.csv').Cst8CA6 readtable('Acclimation28D.csv').Cst8CA7 readtable('Acclimation28D.csv').Cst8CA8 readtable('Acclimation28D.csv').Cst8CA9 readtable('Acclimation28D.csv').Cst8CA10 readtable('Acclimation28D.csv').Cst8CA11 readtable('Acclimation28D.csv').Cst8CA12 readtable('Acclimation28D.csv').Cst8CA13 readtable('Acclimation28D.csv').Cst8CA14 readtable('Acclimation28D.csv').Cst8CA15];
Cst10C = [readtable('Acclimation28D.csv').Cst10CA1 readtable('Acclimation28D.csv').Cst10CA2 readtable('Acclimation28D.csv').Cst10CA3 readtable('Acclimation28D.csv').Cst10CA4 readtable('Acclimation28D.csv').Cst10CA5 readtable('Acclimation28D.csv').Cst10CA6 readtable('Acclimation28D.csv').Cst10CA7 readtable('Acclimation28D.csv').Cst10CA8 readtable('Acclimation28D.csv').Cst10CA9 readtable('Acclimation28D.csv').Cst10CA10 readtable('Acclimation28D.csv').Cst10CA11 readtable('Acclimation28D.csv').Cst10CA12 readtable('Acclimation28D.csv').Cst10CA13 readtable('Acclimation28D.csv').Cst10CA14 readtable('Acclimation28D.csv').Cst10CA15];
Cst12C = [readtable('Acclimation28D.csv').Cst12CA1 readtable('Acclimation28D.csv').Cst12CA2 readtable('Acclimation28D.csv').Cst12CA3 readtable('Acclimation28D.csv').Cst12CA4 readtable('Acclimation28D.csv').Cst12CA5 readtable('Acclimation28D.csv').Cst12CA6 readtable('Acclimation28D.csv').Cst12CA7 readtable('Acclimation28D.csv').Cst12CA8 readtable('Acclimation28D.csv').Cst12CA9 readtable('Acclimation28D.csv').Cst12CA10 readtable('Acclimation28D.csv').Cst12CA11 readtable('Acclimation28D.csv').Cst12CA12 readtable('Acclimation28D.csv').Cst12CA13 readtable('Acclimation28D.csv').Cst12CA14 readtable('Acclimation28D.csv').Cst12CA15];
Cst14C = [readtable('Acclimation28D.csv').Cst14CA1 readtable('Acclimation28D.csv').Cst14CA2 readtable('Acclimation28D.csv').Cst14CA3 readtable('Acclimation28D.csv').Cst14CA4 readtable('Acclimation28D.csv').Cst14CA5 readtable('Acclimation28D.csv').Cst14CA6 readtable('Acclimation28D.csv').Cst14CA7 readtable('Acclimation28D.csv').Cst14CA8 readtable('Acclimation28D.csv').Cst14CA9 readtable('Acclimation28D.csv').Cst14CA10 readtable('Acclimation28D.csv').Cst14CA11 readtable('Acclimation28D.csv').Cst14CA12 readtable('Acclimation28D.csv').Cst14CA13 readtable('Acclimation28D.csv').Cst14CA14 readtable('Acclimation28D.csv').Cst14CA15];
Cst16C = [readtable('Acclimation28D.csv').Cst16CA1 readtable('Acclimation28D.csv').Cst16CA2 readtable('Acclimation28D.csv').Cst16CA3 readtable('Acclimation28D.csv').Cst16CA4 readtable('Acclimation28D.csv').Cst16CA5 readtable('Acclimation28D.csv').Cst16CA6 readtable('Acclimation28D.csv').Cst16CA7 readtable('Acclimation28D.csv').Cst16CA8 readtable('Acclimation28D.csv').Cst16CA9 readtable('Acclimation28D.csv').Cst16CA10 readtable('Acclimation28D.csv').Cst16CA11 readtable('Acclimation28D.csv').Cst16CA12 readtable('Acclimation28D.csv').Cst16CA13 readtable('Acclimation28D.csv').Cst16CA14 readtable('Acclimation28D.csv').Cst16CA15];
Cst18C = [readtable('Acclimation28D.csv').Cst18CA1 readtable('Acclimation28D.csv').Cst18CA2 readtable('Acclimation28D.csv').Cst18CA3 readtable('Acclimation28D.csv').Cst18CA4 readtable('Acclimation28D.csv').Cst18CA5 readtable('Acclimation28D.csv').Cst18CA6 readtable('Acclimation28D.csv').Cst18CA7 readtable('Acclimation28D.csv').Cst18CA8 readtable('Acclimation28D.csv').Cst18CA9 readtable('Acclimation28D.csv').Cst18CA10 readtable('Acclimation28D.csv').Cst18CA11 readtable('Acclimation28D.csv').Cst18CA12 readtable('Acclimation28D.csv').Cst18CA13 readtable('Acclimation28D.csv').Cst18CA14 readtable('Acclimation28D.csv').Cst18CA15];
Cst20C = [readtable('Acclimation28D.csv').P02A1 readtable('Acclimation28D.csv').P02A2 readtable('Acclimation28D.csv').P02A3 readtable('Acclimation28D.csv').P02A4 readtable('Acclimation28D.csv').P02A5 readtable('Acclimation28D.csv').P02A6 readtable('Acclimation28D.csv').P02A7 readtable('Acclimation28D.csv').P02A8 readtable('Acclimation28D.csv').P02A9 readtable('Acclimation28D.csv').P02A10 readtable('Acclimation28D.csv').P02A11 readtable('Acclimation28D.csv').P02A12 readtable('Acclimation28D.csv').P02A13 readtable('Acclimation28D.csv').P02A14 readtable('Acclimation28D.csv').P02A15]; %The P02 experiment is equivalent to a Cst22C experiment
Cst22C = [readtable('Acclimation28D.csv').Cst22CA1 readtable('Acclimation28D.csv').Cst22CA2 readtable('Acclimation28D.csv').Cst22CA3 readtable('Acclimation28D.csv').Cst22CA4 readtable('Acclimation28D.csv').Cst22CA5 readtable('Acclimation28D.csv').Cst22CA6 readtable('Acclimation28D.csv').Cst22CA7 readtable('Acclimation28D.csv').Cst22CA8 readtable('Acclimation28D.csv').Cst22CA9 readtable('Acclimation28D.csv').Cst22CA10 readtable('Acclimation28D.csv').Cst22CA11 readtable('Acclimation28D.csv').Cst22CA12 readtable('Acclimation28D.csv').Cst22CA13 readtable('Acclimation28D.csv').Cst22CA14 readtable('Acclimation28D.csv').Cst22CA15];
Cst24C = [readtable('Acclimation28D.csv').Cst24CA1 readtable('Acclimation28D.csv').Cst24CA2 readtable('Acclimation28D.csv').Cst24CA3 readtable('Acclimation28D.csv').Cst24CA4 readtable('Acclimation28D.csv').Cst24CA5 readtable('Acclimation28D.csv').Cst24CA6 readtable('Acclimation28D.csv').Cst24CA7 readtable('Acclimation28D.csv').Cst24CA8 readtable('Acclimation28D.csv').Cst24CA9 readtable('Acclimation28D.csv').Cst24CA10 readtable('Acclimation28D.csv').Cst24CA11 readtable('Acclimation28D.csv').Cst24CA12 readtable('Acclimation28D.csv').Cst24CA13 readtable('Acclimation28D.csv').Cst24CA14 readtable('Acclimation28D.csv').Cst24CA15];
Cst26C = [readtable('Acclimation28D.csv').Cst26CA1 readtable('Acclimation28D.csv').Cst26CA2 readtable('Acclimation28D.csv').Cst26CA3 readtable('Acclimation28D.csv').Cst26CA4 readtable('Acclimation28D.csv').Cst26CA5 readtable('Acclimation28D.csv').Cst26CA6 readtable('Acclimation28D.csv').Cst26CA7 readtable('Acclimation28D.csv').Cst26CA8 readtable('Acclimation28D.csv').Cst26CA9 readtable('Acclimation28D.csv').Cst26CA10 readtable('Acclimation28D.csv').Cst26CA11 readtable('Acclimation28D.csv').Cst26CA12 readtable('Acclimation28D.csv').Cst26CA13 readtable('Acclimation28D.csv').Cst26CA14 readtable('Acclimation28D.csv').Cst26CA15];
Cst28C = [readtable('Acclimation28D.csv').Cst28CA1 readtable('Acclimation28D.csv').Cst28CA2 readtable('Acclimation28D.csv').Cst28CA3 readtable('Acclimation28D.csv').Cst28CA4 readtable('Acclimation28D.csv').Cst28CA5 readtable('Acclimation28D.csv').Cst28CA6 readtable('Acclimation28D.csv').Cst28CA7 readtable('Acclimation28D.csv').Cst28CA8 readtable('Acclimation28D.csv').Cst28CA9 readtable('Acclimation28D.csv').Cst28CA10 readtable('Acclimation28D.csv').Cst28CA11 readtable('Acclimation28D.csv').Cst28CA12 readtable('Acclimation28D.csv').Cst28CA13 readtable('Acclimation28D.csv').Cst28CA14 readtable('Acclimation28D.csv').Cst28CA15];
Cst30C = [readtable('Acclimation28D.csv').Cst30CA1 readtable('Acclimation28D.csv').Cst30CA2 readtable('Acclimation28D.csv').Cst30CA3 readtable('Acclimation28D.csv').Cst30CA4 readtable('Acclimation28D.csv').Cst30CA5 readtable('Acclimation28D.csv').Cst30CA6 readtable('Acclimation28D.csv').Cst30CA7 readtable('Acclimation28D.csv').Cst30CA8 readtable('Acclimation28D.csv').Cst30CA9 readtable('Acclimation28D.csv').Cst30CA10 readtable('Acclimation28D.csv').Cst30CA11 readtable('Acclimation28D.csv').Cst30CA12 readtable('Acclimation28D.csv').Cst30CA13 readtable('Acclimation28D.csv').Cst30CA14 readtable('Acclimation28D.csv').Cst30CA15];
Cst32C = [readtable('Acclimation28D.csv').Cst32CA1 readtable('Acclimation28D.csv').Cst32CA2 readtable('Acclimation28D.csv').Cst32CA3 readtable('Acclimation28D.csv').Cst32CA4 readtable('Acclimation28D.csv').Cst32CA5 readtable('Acclimation28D.csv').Cst32CA6 readtable('Acclimation28D.csv').Cst32CA7 readtable('Acclimation28D.csv').Cst32CA8 readtable('Acclimation28D.csv').Cst32CA9 readtable('Acclimation28D.csv').Cst32CA10 readtable('Acclimation28D.csv').Cst32CA11 readtable('Acclimation28D.csv').Cst32CA12 readtable('Acclimation28D.csv').Cst32CA13 readtable('Acclimation28D.csv').Cst32CA14 readtable('Acclimation28D.csv').Cst32CA15];
Cst34C = [readtable('Acclimation28D.csv').Cst34CA1 readtable('Acclimation28D.csv').Cst34CA2 readtable('Acclimation28D.csv').Cst34CA3 readtable('Acclimation28D.csv').Cst34CA4 readtable('Acclimation28D.csv').Cst34CA5 readtable('Acclimation28D.csv').Cst34CA6 readtable('Acclimation28D.csv').Cst34CA7 readtable('Acclimation28D.csv').Cst34CA8 readtable('Acclimation28D.csv').Cst34CA9 readtable('Acclimation28D.csv').Cst34CA10 readtable('Acclimation28D.csv').Cst34CA11 readtable('Acclimation28D.csv').Cst34CA12 readtable('Acclimation28D.csv').Cst34CA13 readtable('Acclimation28D.csv').Cst34CA14 readtable('Acclimation28D.csv').Cst34CA15];
% Store all cst experiments in a 3D matrix, and normalize
CstData3D = cat(3, Cst4C, Cst6C, Cst8C, Cst10C, Cst12C, Cst14C, Cst16C,...
    Cst18C, Cst20C, Cst22C, Cst24C, Cst26C, Cst28C, Cst30C, Cst32C, Cst34C);
CstData3D = CstData3D./5; %Normalize
% CstData3D = cat(1,nan(15,15,16),CstData3D); %Add the missing data from
% training phase (not necessary yet)
clear Cst4C Cst6C Cst8C Cst10C Cst12C Cst14C Cst16C Cst18C Cst20C Cst22C Cst24C Cst26C Cst28C Cst30C Cst32C Cst34C;

% REPEAT OF CONSTANT 10C13D over 30 animals
Cst10C13D30Reps = [readtable('Acclimation28D.csv').Cst10C13DA1 readtable('Acclimation28D.csv').Cst10C13DA2 readtable('Acclimation28D.csv').Cst10C13DA3 readtable('Acclimation28D.csv').Cst10C13DA4 readtable('Acclimation28D.csv').Cst10C13DA5 readtable('Acclimation28D.csv').Cst10C13DA6 readtable('Acclimation28D.csv').Cst10C13DA7 readtable('Acclimation28D.csv').Cst10C13DA8 readtable('Acclimation28D.csv').Cst10C13DA9 readtable('Acclimation28D.csv').Cst10C13DA10 readtable('Acclimation28D.csv').Cst10C13DA11 readtable('Acclimation28D.csv').Cst10C13DA12 readtable('Acclimation28D.csv').Cst10C13DA13 readtable('Acclimation28D.csv').Cst10C13DA14 readtable('Acclimation28D.csv').Cst10C13DA15 readtable('Acclimation28D.csv').Cst10C13DA16 readtable('Acclimation28D.csv').Cst10C13DA17 readtable('Acclimation28D.csv').Cst10C13DA18 readtable('Acclimation28D.csv').Cst10C13DA19 readtable('Acclimation28D.csv').Cst10C13DA20 readtable('Acclimation28D.csv').Cst10C13DA21 readtable('Acclimation28D.csv').Cst10C13DA22 readtable('Acclimation28D.csv').Cst10C13DA23 readtable('Acclimation28D.csv').Cst10C13DA24 readtable('Acclimation28D.csv').Cst10C13DA25 readtable('Acclimation28D.csv').Cst10C13DA26 readtable('Acclimation28D.csv').Cst10C13DA27 readtable('Acclimation28D.csv').Cst10C13DA28 readtable('Acclimation28D.csv').Cst10C13DA29 readtable('Acclimation28D.csv').Cst10C13DA30]./5;

%%% DEACCLIMATION DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = readtable('AcclimationPersistenceData.xlsx');
Deacclimation = table2array(D(:,5:end));
Deacclimation(1,:) = [];
HealthNaive = zeros(size(Deacclimation,1), 28); %Stores the first 28 days of Health

%Fist dataset
for i=1:28
    HealthNaive(i,1:28) = Deacclimation(i,i:i+27);
end
%Second dataset (naive animals)
for i=29:size(Deacclimation,1)
    HealthNaive(i,1:28) = Deacclimation(i,1:28);
end
Deacclimation = HealthNaive./5; clear D HealthNaive i;
%%% PILOT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Gradients
Pilot.AcclimGradient0_5c = rmmissing([readtable('AcclimationPilotData.csv').BeforeHab0_5 readtable('AcclimationPilotData.csv').AfterHab0_5])./5;
Pilot.AcclimGradient0_625c = rmmissing([readtable('AcclimationPilotData.csv').BeforeHab0_625 readtable('AcclimationPilotData.csv').AfterHab0_625])./5;
Pilot.AcclimGradient0_75c = rmmissing([readtable('AcclimationPilotData.csv').BeforeHab0_75 readtable('AcclimationPilotData.csv').AfterHab0_75])./5;
Pilot.AcclimGradient1c = rmmissing([readtable('AcclimationPilotData.csv').BeforeHab1 readtable('AcclimationPilotData.csv').AfterHab1])./5;
Pilot.AcclimGradient2c = rmmissing([readtable('AcclimationPilotData.csv').BeforeHab2 readtable('AcclimationPilotData.csv').AfterHab2])./5;
Pilot.NotHabituated = rmmissing([readtable('AcclimationPilotData.csv').BeforeNoHab readtable('AcclimationPilotData.csv').AfterNoHab])./5;
Pilot.dataGradientAll = [Pilot.NotHabituated(:,2) Pilot.AcclimGradient2c(:,2) Pilot.AcclimGradient1c(:,2) Pilot.AcclimGradient0_75c(:,2) Pilot.AcclimGradient0_625c(1:10,2) Pilot.AcclimGradient0_5c(1:10,2)];
%Constant
Pilot.Constant12C = rmmissing([readtable('AcclimationPilotData.csv').BeforeConst12 readtable('AcclimationPilotData.csv').AfterConst12])./5;
Pilot.Constant15C = rmmissing([readtable('AcclimationPilotData.csv').BeforeConst15 readtable('AcclimationPilotData.csv').AfterConst15])./5;
Pilot.Constant18C = rmmissing([readtable('AcclimationPilotData.csv').BeforeConst18 readtable('AcclimationPilotData.csv').AfterConst18])./5;
Pilot.Constant22C = rmmissing([readtable('AcclimationPilotData.csv').BeforeConst22 readtable('AcclimationPilotData.csv').AfterConst22])./5;
%Summary Table for Asymmetric pulses 1 Day Test: CycleType, PulseDuration, CycleDuration, TimeToAcclimation, AverageTemp, Iterations
CycleType = rmmissing(readtable('AcclimationPilotData.csv').CycleType);
PulseDuration = rmmissing(readtable('AcclimationPilotData.csv').PulseDuration);
RecoveryTime = rmmissing(readtable('AcclimationPilotData.csv').RecoveryTime);
CycleDuration = rmmissing(readtable('AcclimationPilotData.csv').CycleDuration);
Iterations = rmmissing(readtable('AcclimationPilotData.csv').Iterations);
TotalTimeAt4C = rmmissing(readtable('AcclimationPilotData.csv').TotalTimeAt4C);
AverageTemp = rmmissing(readtable('AcclimationPilotData.csv').AverageTemp);
TimeToAcclimation = rmmissing(readtable('AcclimationPilotData.csv').TimeToAcclimation);
Pilot.TableSummaryCRAS1DTest = table(CycleType, PulseDuration, RecoveryTime, CycleDuration, Iterations, TotalTimeAt4C, AverageTemp, TimeToAcclimation);
clear CycleType PulseDuration RecoveryTime CycleDuration Iterations TotalTimeAt4C AverageTemp TimeToAcclimation
%Symetric pulses (CRS = Cold Room Symmetric)
Pilot.AcclimLevelCRS24hr = rmmissing([readtable('AcclimationPilotData.csv').Iteration_CRS24hr readtable('AcclimationPilotData.csv').AcclimLevelCRS24hr./5]);
Pilot.AcclimLevelCRS8hr = rmmissing([readtable('AcclimationPilotData.csv').Iteration_CRS8hr readtable('AcclimationPilotData.csv').AcclimLevelCRS8hr./5]);
Pilot.AcclimLevelCRS1hr = rmmissing([readtable('AcclimationPilotData.csv').Iteration_CRS1hr readtable('AcclimationPilotData.csv').AcclimLevelCRS1hr./5]);
Pilot.AcclimLevelCRS0_5hr = rmmissing([readtable('AcclimationPilotData.csv').Iteration_CRS0_5hr readtable('AcclimationPilotData.csv').AcclimLevelCRS0_5hr./5]);
Pilot.AcclimLevelCRS0_1hr = rmmissing([readtable('AcclimationPilotData.csv').Iteration_CRS0_1hr readtable('AcclimationPilotData.csv').AcclimLevelCRS0_1hr./5]);
Pilot.AcclimPulsesSym = rmmissing([readtable('AcclimationPilotData.csv').PulseDurationCRS readtable('AcclimationPilotData.csv').TimeToAcclimationCRS]);
%Asymetric pulses (CRAS = Cold Room ASymmetric)
Pilot.AcclimLevelCRAS1hr1hr = rmmissing([readtable('AcclimationPilotData.csv').Iteration_CRAS1hr1hr readtable('AcclimationPilotData.csv').AcclimLevelCRAS1hr1hr./5]);
Pilot.AcclimLevelCRAS0_5hr1_5hr = rmmissing([readtable('AcclimationPilotData.csv').Iteration_CRAS0_5hr1_5hr readtable('AcclimationPilotData.csv').AcclimLevelCRAS0_5hr1_5hr./5]);
Pilot.AcclimLevelCRAS0_2hr1_8hr = rmmissing([readtable('AcclimationPilotData.csv').Iteration_CRAS0_2hr1_8hr readtable('AcclimationPilotData.csv').AcclimLevelCRAS0_2hr1_8hr./5]);
Pilot.AcclimLevelCRAS0_1hr1_9hr = rmmissing([readtable('AcclimationPilotData.csv').Iteration_CRAS0_1hr1_9hr readtable('AcclimationPilotData.csv').AcclimLevelCRAS0_1hr1_9hr./5]);
Pilot.AcclimLevelCRAS0_05hr1_95hr = rmmissing([readtable('AcclimationPilotData.csv').Iteration_CRAS0_05hr1_95hr readtable('AcclimationPilotData.csv').AcclimLevelCRAS0_05hr1_95hr./5]);
Pilot.AcclimLevelCRAS0_025hr1_975hr = rmmissing([readtable('AcclimationPilotData.csv').Iteration_CRAS0_025hr1_975hr readtable('AcclimationPilotData.csv').AcclimLevelCRAS0_025hr1_975hr./5]);
Pilot.AcclimLevelCRAS0hr2hr = rmmissing([readtable('AcclimationPilotData.csv').Iteration_CRAS0hr2hr readtable('AcclimationPilotData.csv').AcclimLevelCRAS0hr2hr./5]);
Pilot.AcclimLevelCRAS0hr2hrExtTest = rmmissing([readtable('AcclimationPilotData.csv').CRAS0hr2hrExtTestDay0 readtable('AcclimationPilotData.csv').CRAS0hr2hrExtTestDay1 readtable('AcclimationPilotData.csv').CRAS0hr2hrExtTestDay2 readtable('AcclimationPilotData.csv').CRAS0hr2hrExtTestDay3 readtable('AcclimationPilotData.csv').CRAS0hr2hrExtTestDay4 readtable('AcclimationPilotData.csv').CRAS0hr2hrExtTestDay5])./5;
Pilot.AcclimPulsesAsym = rmmissing([readtable('AcclimationPilotData.csv').PulseDurationAsym readtable('AcclimationPilotData.csv').TimeToAcclimationCRAS]);
%Regeneration and Recovery
Pilot.DecayRegen.HDN = rmmissing([readtable('AcclimationPilotData.csv').HDN0 readtable('AcclimationPilotData.csv').HDN1 readtable('AcclimationPilotData.csv').HDN2 readtable('AcclimationPilotData.csv').HDN3 readtable('AcclimationPilotData.csv').HDN4 readtable('AcclimationPilotData.csv').HDN5])./5;
Pilot.DecayRegen.HRN = rmmissing([readtable('AcclimationPilotData.csv').HRN0 readtable('AcclimationPilotData.csv').HRN1 readtable('AcclimationPilotData.csv').HRN2 readtable('AcclimationPilotData.csv').HRN3 readtable('AcclimationPilotData.csv').HRN4 readtable('AcclimationPilotData.csv').HRN5])./5;
Pilot.DecayRegen.HDA = rmmissing([readtable('AcclimationPilotData.csv').HDA0 readtable('AcclimationPilotData.csv').HDA1 readtable('AcclimationPilotData.csv').HDA2 readtable('AcclimationPilotData.csv').HDA3 readtable('AcclimationPilotData.csv').HDA4 readtable('AcclimationPilotData.csv').HDA5])./5;
Pilot.DecayRegen.HRA = rmmissing([readtable('AcclimationPilotData.csv').HRA0 readtable('AcclimationPilotData.csv').HRA1 readtable('AcclimationPilotData.csv').HRA2 readtable('AcclimationPilotData.csv').HRA3 readtable('AcclimationPilotData.csv').HRA4 readtable('AcclimationPilotData.csv').HRA5])./5;
%Temperature threshold
Pilot.TempThreshold = rmmissing([readtable('AcclimationPilotData.csv').TTEAnimal readtable('AcclimationPilotData.csv').TTETemp readtable('AcclimationPilotData.csv').TempThresholdHealthBefore./5 readtable('AcclimationPilotData.csv').TempThresholdHealthAfter./5]);
%Temperature measurement during abrupt transition
Pilot.ATTM.Data = rmmissing([readtable('AcclimationPilotData.csv').ATTMTime readtable('AcclimationPilotData.csv').ATTMDown1 readtable('AcclimationPilotData.csv').ATTMDown2 readtable('AcclimationPilotData.csv').ATTMDown3 readtable('AcclimationPilotData.csv').ATTMUp1 readtable('AcclimationPilotData.csv').ATTMUp2 readtable('AcclimationPilotData.csv').ATTMUp3]);
Pilot.ATTM.Info = {'ATTM = Abrupt transition temperature measurement', 'Columns correspond to Time, Down1, Down2, Down3, Up1, Up2, Up3 where down means temperature from bench to fridge and up means temperature from fridge to bench'};
%Feeding Intervals
Pilot.FeedingInterval4C = rmmissing(readtable('AcclimationPilotData.csv').FeedingInterval4C);
Pilot.FeedingInterval22C = rmmissing(readtable('AcclimationPilotData.csv').FeedingInterval22C);

%%% LAKE TEMPERATURE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('AcclimationLakeMichiganData.mat');
% LakeMichigan.temperature = ncread('glerl_southern_lake_michigan_temperature_mooring_2016_2017.nc', 'sea_water_temperature');
% LakeMichigan.longitude = ncread('glerl_southern_lake_michigan_temperature_mooring_2016_2017.nc', 'lon');
% LakeMichigan.latitude = ncread('glerl_southern_lake_michigan_temperature_mooring_2016_2017.nc', 'lat');
% LakeMichigan.z = ncread('glerl_southern_lake_michigan_temperature_mooring_2016_2017.nc', 'z');
% LakeMichigan.time_coverage_start = ncreadatt('glerl_southern_lake_michigan_temperature_mooring_2016_2017.nc','/','time_coverage_start');
% LakeMichigan.time_coverage_end = ncreadatt('glerl_southern_lake_michigan_temperature_mooring_2016_2017.nc','/','time_coverage_end');

%%% BEHAVIOR AND CA IMAGING DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('BehaviorData.mat');
load('CaImagingData.mat');
%%% COLORMAPS DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('cmaplist.mat')
%%% SAVE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save ('data.mat')
disp('Done');