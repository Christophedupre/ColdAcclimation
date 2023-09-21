function ModelData = GenerateEmptyModelStructure()

ModelData.tspan = [linspace(0,14,15) linspace(15,43,29)]; %Time vectors for training and testing
ModelData.y0 = [1 1]; %Initial conditions for acclimation and health

ModelData.AlphaFull = zeros(length(ModelData.tspan),15,16); %Stores Acclimation level from full simulation
ModelData.AlphaTruncated = zeros(29,15,16); %Stores Acclimation level for same time points as experimental data
ModelData.HealthFull = zeros(length(ModelData.tspan),15,16); %Stores Health level from full simulation
ModelData.HealthTruncated = zeros(29,15,16); %Stores Health level for same time points as experimental data

ModelData.NoVariability.AlphaFull = zeros(length(ModelData.tspan),15,16); %Stores Acclimation level from full simulation
ModelData.NoVariability.AlphaTruncated = zeros(29,15,16); %Stores Acclimation level for same time points as experimental data
ModelData.NoVariability.HealthFull = zeros(length(ModelData.tspan),15,16); %Stores Health level from full simulation
ModelData.NoVariability.HealthTruncated = zeros(29,15,16); %Stores Health level for same time points as experimental data

ModelData.NoiseAmplitude = 1.5; %Coefficient of variability. The smaller the coefficient the larger the variability
ModelData.TimeResolution = 24;
end