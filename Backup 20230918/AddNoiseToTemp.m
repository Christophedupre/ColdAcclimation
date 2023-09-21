function temperature2 = AddNoiseToTemp(temperature2, NoiseAmplitude)
temperature2 = temperature2+((rand(1,length(temperature2))-0.5)*NoiseAmplitude);

end