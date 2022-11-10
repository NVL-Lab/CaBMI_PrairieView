function [vectorHolo, ISI] = create_Vector_random_stim(frameRate, expectedLengthExperiment, IHSImean, IHSIrange, toplot)
%Inputs: 
% expectedLengthExperiment (frames)


averageStimPeriod = IHSImean*frameRate;%In frames:
range = IHSIrange*frameRate; %ISI will vary between [averageStimPeriod-range averageStimPeriod+range]
num_stims = round(expectedLengthExperiment/averageStimPeriod);
ISI = rand(num_stims,1)*2*range +(averageStimPeriod-range); 
vectorHolo = round(cumsum(ISI));

if toplot
    % Debug:
    h = figure;
    plot(vectorHolo, '.-', 'MarkerSize', 15)
    ylabel('Frame Number'); 
    xlabel('Stim Number'); 
    title('Frame Number vs Stim Number'); 

    h = figure;
    plot(vectorHolo/(frameRate*60), '.-', 'MarkerSize', 15)
    ylabel('Time (min)'); 
    title('Stim time (min) vs Stim Number'); 

    h = figure;
    hist(ISI, 100); 
    xlabel('ISI (frames)'); 
    ylabel('count'); 
    title('ISI distribution (unnormalized)'); 
end
end