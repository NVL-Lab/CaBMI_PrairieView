function plot_Neurons_Ensemble(baseActivity, ensembleNeurons, ensembleID)
%{
Function to plot the temporal activity of neurons collected during Baseline
to select the best neurons.
baseActivity -> activity during baseline
CComp -> C_on from holostim period given by onacid
YrA -> background noise of C
totalNeurons -> amount of neurons to be displayed

%}

    plot_b = [0    0.4470    0.7410];
    plot_o = [0.8500    0.3250    0.0980]; 
    E_color = {plot_b, plot_o}; 

    totalNeurons = length(ensembleNeurons);
    
    subplotnmb = ceil(totalNeurons/2);
	figure('Position', [300,300, subplotnmb*200, 400])
    %sgtitle('Std/mean')
    for idx = 1:length(ensembleNeurons)
		subplot(2,subplotnmb,idx)
        plot_color = E_color{ensembleID(idx)};
		plot(baseActivity(ensembleNeurons(idx), :)', 'color', plot_color);
		title(['ROI ' int2str(ensembleNeurons(idx))]);
        label = num2str(ensembleID(idx)); 
        legend(label); 
    end
end
    