%% Hierarchical bootstrapping
clear all
clc
NoOfsamplesInBootstrap = 2000;
close all;
figure
for z = 1:5
   
    if (z == 1)
        load('C:\Users\Bert\Documents\MATLAB\BootstrapData\cond1_SpikeRate.mat')
        currentCell = cond1_sr;
        condText = 'Cond 1 - ctrl';
    end;
     if (z == 2)
         load('C:\Users\Bert\Documents\MATLAB\BootstrapData\cond2_SpikeRate.mat')
        currentCell = cond2_sr;
        condText = 'Cond 2 - NMDA,AMPA';
    end;
     if (z == 3)
         load('C:\Users\Bert\Documents\MATLAB\BootstrapData\cond3_SpikeRate.mat')
        currentCell = cond3_sr;
         condText = 'Cond 3 - GABA,AMPA';
    end;
     if (z == 4)
         load('C:\Users\Bert\Documents\MATLAB\BootstrapData\cond5_SpikeRate.mat')
        currentCell = cond5_sr;
        condText = 'Cond 5 - AMPA';
    end;
     if (z == 5)
         load('C:\Users\Bert\Documents\MATLAB\BootstrapData\cond6_SpikeRate.mat')
        currentCell = cond6_sr;
        condText = 'Cond 6 - NMDA';
    end;
    
    for i = 1:NoOfsamplesInBootstrap %gathering 2000 samples
     FinalResampled_L2 = [];
     arrIndices_L1 = datasample(1:length(currentCell),length(currentCell),'Replace',true);
     for j = 1:length(arrIndices_L1)
         Origsamples_L2 = currentCell{arrIndices_L1(j)};
         Resampled_L2 =  Origsamples_L2(datasample(1:length(Origsamples_L2),length(Origsamples_L2),'Replace',true));
         FinalResampled_L2 = [FinalResampled_L2 Resampled_L2'];
     end;
     SampleMean(i) = mean(FinalResampled_L2);
    end;
    [fi,xi] = ksdensity(SampleMean);
%     figure
     subplot(2,3,z);
     plot(xi,fi);
     title(['Delta - ' condText ',sample size = ' num2str(NoOfsamplesInBootstrap)]);
    display('Mean +- 90_CI');
    display([num2str(xi(5)) ',' num2str(mean(SampleMean)) ',' num2str(xi(95))]);
    clear currentCell FinalResampled_L2 arrIndices_L1 Origsamples_L2 Resampled_L2 SampleMean;
    
end;

 
 

