v = 15.0;                    % UE velocity in km/h
fc = 4e9;                    % carrier frequency in Hz
c = physconst('lightspeed'); % speed of light in m/s
fd = (v*1000/3600)/c*fc;     % UE max Doppler frequency in Hz
 

cdl = nrCDLChannel;
cdl.DelayProfile = 'CDL-B';
cdl.MaximumDopplerShift = 300.0;
cdl.SampleRate = 10e3;
cdl.Seed = 19;

cdl.TransmitAntennaArray.Size = [1 1 1 1 1];
cdl.ReceiveAntennaArray.Size = [1 1 1 1 1];

T = 40; 
in = ones(T,1);

s = [Inf 5 2]; % sample densities

legends = {};
figure; hold on;
SR = cdl.SampleRate;
for i = 1:length(s)
    
    % call channel with chosen sample density
    release(cdl); cdl.SampleDensity = s(i);
    [out,pathgains,sampletimes] = cdl(in);
    chInfo = info(cdl); tau = chInfo.ChannelFilterDelay;
    
    % plot channel output against time
    t = cdl.InitialTime + ((0:(T-1)) - tau).' / SR;
    h = plot(t,abs(out),'o-'); 
    h.MarkerSize = 2; 
    h.LineWidth = 1.5;
    desc = ['Sample Density = ' num2str(s(i))];
    legends = [legends ['Output, ' desc]];
    disp([desc ', Ncs = ' num2str(length(sampletimes))]);
    
    % plot path gains against sample times
    h2 = plot(sampletimes-tau/SR,abs(sum(pathgains,2)),'o');
    h2.Color = h.Color; 
    h2.MarkerFaceColor = h.Color;
    legends = [legends ['Path Gains, ' desc]];    
end

xlabel('Time (s)');
title('Channel Output and Path Gains vs. Sample Density');
ylabel('Channel Magnitude');
legend(legends,'Location','NorthWest');