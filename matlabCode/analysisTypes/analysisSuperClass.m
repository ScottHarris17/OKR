classdef analysisSuperClass
    properties
        data %structure containing cell's data
        spikeDetector
        recordingType
        analysisType
        cellID
        cellType
        epochsSelected
        analysisFigures = gobjects(0)%gobjects array
        analysisFigureExtensions = {};%file name extension for each figure in analysisFigures
        processedData
    end
    
    methods (Hidden = true)
                %think about how to put a series of functions for raw
                %traces here
    end
    
    
    methods %general methods that can be performed on traces independent of their analysis by a speific analysis type
        
        function obj = lowpassfilter(obj, threshold)
            %low pass filter each epoch
            switch nargin
                case 1
                    threshold = 500; %hz
            end
            
            for i = 1:numel(obj.epochsSelected)
                obj.data.(obj.cellID).epochs(obj.epochsSelected(i)).epoch =...
                    lowpass(obj.data.(obj.cellID).epochs(obj.epochsSelected(i)).epoch, threshold,...
                    obj.data.(obj.cellID).epochs(obj.epochsSelected(i)).meta.sampleRate);
            end   
        end
        
        
        function obj = highpassfilter(obj, threshold)
            %high pass filter each epoch
            switch nargin
                case 1
                    threshold = 50; %hz
            end
            
            for i = 1:numel(obj.epochsSelected)
                obj.data.(obj.cellID).epochs(obj.epochsSelected(i)).epoch =...
                    highpass(obj.data.(obj.cellID).epochs(obj.epochsSelected(i)).epoch, threshold,...
                    obj.data.(obj.cellID).epochs(obj.epochsSelected(i)).meta.sampleRate);
            end   
        end
        
        
        %Build a psth
        function psth = buildPSTH(obj, spikes, binSize, startTime, endTime)
            numBins = ceil((endTime - startTime)/binSize);

            fullMatr = zeros(numel(spikes), numBins);

            for i = 1:numel(spikes)
                spikes_i = spikes{i};
                for j = 1:numel(spikes_i)
                    binNum = floor((spikes_i(j)-startTime)/binSize)+1; if binNum > numBins; binNum = numBins; end
                    fullMatr(i, binNum) = fullMatr(i, binNum) + 1;
                end
            end


            psth = mean(fullMatr, 1);
        end

        %rb
        function [freqs, power] = buildPowerSpectrum(obj, timeseries, duration, stepsize)
            %Adpated from "An Introduction to Field Analysis Techniques: The Power Spectrum and Coherence" by Mark A. Kramer, PhD
            xf = fft(timeseries); %1. Compute the Fourier transform.
            Sxx = 2*stepsize^2/duration * xf.*conj(xf); %2. Compute the power spectrum.
            power = Sxx(1:length(timeseries)/2+1); %3. Ignore negative frequencies.
            %powerDB = mag2db(power);%convert to decibels (logscale)
            df = 1/duration; %4. Determine the frequency resolution.
            fNQ=1/stepsize/2; %5. Determine the Nyquistfrequency.
            freqs = (0:df:fNQ); %6. Construct the frequency axis.

            %remove 0 frequency
            freqs = freqs(2:end);
            power = power(2:end);
        end
    end
end