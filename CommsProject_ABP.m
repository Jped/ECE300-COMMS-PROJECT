%
% Communication Theory Projects 1 & 2
% Group: Shifra, Jonny, & Guy
%
% Main function
 
% A skeleton BER script for a wireless link simulation
clear all;close all;clc
% For the final version of this project, you must use these 3
% parameter. You will likely want to set numIter to 1 while you debug your
% link, and then increase it to get an average BER.
numIter = 1;  % The number of iterations of the simulation
nSym = 1000;    % The number of symbols per packet
SNR_Vec = 0:2:16;
lenSNR = length(SNR_Vec);

M = 2;   % The M-ary number, 2 corresponds to binary modulation.
% NOTE: M doesn't currenlty work for size 8 (and 32, and onwards,) due to
% the bits size being a multiple of M and then us deviding by a multiple of
% log2(M) when converting from bits to message (in the M=8 case, the bits
% will be 80000 length, yet we will be trying to devide this by log2(8)=3.
% We can solve this later on by padding the message as necessary, but for
% now I think we should just try and get M=4 and 16 working..

modulation = 2; % choose between PAM, QAM, and PSK (1,2,3 respectively)

%chan = 1;          % No channel
%chan = [1 .2 .4]; % Somewhat invertible channel impulse response, Moderate ISI
chan = [0.227 0.460 0.688 0.460 0.227]';   % Not so invertible, severe ISI


% Time-varying Rayleigh multipath channel, try it if you dare. Or take
% wireless comms.
% ts = 1/1000;
% chan = rayleighchan(ts,1);
% chan.pathDelays = [0 ts 2*ts];
% chan.AvgPathGaindB = [0 5 10];
% chan.StoreHistory = 1; % Uncomment if you want to be able to do plot(chan)

% Choose number of training symbols (max=len(msg)=1000)
numTrain = 200;

% Create a vector to store the BER computed during each iteration
berVec = zeros(numIter, lenSNR);

% Run the simulation numIter amount of times
for i = 1:numIter
    
    %  bits = randint(1, nSym*M, [0 1]);     % Generate random bits
    bits = randi(2,[nSym*log2(M), 1])-1;
    % New bits must be generated at every iteration
    
    % this function is correct for all cases (in particular, also M=2)
    msg = bits2msg(bits, M);

    
    for j = 1:lenSNR % one iteration of the simulation at each SNR Value
            
        % modulation
        if isequal(modulation, 1)
            tx = pammod(msg, M);  % PAM modulation
        elseif isequal(modulation, 2)
            tx = qammod(msg, M);  % QAM modulation
        else
            tx = pskmod(msg, M);  % PSK modulation
        end
        
       % Sequence of Training Symbols
        trainseq = tx(1:numTrain);
        % transmit (convolve) through channel
        if isequal(chan,1)
            txChan = tx;
        elseif isa(chan,'channel.rayleigh')
            reset(chan) % Draw a different channel each iteration
            txChan = filter(chan,tx);
        else
            txChan = filter(chan,1,tx);  % Apply the channel.
        end
        
        % Convert from EbNo to SNR.
        % Note: Because No = 2*noiseVariance^2, we must add ~3 dB to get SNR (because 10*log10(2) ~= 3).
        noise_addition = round(10*log10(2*log2(M)));
        %txNoisy = awgn(txChan, noise_addition+SNR_Vec(j), 'measured'); % Add AWGN
        txNoisy = txChan;
        
        % Equalizer
        stepsize = 0.01;
        forgetfactor = 0.3; %between 0 and 1
        
        % Adaptive Algorithm
        %AdapAlgo = lms(stepsize);
        AdapAlgo = rls(forgetfactor);
        
        % Number of Weights in Adaptive Filter
        %NWeights = length(chan); %Currently matched to size of channel. Maybe works better if more weights/larger filter order?
        NWeights = 10;
        NWEIGHTS_Feedback = 2;
        
        % Equalizer Object
        eqobj = lineareq(NWeights,AdapAlgo); %comparable to an FIR
        %eqobj = dfe(NWeights,NWEIGHTS_Feedback,AdapAlgo); %comparable to an IIR
        eqobj.SigConst = qammod((0:M-1)',M)';
        
        [yd,rx2] = equalize(eqobj,txNoisy,trainseq);
        
        % de-modulation
        if isequal(modulation, 1)
            rx = pamdemod(txNoisy, M);  % PAM
            rx2 = pamdemod(rx2, M);  % PAM
        elseif isequal(modulation, 2)
            rx = qamdemod(txNoisy, M);  % QAM
             rx2 = qamdemod(rx2, M);  % QAM
        else
            rx = pskdemod(txNoisy, M);  % PSK
            rx2 = pskdemod(rx2, M);  % PSK
        end
        rxMSG = msg2bits(rx, M);
        rxMSGE = msg2bits(rx2, M);
         %{       
        % Equalizer
        stepsize = .7;
        forgetfactor = 0.7; %between 0 and 1
        
        % Adaptive Algorithm
        AdapAlgo = lms(stepsize);
        %AdapAlgo = rls(forgetfactor);
        
        % Number of Weights in Adaptive Filter
        %NWeights = length(chan); %Currently matched to size of channel. Maybe works better if more weights/larger filter order?
        NWeights = 10;
        NWEIGHTS_Feedback = 2;
        
        % Equalizer Object
        %eqobj = lineareq(NWeights,AdapAlgo); %comparable to an FIR
        eqobj = dfe(NWeights,NWEIGHTS_Feedback,AdapAlgo); %comparable to an IIR
        eqobj.SigConst = qammod((0:M-1)',M)';
        
        rxE = equalize(eqobj,rx,trainseq);
        
        % converting the recieved, equalized message back into bits
        rxE(rxE<0) = 0;
        rxE(rxE>M-1) = M-1;
        
        rxE = round(rxE);
        rxE = rxE(numTrain+1:end); %no train
        rxMSGE = msg2bits(rxE, M);
        %}
        % Compute and store the EQUALIZER BER for this iteration
        % We're interested in the BER, which is the 2nd output of BITERR
        
        %bits_notrain = bits((numTrain)+1:end);
        rxMSGE_notrain = rxMSGE(numTrain+1:end);
        bits_notrain = bits(numTrain+1:end);
        [zzz, berVecENT(i,j)] = biterr(bits_notrain, rxMSGE_notrain); 
        
        
        % Compute and store the BER for this iteration
        % We're interested in the BER, which is the 2nd output of BITERR
        [zzz, berVecE(i,j)] = biterr(bits, rxMSGE); 
        [zzz, berVec(i,j)] = biterr(bits, rxMSG); 
        
    end  % End SNR iteration
end      % End numIter iteration


% Compute and plot the mean EQUALIZER BER
berE = mean(berVecE,1);

% Compute and plot the mean BER
ber = mean(berVec,1);
berNT = mean(berVecENT,1);

figure
semilogy(SNR_Vec, berE, 'DisplayName', 'Equalize')

hold on
semilogy(SNR_Vec, ber, 'DisplayName', 'BER')
hold on
semilogy(SNR_Vec, berNT, 'DisplayName', 'No TrainBER')

% Compute the theoretical BER for this scenario
% NOTE: there is no theoretical BER when you have a multipath channel
if isequal(modulation, 1) || (M<4) % if M<4, qam berawgn is anyways pam berawng
    berTheory = berawgn(SNR_Vec, 'pam', M); % PAM
elseif isequal(modulation, 2)
    berTheory = berawgn(SNR_Vec, 'qam', M); % QAM
else
    berTheory = berawgn(SNR_Vec, 'psk', M); % PSK
end

hold on
semilogy(SNR_Vec, berTheory, 'r', 'DisplayName', 'Theoretical BER')
xlabel('E_b/N_0(dB)');  ylabel('BER');
legend
grid


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [msg] = bits2msg(bits, M)
    % Convert the message from bits into the correct integer values
    % based on the inputted M-ary modulation.
    % NOTE: M has to be a multiple of 2.
    
    % The length of bits that will be converted into decimal.
    len = log2(M); 
    
    msg = zeros(size(bits,1)/len, 1);
    
    for i = 1:size(bits,1)/len
        msg(i) = bi2de(bits(1+(i-1)*len : 1+(i-1)*len + (len-1))');
    end
    
end


function [bits] = msg2bits(msg, M)
    % Convert the message from integers into the bit values
    % based on the inputted M-ary modulation.
    % NOTE: M has to be a multiple of 2.
    
    % The length of bits that will be converted into decimal.
    len = log2(M); 
    
    bits = zeros(1, size(msg,1)*len);
    
    for i = 1:size(msg,1)
        bits(1+(i-1)*len:1:1+(i-1)*len + (len-1)) = de2bi(msg(i), len);
    end
    bits = bits';
end
