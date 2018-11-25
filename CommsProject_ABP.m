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



%%%%%%%%%%%%%%%%%%% HYPERPARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numIter = 10;  % The number of iterations of the simulation
nSym = 20000;    % The number of symbols per packet
SNR_Vec = 0:2:16;
lenSNR = length(SNR_Vec);

M = 2;  

% The M-ary number, 2 corresponds to binary modulation.
% NOTE: M doesn't currenlty work for size 8 (and 32, and onwards,) due to
% the bits size being a multiple of M and then us deviding by a multiple of
% log2(M) when converting from bits to message (in the M=8 case, the bits
% will be 80000 length, yet we will be trying to devide this by log2(8)=3.
% We can solve this later on by padding the message as necessary, but for
% now I think we should just try and get M=4 and 16 working..

modulation = 2; % choose between PAM, QAM, and PSK (1,2,3 respectively)

%chan = 1;          % No channel
chan = [1 .2 .4]; % Somewhat invertible channel impulse response, Moderate ISI
%chan = [0.227 0.460 0.688 0.460 0.227]';   % Not so invertible, severe ISI


% Time-varying Rayleigh multipath channel, try it if you dare. Or take
% wireless comms.
% ts = 1/1000;
% chan = rayleighchan(ts,1);
% chan.pathDelays = [0 ts 2*ts];
% chan.AvgPathGaindB = [0 5 10];
% chan.StoreHistory = 1; % Uncomment if you want to be able to do plot(chan)

% Choose number of training symbols (max=len(msg)=1000)
numTrain = 100;

% Equalizer
stepsize = 0.01;
forgetfactor = 0.3; %between 0 and 1

% Number of Weights in Adaptive Filter
NWeights = length(chan) + 5;
NWEIGHTS_Feedback = 2;
% Adaptive Algorithm
AdapAlgo = lms(stepsize);
%AdapAlgo = rls(forgetfactor);

% Equalizer Object
eqobj = lineareq(NWeights,AdapAlgo); %comparable to an FIR
%eqobj = dfe(NWeights,NWEIGHTS_Feedback,AdapAlgo); %comparable to an IIR


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a vector to store the BER computed during each iteration
berVec = zeros(numIter, lenSNR);

for i = 1:numIter
    

    bits = randi(2,[nSym*log2(M), 1])-1;
   
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
        noise_addition = round(10*log10(log2(M)));
        txNoisy = awgn(txChan, noise_addition+SNR_Vec(j), 'measured'); % Add AWGN
        %txNoisy = txChan;
        
       
       
        
        
        
        [symbolest,yd] = equalize(eqobj,txNoisy,trainseq);
        
        % de-modulation
        if isequal(modulation, 1)
            rx = pamdemod(yd, M);  % PAM
     
        elseif isequal(modulation, 2)
            rx = qamdemod(yd, M);  % QAM
            rx2 = qamdemod(txNoisy,M);
        else
            rx = pskdemod(yd, M);  % PSK
            
        end
        rxMSGE = msg2bits(rx, M);
        rxMSG2 = msg2bits(rx2,M);
        % Compute and store the EQUALIZER BER for this iteration
        % We're interested in the BER, which is the 2nd output of BITERR
        
        
        
        % Compute and store the BER for this iteration
        % We're interested in the BER, which is the 2nd output of BITERR
        [zzz, berVecE(i,j)] = biterr(bits(numTrain:end), rxMSGE(numTrain:end)); 
        [zzz, berVec(i,j)] = biterr(bits,rxMSG2);
    end  % End SNR iteration
end      % End numIter iteration


% Compute and plot the mean EQUALIZER BER
berE = mean(berVecE,1);
ber2 = mean(berVec,1);
figure
semilogy(SNR_Vec, berE, 'DisplayName', 'Equalize')
hold on
semilogy(SNR_Vec, ber2, 'DisplayName', 'No equalize')

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
