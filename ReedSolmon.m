clear all;
close all;
clc

numIter = 125;  % The number of iterations of the simulation
nSym = 10003;    % The number of symbols per packet
SNR_Vec = 12;
lenSNR = length(SNR_Vec);

M = 8;  

% Modulation
%  - 1 = PAM
%  - 2 = QAM
%  - 3 = PSK
modulation = 3;

%chan = 1;          % No channel
chan = [1 .2 .4]; % Somewhat invertible channel impulse response, Moderate ISI
%chan = [0.227 0.460 0.688 0.460 0.227]';   % Not so invertible, severe ISI
numTrain = 100;

% Adaptive Algorithm
%  - 0 = varlms
%  - 1 = lms
%  - 2 = rls
adaptive_algo = 2;

% Equalizer
%  - 0 = lineareq
%  - 1 = dfe
equalize_val = 0;

% equalizer hyperparameters
NWeights = 6;
NWEIGHTS_Feedback = 5;

numRefTap = 3;

stepsize = 0.005;
forgetfactor = 1; % between 0 and 1


%%%%%%%%%%%%%%%%%%%%%%%%% CREATING EQUALIZER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% adaptive filter algorithm
if isequal(adaptive_algo, 0)
    AdapAlgo = varlms(stepsize,0.01,0,0.01);
elseif isequal(adaptive_algo, 1)
    AdapAlgo = lms(stepsize);
else
    AdapAlgo = rls(forgetfactor);
end

% Equalizer Object
if isequal(equalize_val, 0)
    eqobj = lineareq(NWeights,AdapAlgo); %comparable to an FIR
    eqobj.RefTap = numRefTap;
    delay = (numRefTap - 1)/eqobj.nSampPerSym;
    %eqobj.ResetBeforeFiltering = 0;
else
    eqobj = dfe(NWeights, NWEIGHTS_Feedback, AdapAlgo); %comparable to an IIR
end

X = log2(M);

%CODEWORD LENGTH:
n = 2^X-1; % default is 15, max allowed is 65,535
k =5; %default is 5, Example: 5 specifies a Galois array with five elements, 2^m(second value in gf)
paritypos='end';


% Create a vector to store the BER computed during each iteration
berVecE = zeros(numIter, lenSNR);
berVec2 = zeros(numIter, lenSNR);

for i = 1:numIter
    
    % message to transmit
    bits = randi(2,[(nSym + delay)*log2(M), 1])-1;
    msg = bits2msg(bits, M);
    msg= reshape(msg,[(nSym + delay)/X,X]);
    msg_gf = gf(msg,log2(M));
    msg_RS = rsenc(msg_gf,n,X);
    msg_RS_x = msg_RS.x;
    msg_RS_x = double(msg_RS_x(:));
    % modulation
    if isequal(modulation, 1)
        tx = pammod(msg, M);  % PAM modulation
    elseif isequal(modulation, 2)
        tx = qammod(msg_RS_x, M);% QAM modulation
        tx2 = qammod(msg(:),M);
    else
        tx = pskmod(msg_RS_x, M);  % PSK modulation
        tx2 = pskmod(msg(:),M);
    end
    trainseq = tx(1:numTrain);
    trainseq2= tx2(1:numTrain);
    % transmit (convolve) through channel
    if isequal(chan,1)
        txChan = tx;
        txChan2 = tx2;
    elseif isa(chan,'channel.rayleigh')
        reset(chan) % Draw a different channel each iteration
        txChan = filter(chan,tx);
    else
         txChan = filter(chan,1,tx);  % Apply the channel.
         txChan2 = filter(chan,1,tx2);
    end

    % Convert from EbNo to SNR.
    % Note: Because No = 2*noiseVariance^2, we must add ~3 dB to get SNR (because 10*log10(2) ~= 3).
    noise_addition = 10*log10(log2(M));
    for j = 1:lenSNR % one iteration of the simulation at each SNR Value
        
      
        txNoisy = awgn(txChan, noise_addition+SNR_Vec(j), 'measured'); % Add AWGN
        txNoisy2 = awgn(txChan2,noise_addition+SNR_Vec(j), 'measured');
        yd = equalize(eqobj,txNoisy,trainseq);
        yd2 = equalize(eqobj,txNoisy2,trainseq2);
        % de-modulation
        if isequal(modulation, 1)
            rx = pamdemod(yd, M);  % PAM
        elseif isequal(modulation, 2)
            rx = qamdemod(yd, M);  % QAM
            rx2 = qamdemod(yd2,M);
        else
            rx = pskdemod(yd, M);  % PSK
            rx2 = pskdemod(yd2,M);
        end
        rx = gf(reshape(rx, [size(rx,1)/n,n]),msg_gf.m, msg_gf.prim_poly);
        [rx_decode, cnummerr] = rsdec(rx,n,X);
        rx_decode = rx_decode.x;
        rx_decode = double(rx_decode(:));
        rxMSG = msg2bits(rx2, M);
        rxMSG_de= msg2bits(rx_decode,M);
        % Compute and store the BER for this 
        % We're interested in the BER, which is the 2nd output of BITERR
        numTrainBits = numTrain*log2(M);
        [~, berVecE(i,j)] = biterr(bits(numTrainBits:end-delay), rxMSG_de(numTrainBits+delay:end)); 
        [~, berVec2(i,j)] = biterr(bits(numTrainBits:end), rxMSG(numTrainBits:end)); 
    end  % End SNR iteration
end      % End numIter iteration


% Compute and plot the mean EQUALIZER BER
berE = mean(berVecE,1); %no coding
ber2 = mean(berVec2,1);
figure
semilogy(SNR_Vec, berE, 'DisplayName', 'Equalize with coding')
hold on
semilogy(SNR_Vec, ber2, 'DisplayName', 'Equalize with out coding')
% Compute the theoretical BER for this scenario
% NOTE: there is no theoretical BER when you have a multipath channel
if isequal(modulation, 1) || (M<4) % if M<4, qam berawgn is anyways pam berawng
    berTheory = berawgn(SNR_Vec, 'pam', M); % PAM
elseif isequal(modulation, 2)
    berTheory = berawgn(SNR_Vec, 'qam', M); % QAM
else
    berTheory = berawgn(SNR_Vec, 'psk', M, 'nondiff'); % PSK
end



function [msg] = bits2msg(bits, M)
    % Convert the message from bits into the correct integer values
    % based on the inputted M-ary modulation.
    % NOTE: M has to be a multiple of 2.
    
    % The length of bits that will be converted into decimal.
    len = log2(M)
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
        if msg(i) >= M
            msg(i) = 0;
        end
        bits(1+(i-1)*len:1:1+(i-1)*len + (len-1)) = de2bi(msg(i), len);
    end
    bits = bits';
end
