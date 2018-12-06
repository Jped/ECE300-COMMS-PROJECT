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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HYPERPARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numIter = 10;  % The number of iterations of the simulation
nSym = 100000;    % The number of symbols per packet, worked well for 10,000 and 25 iterations
SNR_Vec = 0:2:16;
lenSNR = length(SNR_Vec);

% The M-ary number. 2 corresponds to binary modulation.
M = 2;  

% Modulation
%  - 1 = PAM
%  - 2 = QAM
%  - 3 = PSK
modulation = 2;

%chan = 1;          % No channel
chan = [1 .2 .4]; % Somewhat invertible channel impulse response, Moderate ISI
%chan = [0.227 0.460 0.688 0.460 0.227]';   % Not so invertible, severe ISI

% Time-varying Rayleigh multipath channel, try it if you dare. Or take wireless comms.
% ts = 1/1000;
% chan = rayleighchan(ts,1);
% chan.pathDelays = [0 ts 2*ts];
% chan.AvgPathGaindB = [0 5 10];
% chan.StoreHistory = 1; % Uncomment if you want to be able to do plot(chan)

% Choose number of training symbols (max=len(msg)=1000)
numTrain = 175;

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

%numRefTap = 2;

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
    %eqobj.RefTap = numRefTap;
else
    eqobj = dfe(NWeights, NWEIGHTS_Feedback, AdapAlgo); %comparable to an IIR
end

%%%%%%%%%%%%%%%%%%% CREATING ERROR CONTROL CODING SCHEME %%%%%%%%%%%%%%%%%%
% The information to be encoded consists of message symbols and the code 
% that is produced consists of codewords. Each block of K message symbols 
% is encoded into a codeword that consists of N message symbols. K is called 
% the message length, N is called the codeword length, and the code is 
% called an [N,K] code.

%You can structure messages and codewords as binary vector signals, where 
% each vector represents a message word or a codeword. At a given time, the 
% encoder receives an entire message word, encodes it, and outputs the 
% entire codeword. The message and code signals operate over the same sample time.

% BCH: For these codes, the codeword length N must have the form 2M-1, 
% where M is an integer from 3 to 16 (default is 15). The message length K is restricted to 
% particular values that depend on N. To see which values of K are valid 
% for a given N, see the comm.BCHEncoder System object? reference page. 
%No known analytic formula describes the relationship among the codeword 
%length, message length, and error-correction capability for BCH codes.
%Message length default is 5.

X = 4; % integer from 3 to 16; the documentation uses the variable M in 
% place of x, but this is confusing because this value is different than
% the modulation value M.

%CODEWORD LENGTH:
n = 2^X-1; % default is 15, max allowed is 65,535
% This function: bchnumerr(n) will return all possible K/message values for
% a particular N and the number of correctable errors in a three column
% matrix

%MESSAGE LENGTH:
k = 5; %default is 5, Example: 5 specifies a Galois array with five elements, 2^m(second value in gf)
paritypos='beginning';
%paritypos = 'end' or 'beginning' specify whether parity bits at end or
%beginning of signal

% The message must be fed into the encoder using the Galois field array of
% symbols over GF(2) Each K-element row of msg represents a message word, 
% where the leftmost symbol is the most significant symbol.
%%%%%%%%%%% THIS IS LATER ON: msgTx = gf(x,1) % Create a Galois array in GF(2^m)
%The elements of x must be integers between 0 and 2^m-1, if only using bits
% (x={0,1}), than m=1/second arguement and we are in GF(2)

%%%%%%%%%%% THIS IS LATER ON:enc = bchenc(msgTx,n,k,paritypos) %This occurs before adding noise after 
%convert from message to bits
%THIS IS HOW NOISE ADDITION IS DONE IN DOC EX:
%Corrupt up to t bits in each codeword where t = bchnumerr(n,k)
%noisycode = enc + randerr(numbits,n,1:t). This is for full
%reconstruction with no errors and so doesnt concern us because we are
%adding AWGN.

%%%%%%%%%%% THIS IS LATER ON:msgRx = bchdec(noisycode,n,k) %decode the noisy message

% Order of Operations
% source:https://www.researchgate.net/figure/System-model_fig2_303940773
% Data Source --> Conversion to Bits --> BCH Encoder --> Modulation -->
% Filter/Channel --> Recieve Filter/AWGN --> Demodulation --> BCH Decoder
% --> Convert back to message

%%%%%%%%%%%%%%%%%%%%%%%%%%% RUNNING SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a vector to store the BER computed during each iteration
berVec = zeros(numIter, lenSNR);
berVecE = zeros(numIter, lenSNR);

for i = 1:numIter
    disp('start iteration')
     % message to transmit
    %bits = randi(2,[nSym*log2(M), 1])-1;
    %msg = bits2msg(bits, M);
    %msg= reshape(msg,[nSym/X,X]);
    %msg_gf = gf(msg,log2(M));
    %msg_RS = rsenc(msg_gf,n,X);
    %msg_RS_x = msg_RS.x;
    %msg_RS_x = double(msg_RS_x(:));
    
    % message to transmit
    bits = randi(2,[nSym*log2(M), 1])-1;
    msg = bits2msg(bits, M);
    msgreshape = reshape(msg,[nSym/k,k]);
    msg_gf = gf(msgreshape,log2(M));
    msg_enc = bchenc(msg_gf,n,k);
    msg_enc = msg_enc.x;
    
    for j = 1:lenSNR % one iteration of the simulation at each SNR Value
        % Not totally sure if encoding should occur inside or outside of
        % this loop
        % BCH encoding
        
        %Must first reshape msg so that each row has k elements
        %X = log2(M);
        %msgreshape = reshape(msg,k,[]).';
        %msg_gf = gf(msgreshape,1);
        %msg_enc = bchenc(msg_gf,n,k,paritypos);
        %msg_enc = msg_enc.x;
        
        %Now must unwrap matrix into vector to input into modulation fns
        %msg_enc = msg_enc_matrix.';
        %msg_enc = msg_enc(:);
        %ISSUE IS FEEDING Galois Field Array into modulation schemes
        
        % modulation
        if isequal(modulation, 1)
            tx = pammod(msg_enc, M);  % PAM modulation
        elseif isequal(modulation, 2)
            tx = qammod(msg_enc(:), M);  % QAM modulation
            txE = qammod(msg, M); %QAM nonencoded mod
        else
            tx = pskmod(msg_enc, M);  % PSK modulation
        end
        
        % Sequence of Training Symbols
        trainseq = tx(1:numTrain);
        trainseqE = txE(1:numTrain); %no encoding
        
        % transmit (convolve) through channel
        if isequal(chan,1)
            txChan = tx;
        elseif isa(chan,'channel.rayleigh')
            reset(chan) % Draw a different channel each iteration
            txChan = filter(chan,tx);
        else
            txChan = filter(chan,1,tx);  % Apply the channel.
            txChanE = filter(chan,1,txE);
        end
        
        % Convert from EbNo to SNR.
        % Note: Because No = 2*noiseVariance^2, we must add ~3 dB to get SNR (because 10*log10(2) ~= 3).
        noise_addition = round(10*log10(2*log2(M)));
        txNoisy = awgn(txChan, 3+SNR_Vec(j), 'measured'); % Add AWGN
        txNoisyE = awgn(txChanE, 3+SNR_Vec(j), 'measured'); %no coding
        
        yd = equalize(eqobj,txNoisy,trainseq);
        ydE = equalize(eqobj,txNoisyE,trainseqE);
        
        % de-modulation
        if isequal(modulation, 1)
            rx = pamdemod(yd, M);  % PAM
        elseif isequal(modulation, 2)
            rx = qamdemod(yd, M);  % QAM
            rxE = qamdemod(ydE,M);
        else
            rx = pskdemod(yd, M);  % PSK
        end
        
        % BCH decoder
        %{
        rx = gf(reshape(rx, [size(rx,1)/n,n]),msg_gf.m, msg_gf.prim_poly);
        [rx_decode, cnummerr] = rsdec(rx,n,X);
        rx_decode = rx_decode.x;
        rx_decode = double(rx_decode(:));
        rxMSG = msg2bits(rx2, M);
        rxMSG_de= msg2bits(rx_decode,M);
        %}
        
        rx = gf(reshape(rx,[size(rx,1)/n,n]),msg_gf.m,msg_gf.prim_poly);
        msgRx = bchdec(rx,n,k);
        msgRx = msgRx.x;
        
        rxMSG = msg2bits(msgRx(:), M);
        rxMSGE = msg2bits(rxE,M);
        
        % Compute and store the BER for this iteration
        % We're interested in the BER, which is the 2nd output of BITERR
        [~, berVecE(i,j)] = biterr(bits(numTrain:end), rxMSGE(numTrain:end)); %no coding
        [~, berVec(i,j)] = biterr(bits,rxMSG); %with coding
        
    end  % End SNR iteration
end      % End numIter iteration


% Compute and plot the mean EQUALIZER BER
berE = mean(berVecE,1); %no coding
ber2 = mean(berVec,1); %with coding
figure
semilogy(SNR_Vec, berE, 'DisplayName', 'Equalize no coding')
hold on
semilogy(SNR_Vec, ber2, 'DisplayName', 'Equalize with coding')

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
title('BER Rates of BPSK')
legend
grid


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%adaptives = 0.9:0.005:1;
%weights = 5:9;
%trains = 25:25:300;
%[best_adap_val, best_n_weight, best_train_len] = gridSearch(2, 0, adaptives, weights, trains);
%%

% for lms and lineareq, best vals were: adap_val=0.005, n_weight=8, best_train_len=75
% for rls and lineareq, best vals were: adap_val=1, n_weight=6, best_train_len=175
% for lms and def, best vals were: adap_val=, n_weight=, feed_back_weights =, best_train_len=
% for rls and def, best vals were: adap_val=, n_weight=, , feed_back_weights =, best_train_len=

function [optimal_adaptive_algo, optimal_n_weight, optimal_train_len] = gridSearch(adaptive_algo, equalizer, adaptive_algo_contenders, n_weight_contenders, num_train_contenders)
    % Performs a grid search across the specified parameters, returning
    % the optimal values find.
    
    M = 2;
    nSym=20000;
    modulation = 2;
    chan = [1 .2 .4];
    
    NWEIGHTS_Feedback = 5;
    optimal_adaptive_algo = -1;
    optimal_n_weight = -1;
    optimal_train_len = -1;
    best_score = 1;
    for cur_adaptive_algo = adaptive_algo_contenders
        for cur_n_weight = n_weight_contenders
            for cur_train_val = num_train_contenders
                %fprintf('Current run: adaptive: %f, weight: %d, train: %d\n', cur_adaptive_algo, cur_n_weight, cur_train_val);
                
                % adaptive filter algorithm
                if isequal(adaptive_algo, 0)
                    AdapAlgo = varlms(cur_adaptive_algo,0.01,0,0.01);
                elseif isequal(adaptive_algo, 1)
                    AdapAlgo = lms(cur_adaptive_algo);
                else
                    AdapAlgo = rls(cur_adaptive_algo);
                end

                % Equalizer Object
                if isequal(equalizer, 0)
                    eqobj = lineareq(cur_n_weight,AdapAlgo); %comparable to an FIR
                else
                    eqobj = dfe(cur_n_weight,NWEIGHTS_Feedback,AdapAlgo); %comparable to an IIR
                end

                % message to transmit
                bits = randi(2,[nSym*log2(M), 1])-1;
                msg = bits2msg(bits, M);

                % modulation
                if isequal(modulation, 1)
                    tx = pammod(msg, M);  % PAM modulation
                elseif isequal(modulation, 2)
                    tx = qammod(msg, M);  % QAM modulation
                else
                    tx = pskmod(msg, M);  % PSK modulation
                end

                % Sequence of Training Symbols
                trainseq = tx(1:cur_train_val);

                % transmit (convolve) through channel
                txChan = filter(chan,1,tx);  % Apply the channel.

                % Convert from EbNo to SNR.
                % we are specifically intersted in the SNR value of 12dB.
                txNoisy = awgn(txChan, 3+12, 'measured'); % Add AWGN     

                yd = equalize(eqobj,txNoisy,trainseq);

                % de-modulation
                if isequal(modulation, 1)
                    rx = pamdemod(yd, M);  % PAM

                elseif isequal(modulation, 2)
                    rx = qamdemod(yd, M);  % QAM
                else
                    rx = pskdemod(yd, M);  % PSK

                end
                rxMSGE = msg2bits(rx, M);

                % Compute and store the BER for this iteration
                % We're interested in the BER, which is the 2nd output of BITERR
                [~, ber] = biterr(bits(cur_train_val:end), rxMSGE(cur_train_val:end)); 
                
                if ber ~= 0 && ber < best_score
                    best_score = ber;
                    optimal_train_len = cur_train_val;
                    optimal_adaptive_algo = cur_adaptive_algo;
                    optimal_n_weight = cur_n_weight;
                    fprintf('New best ber achieved: %e.\n', ber);
                end
            end
        end
    end

end


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
