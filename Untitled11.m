clear all

K = 16;     % Sub-carriers                                             &此檔為LMMSE調變的IDFHE
M = 32;      % Sub-Symbols                                              
N = K*M ;  % N個samples                                                 
N_GFDM = 10; % number of GFDM blocks                                    
M_mod = 4;    % M-QAM                                                  
mm = log2(M_mod); % Per Symbol of bits                                  
N_bits = N_GFDM*N*mm ; % number of total bits                       
Ncp = N/16;
P = 2*(M_mod-1)/3 ; % M-QAM 的Signal Power
 
% 生成Ａ矩陣
%  Raised Cosine filter (time domain)
a = 0.1;  % rollof factor
t = linspace(-M/2, M/2, M*K+1); % 從-M/2到M/2，分成M*K+1個值 
t = t(1:end-1); % 減去向量的最後一列
t = t';
g = (sinc(t) .* cos(pi*a*t) ./ (1-4*a*a*t.*t)); % Raise cosine filter的脈衝響應

g = fftshift(g); % Y = fftshift(X) 重新排列的輸出fft，fft2和fftn通過零的頻率成分移動到陣列的中心。這對於可視化頻譜中間具有零頻率分量的傅立葉變換很有用。%對於向量，fftshift(X)交換的左右兩半X。對於矩陣，fftshift(X)將第一個像限與第三個像限交換，將第二個像限與第四個像限交換。
g(K+1:K:end) = 0; % 從K+1開始每隔K列到最後end為零
g = g / sqrt(sum(g.*g));
A = zeros(N, N);
for k=0:K-1
    for m=0:M-1
        for n=0:N-1
        A(:,m*K+k+1) = circshift(g , mod(n-m*K,N)) .* exp(1j*2*pi*(k/K)*[0:(N-1)])' ; % 最右邊的『　[0:(N-1)].'　』，代表每列每行的小儲存格都分別乘以它，才不會讓整列都乘以相同的數字
       end
    end
end



% (Circulant Channel Matrix)
r1 = zeros(1, N);
H = zeros(N, N);  % GFDM通道

%%%%%%%%%%%%% ＢＥＲ　%%%%%%%%%%%%%%%

SNR_range = 10:2:40;
BER_Er =  zeros(1,length(SNR_range));
BER_Er2 =  zeros(1,length(SNR_range));
BER_Er3 =  zeros(1,length(SNR_range));
BER_Er4 =  zeros(1,length(SNR_range));
BER_Er5 =  zeros(1,length(SNR_range));
last_BER_Er = BER_Er;
last_BER_Er2 = BER_Er2;
last_BER_Er3 = BER_Er3;
last_BER_Er4 = BER_Er4;
last_BER_Er5 = BER_Er5;

Tau = Ncp/5;
PDP = exp(-((0:Ncp-1)/Tau)).' ;
PDP = PDP/sum(PDP);
PDP2 = diag(sqrt(PDP));
N_ch = 20000;
h_all = PDP2*((randn(Ncp,N_ch)+1i*randn(Ncp,N_ch))/sqrt(2));
 

for abc = 1:N_ch % 通道跑的次數
    h = h_all(:,abc);
    r1(1:length(h)) = h;            
    for x = 1:N
        H(:, x) = circshift(r1,x-1);  %Ｈ矩陣的第ｒ行，循環移位 
    end      
    noise = (randn(N,N_GFDM) + 1i*randn(N,N_GFDM))/sqrt(2); % 除以根號２是為了讓noise實部跟虛部的power總和為１
     
    % 調製ＱＡＭ　Ｓｙｍｂｏｌ
    dataIn = randi([0 1],N_bits/mm,mm); % 生成二進制的 Data Sequence
    dataSymbol = bi2de(dataIn); % 將二進制的 Data Sequence 轉換成整數信號
    data_qam = (qammod(dataSymbol,M_mod))/sqrt(P); %除以sqrt(P)，正規化Signal Power
    data_qam = reshape(data_qam, N, N_GFDM); %最終D產生
       
    for  SNR_bit_dB = SNR_range
        nErrors_Er1 = 0 ;
        nErrors_Er2 = 0 ;
        nErrors_Er3 = 0 ;
        nErrors_Er4 = 0 ;         
        nErrors_Er5 = 0 ;
        
        SNR_bit = 10^(SNR_bit_dB/10);
        SNR_symbol = SNR_bit*mm ;      
        %%%%%%%%1st iteration%%%%%%%%
        H_eq = H*A;
        r = H_eq*sqrt(SNR_symbol)*data_qam+noise;%接收端收到的訊號
%         E = H_eq'*(H_eq*H_eq'+SNR_symbol^(-1)*eye(N))^(-1);
        E = (H_eq'*H_eq+SNR_symbol^(-1)*eye(N))^(-1)*H_eq';%LMMSE解調
        Er = E*r;
        Sig_Er = qamdemod(Er,4);% 做QPSK解調
        data_est_1 = (qammod(Sig_Er,M_mod))/sqrt(P).*(sqrt(SNR_symbol));%下一次疊代要扣除的ISI,先算好
        Sig_Er = reshape(Sig_Er,N*N_GFDM,1);
        dataOut_Er = de2bi(Sig_Er,mm); 
        nErrors_Er = biterr(dataIn,dataOut_Er); 
        BER_Er(SNR_bit_dB/2-SNR_range(1)/2+1) = nErrors_Er/(N_bits);
      
        %%%%%%%%2st iteration%%%%%%%%
        B = E*H_eq;
        BD = diag(diag(B));
        BI = B-BD;
      
        Er2 = Er-BI*data_est_1;
        Sig_Er2 = qamdemod(Er2,4);% 做QPSK解調
        data_est_2 = (qammod(Sig_Er2,M_mod))/sqrt(P).*(sqrt(SNR_symbol));%下一次疊代要扣除的ISI,先算好
        Sig_Er2 = reshape(Sig_Er2,N*N_GFDM,1);
        dataOut_Er2 = de2bi(Sig_Er2,mm); 
        nErrors_Er2 = biterr(dataIn,dataOut_Er2); 
        BER_Er2(SNR_bit_dB/2-SNR_range(1)/2+1) = nErrors_Er2/(N_bits);
    
%         %%%%%%%%%%3st iteration%%%%%%%%%
        Er3 = Er-BI*data_est_2;
        Sig_Er3 = qamdemod(Er3,4);
        data_est_3 = (qammod(Sig_Er3,M_mod))/sqrt(P).*(sqrt(SNR_symbol));%下一次疊代要扣除的ISI,先算好
        Sig_Er3 = reshape(Sig_Er3,N*N_GFDM,1);
        dataOut_Er3 = de2bi(Sig_Er3,mm); 
        nErrors_Er3 = biterr(dataIn,dataOut_Er3); 
        BER_Er3(SNR_bit_dB/2-SNR_range(1)/2+1) = nErrors_Er3/(N_bits) ;
     
       %%%%%%%%%E_CLS%%%%%%%%%%%%%%%%%
        [U,S,V] =svd(H_eq);%奇異值分解
        Ecls = V*U';
        B2 = Ecls*H_eq;
        BD2 = diag(diag(B2));%求出對角矩陣
        BI2 = B2-BD2;%求出非對角的部分
        Eclsr = Ecls*r;
     
        Er4 = Eclsr-BI2*data_est_2;
        Sig_Er4 = qamdemod(Er4,4);% 做QPSK解調
        data_est_4 = (qammod(Sig_Er4,M_mod))/sqrt(P).*(sqrt(SNR_symbol));%下一次疊代要扣除的ISI,先算好
        Sig_Er4 = reshape(Sig_Er4,N*N_GFDM,1);
        dataOut_Er4 = de2bi(Sig_Er4,mm); 
        nErrors_Er4 = biterr(dataIn,dataOut_Er4); 
        BER_Er4(SNR_bit_dB/2-SNR_range(1)/2+1) = nErrors_Er4/(N_bits);
        
        %%%%%%%%%E_CLS_2%%%%%%%%%%%%%%%%%
         Er5 = Eclsr-BI2*data_est_4;
         Sig_Er5 = qamdemod(Er5,4);% 做QPSK解調
         Sig_Er5 = reshape(Sig_Er5,N*N_GFDM,1);
         dataOut_Er5 = de2bi(Sig_Er5,mm); 
         nErrors_Er5 = biterr(dataIn,dataOut_Er5); 
         BER_Er5(SNR_bit_dB/2-SNR_range(1)/2+1) = nErrors_Er5/(N_bits);
     end
       last_BER_Er =  BER_Er+last_BER_Er;
       last_BER_Er2 =  BER_Er2+last_BER_Er2;
       last_BER_Er3 =  BER_Er3+last_BER_Er3;
       last_BER_Er4 =  BER_Er4+last_BER_Er4;
       last_BER_Er5 =  BER_Er5+last_BER_Er5;       
end      
BER_Er = last_BER_Er/N_ch ;
BER_Er2 = last_BER_Er2/N_ch ;
BER_Er3 = last_BER_Er3/N_ch ;
BER_Er4 = last_BER_Er4/N_ch ;
BER_Er5 = last_BER_Er5/N_ch ;
semilogy(SNR_range,BER_Er,'-o') 
xlabel('E_b/N_0(dB)','FontSize',15)
ylabel('BER','FontSize',15)
set(gca,'FontSize',12)
hold on
grid on
semilogy(SNR_range,BER_Er2,'-s') 
semilogy(SNR_range,BER_Er3,'-^') 
semilogy(SNR_range,BER_Er4,'-x') 
semilogy(SNR_range,BER_Er5,'-*') 
legend('non-iterative LMMSE','iterative LMMSE (I=2 iteration)','iteration LMMSE(I=3 iteration)','LMMSE-based IDFHE (I=3 iteration)','LMMSE-based IDFHE (I=4 iteration)')
ylim([1e-04 1e-1]) % 1e為10 10e為100
xlim([10 30])