input = 115500; T = 8;
enc7 = fec.bchenc(15,7);
enc11 = fec.bchenc(15,11);
decoder7 = fec.bchdec(enc7);
decoder11 = fec.bchdec(enc11);
fin_brr_sq = [0,0,0,0,0,0,0,0];
fin_brr_t = [0,0,0,0,0,0,0,0];
fin_brr_sq7 = [0,0,0,0,0,0,0,0];
fin_brr_t7 = [0,0,0,0,0,0,0,0];
fin_brr_sq11 = [0,0,0,0,0,0,0,0];
fin_brr_t11 = [0,0,0,0,0,0,0,0];
snr=[-5 -2 0 2 5 8 12 15];
for f1=1:5
    msg = randi([0 1],input,1);
    bch7 = encode(enc7,msg);
    bch11  = encode(enc11,msg);
    triangle_weights = [0,1/4,1/2,3/4,1,3/4,1/2,1/4];
    onwa = [1,1,1,1,1,1,1,1];
    msg(msg==0)= -1;
    bch7(bch7==0) = -1;
    bch11(bch11==0) = -1;
    
    temp_msg7 = kron(bch7,ones(T,1));
    square_msg7 = reshape(temp_msg7,T,size(bch7,1));
    triangle_msg7 = diag(triangle_weights)*square_msg7;
    
    
    temp_msg11 = kron(bch11,ones(T,1));
    square_msg11 = reshape(temp_msg11,T,size(bch11,1));
    triangle_msg11 = diag(triangle_weights)*square_msg11;
    
    
    temp_msg = kron(msg,ones(T,1));
    square_msg = reshape(temp_msg,T,input);
    triangle_msg = diag(triangle_weights)*square_msg;
    
    brror_sq_snr=[];
    brror_t_snr=[];
    brror_sq_snr7=[];
    brror_t_snr7=[];
    brror_sq_snr11=[];
    brror_t_snr11=[];
    
    for f2=1:numel(snr)
        y_triangle = awgn(triangle_msg,snr(f2),'measured');
        y_square = awgn(square_msg,snr(f2),'measured');
        y_triangle7 = awgn(triangle_msg7,snr(f2),'measured');
        y_square7 = awgn(square_msg7,snr(f2),'measured');
        y_triangle11 = awgn(triangle_msg11,snr(f2),'measured');
        y_square11 = awgn(square_msg11,snr(f2),'measured');
        
        sze = size(y_square,2);
        sze7 = size(y_square7,2);
        sze11 = size(y_square11,2);
        
        
        
        sig=[];
        sig7=[];
        sig11=[];
        sig_t=[];
        sig7_t=[];
        sig11_t=[];
        
        
        
        
        
        %-----------------------------------------------------%
        for i=1:sze
            a = mean(xcorr(onwa,y_square(:,i)));
            b = mean(xcorr(triangle_weights,y_triangle(:,i)));
            if a>=0
                sig=[sig,1];
            else
                sig=[sig,0];
            end
            if b>=0
                sig_t=[sig_t,1];
            else
                sig_t=[sig_t,0];
            end
        end
        msg(msg==-1)=0;
        brror_sq=biterr(sig,msg','overall');
        brror_t=biterr(sig_t,msg','overall');
        brror_sq_snr=[brror_sq_snr,brror_sq];
        brror_t_snr=[brror_t_snr,brror_t];
        %------------------------------------------------------%
        for i=1:sze7
            a = mean(xcorr(onwa,y_square7(:,i)));
            b = mean(xcorr(triangle_weights,y_triangle7(:,i)));
            if a>=0
                sig7=[sig7,1];
            else
                sig7=[sig7,0];
            end
            if b>=0
                sig7_t=[sig7_t,1];
            else
                sig7_t=[sig7_t,0];
            end
        end
        fin_msg7_sq=decode(decoder7,sig7');
        fin_msg7_t=decode(decoder7,sig7_t');
        brror7_sq=biterr(fin_msg7_sq,msg,'overall');
        brror7_t=biterr(fin_msg7_t,msg,'overall');
        brror_sq_snr7=[brror_sq_snr7,brror7_sq];
        brror_t_snr7=[brror_t_snr7,brror7_t];
        
        %-------------------------------------------------------%
        for i=1:sze11
            a = mean(xcorr(onwa,y_square11(:,i)));
            b = mean(xcorr(triangle_weights,y_triangle11(:,i)));
            if a>=0
                sig11=[sig11,1];
            else
                sig11=[sig11,0];
            end
            if b>=0
                sig11_t=[sig11_t,1];
            else
                sig11_t=[sig11_t,0];
            end
        end
        fin_msg11_sq=decode(decoder11,sig11');
        fin_msg11_t=decode(decoder11,sig11_t');
        brror11_sq=biterr(fin_msg11_sq,msg,'overall');
        brror11_t=biterr(fin_msg11_t,msg,'overall');
        brror_sq_snr11=[brror_sq_snr11,brror11_sq];
        brror_t_snr11=[brror_t_snr11,brror11_t];
    end
    fin_brr_sq = fin_brr_sq+(brror_sq_snr/input);
    fin_brr_t = fin_brr_t+(brror_t_snr/input);
    fin_brr_sq11 = fin_brr_sq11+(brror_sq_snr11/input);
    fin_brr_t11 = fin_brr_t11+(brror_t_snr11/input);
    fin_brr_sq7 = fin_brr_sq7+(brror_sq_snr7/input);
    fin_brr_t7 = fin_brr_t7+(brror_t_snr7/input);
    
end
fin_brr_sq = fin_brr_sq/5;
fin_brr_t = fin_brr_t/5;
fin_brr_sq7 = fin_brr_sq7/5;
fin_brr_t7 = fin_brr_t7/5;
fin_brr_sq11 = fin_brr_sq11/5;
fin_brr_t11 = fin_brr_t11/5;

figure;
semilogy(snr,fin_brr_sq,snr,fin_brr_sq11,snr,fin_brr_sq7)
xlabel('SNR in dB');
ylabel('BER @semi-log scale');
grid on;
legend('Without error correction','with BCH(15,11)','with BCH(15,7)');
title('BER (semi-log scale) vs SNR in dB over 5 runs of size 10^5 (BPSK square wave)');
figure;
semilogy(snr,fin_brr_t,snr,fin_brr_t11,snr,fin_brr_t7)
xlabel('SNR in dB');
ylabel('BER @semi-log scale');
grid on;
legend('Without error correction','with BCH(15,11)','with BCH(15,7)');
title('BER (semi-log scale) vs SNR in dB over 5 runs of size 10^5 (BPSK triangular wave)');
