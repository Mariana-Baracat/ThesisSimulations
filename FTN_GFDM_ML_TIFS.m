% This simulation was developed by Mariana Baracat as part of her PhD research
% and was employed in Chapter 2.

clear all
close all
clc


S=5;
P=3;
N_hat=P*S;

K_hat=5;
M_hat=2.4;

vt=K_hat/S;
vf=M_hat/P;

K=floor(N_hat/M_hat);
M=floor(N_hat/K_hat);
N=K*M;

vf_real=S/K;
vt_real=P/M; 


T=2; 


hfsc=[1 0.4 0.2 0.08]
h=[hfsc zeros(1,N_hat-length(hfsc))];
H=zeros(N_hat,N_hat);
   

 for i_CH=1:(N_hat+Ncp)
     H(:,i_CH)= circshift(h,[0 (i_CH-1)]);
 end

g = zeros(P*S, 1);
g(1:S) = 1;
g = sqrt(vt_real*vf_real)*g.'/ sqrt(sum(g.*g));

E=vf_real*vt_real;


for m=1:M

    p=mod((m-1)*(vt*S),N_hat);
    for k=1:K
        
        A(:,(m-1)*K+(k-1)+1)= circshift(g,[0 p]).*exp((0:N_hat-1)*1j*2*pi*(k-1)*vf/S);
        
    end
end

 c = de2bi((0:T^(M*K)-1)',log2(T)*(M*K),'left-msb');
 cm= qammod(c.',T);
 f=A'*A*cm;

Bh=A';
n_erros = 100;
Erro=0;
ns=0;

for SNR=0:1:10

       
    Erro=0;
    erro_min=0;
    ns=0;
    
    
    SNR_L = 10.^(SNR./10); 
    sigma = sqrt(E./(2.*SNR_L)); 
   
     
    while(erro_min <= n_erros)
     
      
      s = randi(T,N,1)-1;
      d=qammod(s,T);
      
     
     
      x=A*d;

      w=sigma*(randn(size(x))+1j*randn(size(x)));
      y=H*x+w;

      y_eq= ifft(fft(y)./(fft(H(:,1)))); 
     
      r=Bh*y_eq;

          for u=1:T^(M*K)
                dmin(:,u)=sum((abs(f(:,u)-r)).^2);
          end
            
        [val ix] = min(dmin);

      erro_min=erro_min+sum(xor(s',de2bi(ix-1,log2(T)*(M*K),'left-msb')));
        
      ns=ns+N;
         
      
    end
    
SER_h = erro_min/ns;
semilogy(SNR,SER_h,'ob');
pause (0.1)
SER_h
SNR
hold on
end


gain=hfsc;
hold on

H       = fft(gain,S);
magH    = abs(H).^2;

SNR_L   = 10.^((0:10)./10);
pe      = zeros(size(SNR_L));
Rt=N_hat/(N_hat);

for k = 1:S
    
    eta   =(magH(k)*SNR_L*Rt);
    pe   =pe + (1/2*erfc(sqrt(eta)));
end
pb = pe/S;
semilogy(0:10,pb,'-r','Linewidth',2)
grid on

