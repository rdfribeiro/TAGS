#Numerical Laplace Transform (NLT)
#This toolbox was created by Federal University of São João del-Rei (UFSJ) - BRAZIL
#more details at www.ufsj.edu.br/gatci
#Auxiliar Functions
using FFTW
#MATLAB linspace
function linspace(a::Float64,b::Float64,c::Int64)
    return collect(range(a,stop=b,length=c));
end
#s apresenta característica do tipo s = -i*c + 2*pi*freq
function NLT(f_t,t)
   n = length(t);# Number of SAMPLES
   if rem(n,2) != 0;
      t = t[1:end-1];
      f_t = f_t[1:end-1];
   end
   n = length(t);# Number of SAMPLES
   tmax = t[end];# Simulation Total Time
   T = tmax*1.0; #Maximum Simulation Time
   c = -log(0.001)/T;
   dt = T/n;# Time interval
   dw = 2*pi/(n*dt);# Frequency Interval
   tfreq = dt*linspace(0.,convert(Float64,n-1),n);#Vector containing the frequency samples
   kk = convert(Array{Float64},collect(0:n/2));
   sk = zeros(length(kk),1)+1im*zeros(length(kk),1);
   for linha = 1:length(kk)
      sk[linha] = -1im*c + dw*kk[linha]; #também está no pdf
   end
   #Fourier Parameters
   #Using Fast Fourier Transform (FFT)
   #Frequency Domain Signal
   input_freq = fft(exp.(-c*tfreq).*f_t).*dt;
   #SINGLE-SIDED AMPLITUDE
   input_freq = input_freq[1:(convert(Int64,floor(n/2)))+1];
   input_freq[2:end-1] = 2*input_freq[2:end-1];
   F_s = input_freq;
   s = sk;
   return F_s, s;
end
#Inverse NLT
function INLT(F_s,s,Filter::String)
   #Filter. If Filter = ""; without filter;
   #If filter = "Hamming"; HAMMING Filter
   #If filter = "Blackman"; BLACKMAN FILTER
   ns = 2*(length(s)-1);#NUMERO DE AMOSTRAS NA FREQUENCIA
   #COLOCANDO UM FILTRO ...
   k = linspace(0.,ns/2,convert(Int64,ns/2+1));
   #HAMMING FILTER
   if Filter == "Hamming"
      alfa = 0.54;
      w = alfa .+ (1-alfa)*cos.(k*pi/(ns/2));
      F_s = w.*F_s;
   end
   #BLACKMAN FILTER
   if Filter == "Blackman"
      alfa = 0.5;
      beta = 0.08;
      w = (alfa-beta) .+ alfa*cos.(pi*k/(ns/2))+beta*cos.(2*pi*k/(ns/2));
      F_s = w.*F_s;
   end
   F_s = [F_s[1]; F_s[2:convert(Int64,ns/2)+1]/2; (real(F_s[convert(Int64,ns/2):-1:2]/2)-1im*imag(F_s[convert(Int64,ns/2):-1:2]/2))];
   c = - imag(s[1]);
   T = -log(0.001)/c;#Total time
   dt = T/ns; #Time interval
   t = dt*linspace(0.,convert(Float64,ns-1),ns);#Time information vector
   f_t = exp.(c*t)/dt.*ifft(F_s);# Time response (inverse transform)
   f_t = real(f_t); # ensuring that the values are all real
   return f_t, t;
end
