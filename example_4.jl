cd(dirname(@__FILE__));
include("TAGS.jl") #don't forget to let the TAGS.jl in the same directory
include("NLT.jl") #don't forget to let the TAGS.jl in the same directory
using MAT
using Dates
using Plots
#Configuring plot font
font = Plots.font("Arial", 16)
font2 = Plots.font("Arial", 16)
pyplot(guidefont=font, xtickfont=font, ytickfont=font, legendfont=font2);
plot();
#BEGIN -------------------------------------------------------------------------
r = 7.e-3; #Electrode radius
L = 11.; #length of horizontal electrode
prof = 0.5; #Depth
CONDX, CONDY, CONDZ, Node, A, S = MC_Horizontal(r,L,prof,100*r); #getting data
ρ = 100.0; #soil resistivity
ϵr = 9.0; #soil permitivity
nInt = 100;#points in gauss legandre
#Excitation data   -------------------------------------------------------------
tmax = 40.e-6; #Simulation total time
t = linspace(0.,tmax,2^8); #time sample
#San Salvatore Station data
#Heidler's Function
#First Stroke
i0 = [3.00 4.50 3.00 3.80 13.6 11.0 05.7]*10^3;
nkk  = [2.00 3.00 5.00 7.00 44.0 2.00 15.0];
t1 = [3.00 3.50 5.20 6.00 6.60 100.0 11.7]*10^-6;
t2 = [76.0 25.0 20.0 60.0 60.0 600. 48.5]*10^-6;
#Subsequent Stroke
# i0 = [10.7 6.5]*10^3;
# nkk  = [2.00 2.00];
# t1 = [.25 2.1]*10^-6;
# t2 = [2.5 230]*10^-6;
current = zeros(length(t));
nH = length(nkk);
for linha = 1:nH
    global current;
    nk = nkk[linha];
    t1k = t1[linha];
    t2k = t2[linha];
    i0k = i0[linha];
    NiK = exp(-(t1k/t2k)*((nk*t2k/t1k)^(1. /nk)));
    currentb = (i0k/NiK)*exp.(-t/t2k).*(((t/t1k).^nk)./(1. .+((t/t1k).^nk)));
    current = current + currentb;
end
#Using NLT
f_t = current; F_s, s = NLT(f_t,t);
f_t_regenerada, t = INLT(F_s,s,"");
plot(10^6*t,10^-3*f_t, label = "Original Function", linecolor = :black,lw = 2);
plot!(10^6*t,10^-3*f_t_regenerada, label = "Regenerated Function", linecolor = :red,lw = 2,linestyle = :dash, xlabel = "Time [μs]",ylabel = "Current [kA]");

savefig("Injected_Current.png");
nfreq = length(s); V_w = zeros(nfreq) + im*zeros(nfreq);
ncores = Threads.nthreads(); #using paralel computing (number of cores)
println("Inicializating the so-called \"HEM\"")
begin time_original = @elapsed Threads.@threads for linha = 1:nfreq
        Yn = HEM(ρ,ϵr,r,CONDX,CONDY,CONDZ,A,S,s[linha]/(2*pi),nInt,"");
        if Threads.threadid() == 1
            println("$(100.0*linha/nfreq*ncores) % concluído ... Aguarde")
        end
        Z = (inv(Yn)[1,1]);
        V_w[linha] = Z*F_s[linha];
        global V_w
    end
end
v_t,t = INLT(V_w,s,"");
begin time_modified = @elapsed Threads.@threads for linha = 1:nfreq
        Yn = MHEM(ρ,ϵr,r,CONDX,CONDY,CONDZ,A,S,s[linha]/(2*pi),nInt,"");
        if Threads.threadid() == 1
            println("$(100.0*linha/nfreq*ncores) % concluído ... Aguarde")
        end
        Z = (inv(Yn)[1,1]);
        V_w[linha] = Z*F_s[linha];
        global V_w
    end
end
v_t_m,t = INLT(V_w,s,"");
#saving data
X = string(now());#string used to save the data
mkdir("Results_"*X[1:13]*X[15:16]) #creating directory to save data
cd(pwd()*"/Results_"*X[1:13]*X[15:16])
Nome = string(convert(Int64, ρ))*"_"* string(convert(Int64, L));
plot(10^6*t,10^-3*v_t, label = "HEM", linecolor = :black,lw = 2);
plot!(10^6*t,10^-3*v_t_m, label = "MHEM", linecolor = :red,lw = 2, xlabel = "Time [μs]",ylabel = "GPR [kV]");
savefig("Frequency_Response"*Nome*".png");
cd(dirname(@__FILE__))
println("\n");
println("Total time (HEM): $(time_original)\n")
println("Total time (MHEM): $(time_modified)\n")
