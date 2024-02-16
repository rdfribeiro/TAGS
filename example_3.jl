cd(dirname(@__FILE__));
include("TAGS.jl") #don't forget to let the TAGS.jl in the same directory
using MAT
using Dates
using Plots
#Configuring plot font
#font = Plots.font("Arial", 16)
#font2 = Plots.font("Arial", 16)
#pyplot(guidefont=font, xtickfont=font, ytickfont=font, legendfont=font2);
plot();

#BEGIN -------------------------------------------------------------------------
r = 7.e-3; #Electrode radius
L = 11.; #length of horizontal electrode
prof = 0.5; #Depth
CONDX, CONDY, CONDZ, Node, A, S = MC_Horizontal(r,L,prof,100*r); #getting data
ρ = 100.0; #soil resistivity
ϵr = 9.0; #soil permitivity
freq = logspace(log10(100),log10(2e6),200); #frequency points
nfreq = length(freq); #number of points
nInt = 100;#points in gauss legandre
Z = zeros(nfreq) + im*zeros(nfreq); #inicializating the impedance
Zm = zeros(nfreq) + im*zeros(nfreq); #inicializating the impedance
ncores = Threads.nthreads(); #using paralel computing (number of cores)
println("Inicializating the so-called \"HEM\"")
begin time_original = @elapsed Threads.@threads for linha = 1:nfreq
        Yn = HEM(ρ,ϵr,r,CONDX,CONDY,CONDZ,A,S,freq[linha],nInt,"");
        if Threads.threadid() == 1
            println("$(100.0*linha/nfreq*ncores) % concluído ... Aguarde")
        end
        Z[linha] = (inv(Yn)[1,1]);
        global Z
    end
end
begin time_modified = @elapsed Threads.@threads for linha = 1:nfreq
        Yn = MHEM(ρ,ϵr,r,CONDX,CONDY,CONDZ,A,S,freq[linha],nInt,"");
        if Threads.threadid() == 1
            println("$(100.0*linha/nfreq*ncores) % concluído ... Aguarde")
        end
        Zm[linha] = (inv(Yn)[1,1]);
        global Zm
    end
end
#saving data
X = string(now());#string used to save the data
mkdir("Results_"*X[1:13]*X[15:16]) #creating directory to save data
cd(pwd()*"/Results_"*X[1:13]*X[15:16])
Nome = string(convert(Int64, ρ))*"_"* string(convert(Int64, L));
plot(freq,abs.(Z), label = "HEM", xscale = :log10, linecolor = :black,lw = 2,xlabel = "Frequency [Hz]",ylabel = "|z(ω)|");
plot!(freq,abs.(Zm), label = "MHEM", xscale = :log10, linecolor = :red,lw = 2,xlabel = "Frequency [Hz]",ylabel = "|z(ω)|");
savefig("Frequency_Response"*Nome*".png");
#saving data in a .mat structure
file = matopen("Grounding_Data.mat", "w"); write(file, "Zm", Zm);
write(file, "Z", Z); write(file, "freq", freq); close(file);
cd(dirname(@__FILE__))
println("\n");
println("Total time (HEM): $(time_original)\n")
println("Total time (MHEM): $(time_modified)\n")
