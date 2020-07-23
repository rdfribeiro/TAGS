cd(dirname(@__FILE__));
include("TAGS.jl") #don't forget to let the TAGS.jl in the same directory
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
L = 11.; #length of each counterpoise (must be greater than 10.0)
prof = 0.5; #Grounding grid depth
CONDX, CONDY, CONDZ, Node, A, S = MC_Torre(r,L,prof,100*r); #getting data
ρ = 100.0; #soil resistivity
ϵr = 9.0; #soil permitivity
freq = logspace(log10(100),log10(2e6),200); #frequency points
nfreq = length(freq); #number of points
nInt = 300;#points in gauss legandre
Z = zeros(nfreq) + im*zeros(nfreq); #inicializating the impedance
ZA = zeros(nfreq) + im*zeros(nfreq); #inicializating the impedance
ncores = Threads.nthreads(); #using paralel computing (number of cores)
Fonte = zeros(Node[end,end],1); #current source (ensuring that 25% of the current go to each counterpoise)
Fonte[Node[1,1]] = 0.25;
Fonte[Node[3,1]] = 0.25;
Fonte[Node[5,1]] = 0.25;
Fonte[Node[7,1]] = 0.25;
println("Inicializating the so-called \"HEM\"")
begin time_original = @elapsed Threads.@threads for linha = 1:nfreq
        Yn = HEM(ρ,ϵr,r,CONDX,CONDY,CONDZ,A,S,freq[linha],nInt,"");
        Z[linha] = ((inv(Yn)*Fonte)[1,1]); # obtaining the Z of the first point, corner
        Yn = HEM(ρ,ϵr,r,CONDX,CONDY,CONDZ,A,S,freq[linha],nInt,"Alipio");
        ZA[linha] = ((inv(Yn)*Fonte)[1,1]); # obtaining the Z of the first point, corner
        if Threads.threadid() == 1
            println("$(100.0*linha/nfreq*ncores) % concluído ... Aguarde")
        end
        global Z, ZA
    end
end
#saving data
X = string(now());#string used to save the data
mkdir("Results_"*X[1:13]*X[15:16]) #creating directory to save data
cd(pwd()*"/Results_"*X[1:13]*X[15:16])
Nome = string(convert(Int64, ρ))*"_"* string(convert(Int64, L));
plot(freq,abs.(Z), label = "HEM", xscale = :log10, linecolor = :black,lw = 2,xlabel = "Frequency [Hz]",ylabel = "|z(ω)|");
plot!(freq,abs.(ZA), label = "HEM-Alipio", xscale = :log10, linecolor = :red,lw = 2,xlabel = "Frequency [Hz]",ylabel = "|z(ω)|");
savefig("Frequency_Response"*Nome*".png");
#saving data in a .mat structure
file = matopen("Grounding_Data.mat", "w");
write(file, "ZA", ZA); write(file, "Z", Z); write(file, "freq", freq); close(file);
cd(dirname(@__FILE__))
