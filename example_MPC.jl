cd(dirname(@__FILE__));
include("TAGS.jl") #don't forget to let the TAGS.jl in the same directory
using Plots
#BEGIN -------------------------------------------------------------------------
r = 10.e-3; #Electrode radius
L = 10.; #length of horizontal electrode
prof = 0.5; #Depth
CONDX, CONDY, CONDZ, Node, A, S = MC_Malha(r, L, 5, prof, 100*r); #getting data
#CONDX, CONDY, CONDZ, Node, A, S = MC_Horizontal(r, L,0.5, 100*r);
#CONDX, CONDY, CONDZ, Node, A, S = MC_Vertical_Rods(r, 3., 1, 3., 0.5, 100*r);    
ρ = 1000.0; #soil resistivity
ϵr = 9.0; #soil permitivity
nInt = 100;#points in gauss legandre
println("================== Inicializando o HEM ===================")
tempohem = @elapsed Yn = HEM(ρ,ϵr,r,CONDX,CONDY,CONDZ,A,S,0.00001,nInt,"");
Z = (inv(Yn)[1,1]);
println("Impedância calculada por meio do HEM: ", Z, " Ω")
println("Tempo de execução do HEM: ", tempohem, " s")
println("=== Inicializando o método do Potencial Constante (PC) ===")
tempoPC = @elapsed Rt = PotencialConstante(ρ,r,CONDX,CONDY,CONDZ,nInt);
println("Impedância calculada por meio do PC: ", Rt, " Ω")
println("Tempo de execução do PC: ", tempoPC, " s")
