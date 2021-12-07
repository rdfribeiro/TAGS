#Transients Analysis of Grounding Systems (TAGS) PROJECT
#This toolbox was created by Federal University of São João del-Rei (UFSJ) - BRAZIL
#more details at www.ufsj.edu.br/gatci
#Auxiliar Functions
#MATLAB linspace
function linspace(a::Float64,b::Float64,c::Int64)
    return collect(range(a,stop=b,length=c));
end
#MATLAB logspace
function logspace(a::Float64,b::Float64,c::Int64)
    return exp10.(range(a, stop=b, length=c))
end
#MATLAB max function
function Base.max(a)
    return maximum(a)
end
#Gauss Legandre weight
function Pesos_GL(n::Int)
    if (rem(n,2) != 0)
        print("o número de pontos de integração deve ser PAR")
        n = n+1;
    end
    T = Array{Float64}(undef,n,);
    A = Array{Float64}(undef,n,);
    z1 = Inf;
    tol = 1.e-15;
    pp = 0.0;
    m = Int(floor(0.5*(n+1)));
    for linha = 1:m
        z = cos(pi*(linha - 0.25)/(n+0.5));
        while abs(z - z1) > tol
            p1 = 1;
            p2 = 0;
            for coluna = 1:n
                p3 = p2;
                p2 = p1;
                p1 = ((2*coluna-1)*z*p2-(coluna-1)*p3)/coluna;
            end
            pp = n*(z*p1-p2)/(z^2-1);
            z1 = z;
            z = z1 - p1/pp;
        end
        T[Int(n/2)+m+1-linha] = z;
        A[Int(n/2)+m+1-linha] = 2/((1-z^2)*pp^2);
        T[linha] = -z;
        A[linha] = 2/((1-z^2)*pp^2);
    end
    return T, A;
end
#Alipio's soil model ("Alipio")
function ps_alipio(ro::Float64,f)
    eps0 = 8.854*10^-12;
    sigma0 = 1/ro;
    sigma0 = sigma0*1000;
    hSig = 1.26*sigma0^-0.73;
    Gama = 0.54;
    Erinf = 12;
    Sig = sigma0 + sigma0*hSig*(f/1e6)^Gama;
    Sig = Sig/1000;
    #=Galip=1/(Roalip);
    #CALCULO DA PERMISSIVIDADE RELATIVA SEGUNDO METODOLOGIA DE ALIPIO
    if f<=10000
        f=10000; %so é valido para f >= 10 kHz
    end=#
    Epsr = Erinf + tan(pi*Gama/2)*1e-3/(2*pi*eps0*(1e6)^Gama)*sigma0*hSig*f^(Gama-1);
    roc = 1/Sig;
    roc, Epsr
end
#HEM Library
#HEM function - Tradicional - Numerical Solution by Gauss Legandre
function HEM(ro::Float64,epsr::Float64,rf::Float64,CONDUTORESx,CONDUTORESy,CONDUTORESz,A,S,freq,nInt::Int64,SoilParameter::String)
    #corrigindo os parâmetros do solo c/freq.
	if SoilParameter == "Alipio"
		ro, epsr = ps_alipio(ro,freq);
	end
    #variáveis eletromagnéticas
    mi0 = 4*pi*10^-7;
	eps0 = 8.854*10^-12;
	sigma = 1/ro;
    #variaveis espaciais e inicialização de variaveis
    nseg = size(CONDUTORESx,1); pfreq = length(freq);
    Zl = zeros(nseg,nseg) + zeros(nseg,nseg)*im;
    Zt = zeros(nseg,nseg) + zeros(nseg,nseg)*im;
    line = 0; colum = 0; lin = 0; col = 0;
    aux1r = 0.; aux1i = 0.; aux2r = 0.; aux2i = 0.;
    Dot = 0.; delefont = 0.; delerece = 0.;
    XR0 = 0.; YR0 = 0.; ZR0 = 0.; XS0 = 0.; YS0 = 0.; ZS0 = 0.;
    XR1 = 0.; YR1 = 0.; ZR1 = 0.; XS1 = 0.; YS1 = 0.; ZS1 = 0.;
    xfonte = 0.; yfonte = 0.; zfonte = 0.; yfontei = 0.;
    w = 0.; funcaor = 0.; funcaoi = 0.; Intl = 0.; Intt = 0.; Intli = 0.; Intti = 0.; M =0.;
    K = 0.0 + 0.0im;
    #Pesos do Gauss Legandre
	if rem(nInt,2) != 0
		nInt = nInt + 1;
	end
    tgl1, Agl1 = Pesos_GL(nInt);
	ntgl2 = Int(round(nInt/4));
	#ensuring it is even number
	if rem(ntgl2,2) != 0
		ntgl2 = ntgl2 + 1;
	end
	tgl2, Agl2 = Pesos_GL(ntgl2);
	ntgl1 = length(tgl1); X = size(S,2);
    Ynodal = Array{Complex{Float64}}(undef,X,X);
    #Loop principal
    w = 2*pi*freq;
    eps = eps0*epsr;
    K = sqrt(1im*w*mi0*(sigma+1im*w*eps));
    Tal_rl = 1 + 0.0im;
    Tal_rt = 1 + 0.0im;

    for line = 1:(nseg)
        XR0 = CONDUTORESx[line,1]; XR1 = CONDUTORESx[line,2];
        YR0 = CONDUTORESy[line,1]; YR1 = CONDUTORESy[line,2];
        ZR0 = CONDUTORESz[line,1]; ZR1 = CONDUTORESz[line,2];

        for colum = line:(nseg)
            XS0 = CONDUTORESx[colum,1]; XS1 = CONDUTORESx[colum,2]; XSmed = (XS0 + XS1)/2;
            YS0 = CONDUTORESy[colum,1]; YS1 = CONDUTORESy[colum,2]; YSmed = (YS0 + YS1)/2;
            ZS0 = CONDUTORESz[colum,1]; ZS1 = CONDUTORESz[colum,2]; ZSmed = (ZS0 + ZS1)/2;

            delefont = sqrt((XS1-XS0)^2 + (YS1-YS0)^2 + (ZS1-ZS0)^2); #DIFERENCIAL DO dl sqrt(x'(t)^2+y'(t)^2 ... CONFERIR EM COMO RESOLVER INTEGRAL DE LINHA PARA COMPREENDER
			delerece = sqrt((XR1-XR0)^2 + (YR1-YR0)^2 + (ZR1-ZR0)^2);
            Dot = ((XS1-XS0)*(XR1-XR0)+(YS1-YS0)*(YR1-YR0)+(ZS1-ZS0)*(ZR1-ZR0))/(delefont*delerece); #Projeção (associado com produtor vetorial entre "dl's")

            aux1r = 0.; aux1i = 0.;
			if line == colum
				ntgl = ntgl1;
				tgl = tgl1;
				Agl = Agl1;
			else
				ntgl = ntgl2;
				tgl = tgl2;
				Agl = Agl2;
			end
            for lin = 1:ntgl
                aux2r = 0.; aux2i = 0.;

                xfonte = (XS0+XS1)/2+((-XS0+XS1))/2*tgl[lin];
                yfonte = (YS0+YS1)/2+((-YS0+YS1))/2*tgl[lin];
                zfonte = (ZS0+ZS1)/2+((-ZS0+ZS1))/2*tgl[lin];
                yfontei = (-(YS0+YS1))/2+((-(-YS0+YS1)))/2*tgl[lin];
                for col = 1:ntgl
                    XR0ece = (XR0+XR1)/2+((-XR0+XR1))/2*tgl[col];
                    YR1ece = (YR0+YR1)/2+((-YR0+YR1))/2*tgl[col];
                    zrece = (ZR0+ZR1)/2+((-ZR0+ZR1))/2*tgl[col];

                    rreal = sqrt((xfonte-XR0ece)^2+(yfonte-YR1ece)^2+(abs(zfonte-zrece) + rf)^2);
                    rimag = sqrt((xfonte-XR0ece)^2+(yfontei-YR1ece)^2+(abs(zfonte-zrece) + rf)^2);

                    funcaor = exp(-K*rreal)/rreal;
                    funcaoi = exp(-K*rimag)/rimag;

                    aux2r = aux2r + Agl[col]*funcaor;
                    aux2i = aux2i + Agl[col]*funcaoi;
                end
                aux1r = aux1r + Agl[lin]*aux2r;
                aux1i = aux1i + Agl[lin]*aux2i;
            end
            Intl = 0.25*(aux1r + Tal_rl*aux1i)*(delerece*delefont)*Dot; #AINDA TEM QUE COLOCAR CONSTANTE DE REFLEXÃO
            Intt = 0.25*(aux1r + Tal_rt*aux1i)*(delerece*delefont);

            Zl[line,colum] = (1im*w*mi0/(4*pi))*Intl; #ARMAZENANDO IMPEDANCIA LONGITUDINAL
            Zt[line,colum] = (1/(4*pi*(sigma+1im*w*eps)*(delerece*delefont)))*Intt; #ARMAZENANDO IMPEDANCIA TRANSVERSAL

            Zl[colum,line] = Zl[line,colum];
            Zt[colum,line] = Zt[line,colum];
        end
    end
    Ynodal = (transpose(A)/(Zt))*A + (transpose(S)/(Zl))*S;
    return Ynodal;
end
#MHEM function - Modofied version - Numerical Solution by Gauss Legandre
#função do MHEM
function MHEM(ro::Float64,epsr::Float64,rf::Float64,CONDUTORESx,CONDUTORESy,CONDUTORESz,A,S,freq,nInt::Int64,SoilParameter::String)
    #corrigindo os parâmetros do solo c/freq.
	if SoilParameter == "Alipio"
		ro, epsr = ps_alipio(ro,freq);
	end
    #variáveis eletromagnéticas
    mi0 = 4*pi*10^-7;
    eps0 = 8.854*10^-12;
    sigma = 1/ro;
    #variaveis espaciais e inicialização de variaveis
    XX = size(CONDUTORESx);
    nseg = XX[1]; pfreq = length(freq);

    Zl = Array{Complex{Float64}}(undef,nseg,nseg);
    Zt = Array{Complex{Float64}}(undef,nseg,nseg);

    line = 0; colum = 0; lin = 0;
    auxr = 0.; auxi = 0.; Dot = 0.; delefont = 0.; delerece = 0.;
    XR0 = 0.; YR0 = 0.; ZR0 = 0.; XS0 = 0.; YS0 = 0.; ZS0 = 0.;
    xfonte = 0.; yfonte = 0.; zfonte = 0.; R1 = 0.; R2 = 0.; R1i = 0.;
    R2i = 0.; funcaor = 0.; funcaoi = 0.; Intl = 0.; Intt = 0.; Intli = 0.; Intti = 0.; M =0.;
    K = 0.0 + 0.0im;
	#Pesos do Gauss Legandre
	if rem(nInt,2) != 0
		nInt = nInt + 1;
	end
    tgl1, Agl1 = Pesos_GL(nInt);
	ntgl2 = Int(round(nInt/4));
	#ensuring it is even number
	if rem(ntgl2,2) != 0
		ntgl2 = ntgl2 + 1;
	end
	tgl2, Agl2 = Pesos_GL(ntgl2);
	ntgl1 = length(tgl1); X = size(S,2);
    Ynodal = Array{Complex{Float64}}(undef,X,X);
    #Loop principal
    w = 2*pi*freq;
    eps = epsr*eps0;
    K = sqrt(1im*w*mi0*(sigma+1im*w*eps));
    for line = 1:(nseg)
        XR0 = CONDUTORESx[line,1]; XR1 = CONDUTORESx[line,2]; XRmed = (XR0 + XR1)/2;
        YR0 = CONDUTORESy[line,1]; YR1 = CONDUTORESy[line,2]; YRmed = (YR0 + YR1)/2;
        ZR0 = CONDUTORESz[line,1]; ZR1 = CONDUTORESz[line,2]; ZRmed = (ZR0 + ZR1)/2;

        for colum = line:(nseg)
            XS0 = CONDUTORESx[colum,1]; XS1 = CONDUTORESx[colum,2]; XSmed = (XS0 + XS1)/2;
            YS0 = CONDUTORESy[colum,1]; YS1 = CONDUTORESy[colum,2]; YSmed = (YS0 + YS1)/2;
            ZS0 = CONDUTORESz[colum,1]; ZS1 = CONDUTORESz[colum,2]; ZSmed = (ZS0 + ZS1)/2;

            RBB = sqrt((XRmed-XSmed)^2+(YRmed-YSmed)^2+(abs(ZRmed-ZSmed)+rf)^2);
            RBBi = sqrt((XRmed-XSmed)^2+(-YRmed-YSmed)^2+(abs(ZRmed-ZSmed)+rf)^2);
            delefont = sqrt((XS1-XS0)^2 + (YS1-YS0)^2 + (ZS1-ZS0)^2); #DIFERENCIAL DO dl sqrt(x'(t)^2+y'(t)^2 ... CONFERIR EM COMO RESOLVER INTEGRAL DE LINHA PARA COMPREENDER
            delerece = sqrt((XR1-XR0)^2 + (YR1-YR0)^2 + (ZR1-ZR0)^2);
            Dot = ((XS1-XS0)*(XR1-XR0)+(YS1-YS0)*(YR1-YR0)+(ZS1-ZS0)*(ZR1-ZR0))/(delefont*delerece); #Projeção (associado com produtor vetorial entre "dl's")
            auxr = 0.0; auxi = 0.0;
			if line == colum
				ntgl = ntgl1;
				tgl = tgl1;
				Agl = Agl1;
			else
				ntgl = ntgl2;
				tgl = tgl2;
				Agl = Agl2;
			end
            for lin = 1:length(tgl)
                xfonte = (XS0+XS1)/2+((-XS0+XS1))/2*tgl[lin];
                yfonte = (YS0+YS1)/2+((-YS0+YS1))/2*tgl[lin];
                zfonte = (ZS0+ZS1)/2+((-ZS0+ZS1))/2*tgl[lin];

                R1 = sqrt((xfonte-XR0)^2+(yfonte-YR0)^2+(abs(zfonte-ZR0)+rf)^2);
                R2 = sqrt((xfonte-XR1)^2+(yfonte-YR1)^2+(abs(zfonte-ZR1)+rf)^2);

                R1i = sqrt((xfonte-XR0)^2+(yfonte+YR0)^2+(abs(zfonte-ZR0)+rf)^2);
                R2i = sqrt((xfonte-XR1)^2+(yfonte+YR1)^2+(abs(zfonte-ZR1)+rf)^2);

                funcaor = log(Complex((R1+R2+delefont)/(R1+R2-delefont)));
                funcaoi = log(Complex((R1i+R2i+delefont)/(R1i+R2i-delefont)));

                auxr = auxr + Agl[lin]*funcaor;
                auxi = auxi + Agl[lin]*funcaoi;
            end
            Intl = 0.5*(auxr)*delerece*Dot;
            Intt = 0.5*(auxr)*delerece;
            Intli = 0.5*(auxi)*delerece*Dot;
            Intti = 0.5*(auxi)*delerece;

            if line == colum
                M = 2*(log((1+sqrt(1+(rf/delefont)^2))/(rf/delefont)) - sqrt(1+(rf/delefont)^2) + rf/delefont);
                Intl = delefont*M;
                Intt = Intl;
            end

            Zl[line,colum] = (1im*w*mi0/(4*pi))*(Intl*exp(-K*RBB) + Intli*exp(-K*RBBi));
            Zl[colum,line] = Zl[line,colum];

            Zt[line,colum] = (1/(4*pi*(sigma+1im*w*eps)*delefont*delerece))*(Intt*exp(-K*RBB) + Intti*exp(-K*RBBi));
            Zt[colum,line] = Zt[line,colum];
        end
    end
    Ynodal[:,:] = (transpose(A)/(Zt))*A + (transpose(S)/(Zl))*S;
    return Ynodal;
end
"""
Contribution from Fernando Lima Viana Costa
PotencialConstante(ρ::Float64,rf::Float64,CONDUTORESx::Array{Float64,2},CONDUTORESy::Array{Float64,2},CONDUTORESz::Array{Float64,2},nInt::Int64)

Calcula a impedância equivalente de uma estrutura de aterramento utilizando o método do
potencial constante

## Entradas
1. ρ: resistividade do solo
2. raio: raio dos eletrodos do aterramento
3. CONDUTORESx: matriz contendo as coordenadas na direção x inicial [1] e final [2] de cada
eletrodo do aterramento
4. CONDUTORESy: análogo a CONDUTORESx, para as coordenadas da direção y
5. CONDUTORESz: análogo a CONDUTORESx, para as coordenadas da direção z
6. nInt: número de pontos de integração

## Saídas
A resistência equivalente da estrutura de aterramento.

"""
function PotencialConstante(ρ::Float64,rf::Float64,
                            CONDUTORESx::Array{Float64,2},
                            CONDUTORESy::Array{Float64,2},
                            CONDUTORESz::Array{Float64,2},
                            nInt::Int64)
    #variáveis eletromagnéticas
    #μ0 = 4*pi*10^-7;
    #ϵ0 = 8.854*10^-12;
    #σ = 1/ro;
    #variaveis espaciais e inicialização de variaveis
    nseg = size(CONDUTORESx,1);
    Rnm = zeros(nseg,nseg);
    #Pesos do Gauss Legandre
	if rem(nInt,2) != 0
		nInt = nInt + 1;
	end
    tgl1, Agl1 = Pesos_GL(nInt);
    ntgl2 = Int(round(nInt/4));
    #ensuring it is even number
    if rem(ntgl2,2) != 0
        ntgl2 = ntgl2 + 1;
    end
    tgl2, Agl2 = Pesos_GL(ntgl2);
    ntgl1 = length(tgl1);

    #Loop principal
    Tal_rt = 1;

    for line = 1:(nseg)
        XR0 = CONDUTORESx[line,1]; XR1 = CONDUTORESx[line,2];
        YR0 = CONDUTORESy[line,1]; YR1 = CONDUTORESy[line,2];
        ZR0 = CONDUTORESz[line,1]; ZR1 = CONDUTORESz[line,2];

        for colum = line:(nseg)
            XS0 = CONDUTORESx[colum,1]; XS1 = CONDUTORESx[colum,2]; XSmed = (XS0 + XS1)/2;
            YS0 = CONDUTORESy[colum,1]; YS1 = CONDUTORESy[colum,2]; YSmed = (YS0 + YS1)/2;
            ZS0 = CONDUTORESz[colum,1]; ZS1 = CONDUTORESz[colum,2]; ZSmed = (ZS0 + ZS1)/2;

            aux1r = 0.; aux1i = 0.;
            if line == colum
                ntgl = ntgl1;
                tgl = tgl1;
                Agl = Agl1;
            else
                ntgl = ntgl2;
                tgl = tgl2;
                Agl = Agl2;
            end
            for lin = 1:ntgl
                aux2r = 0.; aux2i = 0.;

                xfonte = (XS0+XS1)/2+((-XS0+XS1))/2*tgl[lin];
                yfonte = (YS0+YS1)/2+((-YS0+YS1))/2*tgl[lin];
                zfonte = (ZS0+ZS1)/2+((-ZS0+ZS1))/2*tgl[lin];
                yfontei = (-(YS0+YS1))/2+((-(-YS0+YS1)))/2*tgl[lin];
                for col = 1:ntgl
                    XR0ece = (XR0+XR1)/2+((-XR0+XR1))/2*tgl[col];
                    YR1ece = (YR0+YR1)/2+((-YR0+YR1))/2*tgl[col];
                    zrece = (ZR0+ZR1)/2+((-ZR0+ZR1))/2*tgl[col];

                    rreal = sqrt((xfonte-XR0ece)^2+(yfonte-YR1ece)^2+(abs(zfonte-zrece) + rf)^2);
                    rimag = sqrt((xfonte-XR0ece)^2+(yfontei-YR1ece)^2+(abs(zfonte-zrece) + rf)^2);

                    funcaor = 1/rreal;
                    funcaoi = 1/rimag;

                    aux2r = aux2r + Agl[col]*funcaor;
                    aux2i = aux2i + Agl[col]*funcaoi;
                end
                aux1r = aux1r + Agl[lin]*aux2r;
                aux1i = aux1i + Agl[lin]*aux2i;
            end
            ∫ = 0.25*(aux1r + Tal_rt*aux1i);
            Rnm[line,colum] = ∫*ρ/(4*pi); 
            Rnm[colum,line] = Rnm[line,colum];
        end
    end
    return inv(sum(inv(Rnm)));
end
#Grid Data Generator
#Function that generate the counterpoise data
function MC_Torre(r::Float64,L::Float64,h::Float64,segR::Float64)
    n = floor(Int, L/(segR)); #Número de segmentações de cada "eletrodo"
    n_aux = 8;
    Pos = Array{Float64}(undef,n_aux,10); #inicializando matriz auxiliar
    CONDX = Array{Float64}(undef,n*(n_aux),2); #inicialização posicao do vetor x
    CONDY = Array{Float64}(undef,n*(n_aux),2); #inicialização posicao do vetor y
    CONDZ = Array{Float64}(undef,n*(n_aux),2); #inicialização posicao do vetor z
    Nos = Array{Int64}(undef,n*(n_aux),2); #inicialização posicao dos nós que conectam cada segmento
    Node = Array{Int64}(undef,n_aux,2);
    Nos = Array{Int64}(undef,n*(n_aux),2);
    Condutores = Array{Float64}(undef,n+1,);
    #Posicionando informações da malha na matriz
    #inix fimx iniy fimy iniz fimz no (ini/fim) no (ini/fim) ini = 1 fim = 2
    Pos = [ -3. -10. -h -h   3.  10. 99 2  2 1;#1
           -10.   -L -h -h  10.  10.  1 2 99 1;#2
             3.  10. -h -h   3.  10. 99 2  4 1;#3
            10.    L -h -h  10.  10.  3 2 99 1;#4
            -3. -10. -h -h  -3. -10. 99 2  6 1;#5
           -10.   -L -h -h -10. -10.  5 2 99 1;#6
             3.  10. -h -h  -3. -10. 99 2  8 1;#7
            10.    L -h -h -10. -10.  7 2 99 1];#8
    Last_node = 1; #inicializando variável que computa o último nó do Loop
    for linha = 1:length(Pos[:,1])
        #direção x
        Condutores = linspace(Pos[linha,1],Pos[linha,2],n + 1);
        CONDX[((linha-1)*n) + 1:(linha)*n,:] = [Condutores[1:end-1] Condutores[2:end]];
        #CONDX = [CONDX; [CondutoresX(1:end-1).' CondutoresX(2:end).']];
        #direção y
        Condutores = linspace(Pos[linha,3],Pos[linha,4],n + 1);
        CONDY[((linha-1)*n) + 1:(linha)*n,:] = [Condutores[1:end-1] Condutores[2:end]];
        #CONDY = [CONDY; [CondutoresY(1:end-1).' CondutoresY(2:end).']];
        #direção z
        Condutores = linspace(Pos[linha,5],Pos[linha,6],n + 1);
        CONDZ[((linha-1)*n) + 1:(linha)*n,:] = [Condutores[1:end-1] Condutores[2:end]];
        #CONDZ = [CONDZ; [CondutoresZ(1:end-1).' CondutoresZ(2:end).']];
        if Pos[linha,7] > linha && Pos[linha,9] > linha #PRIMEIRO ELETRODO
            No = collect(Last_node:1:(Last_node + n));
        end
        if Pos[linha,7] < linha && Pos[linha,9] > linha
            No = collect(Last_node:1:(Last_node + n-1));
            No = [Node[Int(Pos[linha,7]),Int(Pos[linha,8])] No'];
        end
        if Pos[linha,7] > linha && Pos[linha,9] < linha
            No = collect(Last_node:1:(Last_node + n-1));
            No = [No Node[Pos[linha,9],Pos[linha,10]]];
        end
        if Pos[linha,7] < linha && Pos[linha,9] < linha
            No = collect(Last_node:1:(Last_node + n-2));
            No = [Node[Pos[linha,7],Pos[linha,8]] No Node[Pos[linha,9],Pos[linha,10]]];
        end
        Nos[((linha-1)*n)+1:(linha)*n,:] = [No[1:end-1] No[2:end]];
        #Nos = [Nos; [No[1:end-1]' No[2:end]']];
        Node[linha,:] = [No[1] No[end]];
        Last_node = max(No) + 1;
    end
    #MATRIZES DE CONECTIVIDADE
    A = zeros(n*(n_aux),Last_node - 1);
    S = zeros(n*(n_aux),Last_node - 1);
    #S = A;
    for linha = 1:length(CONDX[:,1]) #NUMERO DE ELETRODOS
        A[linha,Nos[linha,1]] = 0.5;
        S[linha,Nos[linha,1]] = 1.0;
        A[linha,Nos[linha,2]] = 0.5;
        S[linha,Nos[linha,2]] = -1.0;
    end
    #=writedlm("CONDX.csv",CONDX,',');
    writedlm("CONDY.csv",CONDY,',');
    writedlm("CONDZ.csv",CONDZ,',');
    writedlm("Node.csv",Node,',');
    writedlm("A.csv",A,',');
    writedlm("S.csv",S,',');=#

    plot([CONDX[:,1]' CONDX[end,2]]',[CONDZ[:,1]' CONDZ[end,2]]',linestyle =
        :dot,seriestype=:scatter,title="Grounding Grid",legend = false, aspect_ratio = 1);
    savefig("Grounding_Grid_Design.pdf");
    return CONDX, CONDY, CONDZ, Node, A, S
end
#Function that generate the counterpoise data (updated)
function MC_Torre_2(r::Float64,L::Float64,h::Float64,segR::Float64)
    n = floor(Int, L/(segR)); #Número de segmentações de cada "eletrodo"
    n_aux = 4; #número de hastes horizontais (cabos contrapeso)
    Pos = Array{Float64}(undef,n_aux,10); #inicializando matriz auxiliar
    CONDX = Array{Float64}(undef,n*(n_aux),2); #inicialização posicao do vetor x
    CONDY = Array{Float64}(undef,n*(n_aux),2); #inicialização posicao do vetor y
    CONDZ = Array{Float64}(undef,n*(n_aux),2); #inicialização posicao do vetor z
    Nos = Array{Int64}(undef,n*(n_aux),2); #inicialização posicao dos nós que conectam cada segmento
    Node = Array{Int64}(undef,n_aux,2);
    Nos = Array{Int64}(undef,n*(n_aux),2);
    Condutores = Array{Float64}(undef,n+1,);
    #Posicionando informações da malha na matriz
    #      inix prof fimz
    Pos = [-3.  -h   3.  10.;#1
            3.  -h   3.  10.;#2
		   -3.  -h  -3. -10.;#3
 	        3.  -h  -3. -10.]#4
    Last_node = 1; #inicializando variável que computa o último nó do Loop
	ΔPos = L/n;
	Δx = 0.0; Δz = 0.0;
	mm = 1;
	for col = 1:4
		posx = Pos[col,1];
		posy = Pos[col,2];
		posz = Pos[col,3];
		poszf = Pos[col,4];
		if posx > 0.0;
			Δx = ΔPos*sqrt(2)/2;
		else
			Δx = -ΔPos*sqrt(2)/2;
		end
		if posz > 0.0
			Δz = ΔPos*sqrt(2)/2;
		else
			Δz = -ΔPos*sqrt(2)/2;
		end
		for linha = 1:n
			CONDY[mm,1] = posy;
			CONDY[mm,2] = posy;
			CONDX[mm,1] = posx;
			CONDZ[mm,1] = posz;
			posx = posx + Δx;
			posz = posz + Δz;
			CONDX[mm,2] = posx;
			CONDZ[mm,2] = posz;
			if abs(posz) >= abs(poszf) && Δz != 0.0
				Δx = Δx*2/sqrt(2);
				Δz = 0.0;
			end
			Nos[mm,1] = mm + col - 1;
			Nos[mm,2] = mm + col;
			if linha == 1
				Node[col,1] = Nos[mm,1]
			end
			if linha == n
				Node[col,2] = Nos[mm,2]
			end
			mm = mm + 1;
		end
	end
	Last_node = Nos[end,2] + 1;
    #MATRIZES DE CONECTIVIDADE
    A = zeros(n*(n_aux),Last_node - 1);
    S = zeros(n*(n_aux),Last_node - 1);
    #S = A;
    for linha = 1:length(CONDX[:,1]) #NUMERO DE ELETRODOS
        A[linha,Nos[linha,1]] = 0.5;
        S[linha,Nos[linha,1]] = 1.0;
        A[linha,Nos[linha,2]] = 0.5;
        S[linha,Nos[linha,2]] = -1.0;
    end
    #=writedlm("CONDX.csv",CONDX,',');
    writedlm("CONDY.csv",CONDY,',');
    writedlm("CONDZ.csv",CONDZ,',');
    writedlm("Node.csv",Node,',');
    writedlm("A.csv",A,',');
    writedlm("S.csv",S,',');=#

    plot([CONDX[:,1]' CONDX[end,2]]',[CONDZ[:,1]' CONDZ[end,2]]',linestyle =
        :dot,seriestype=:scatter,title="Grounding Grid",legend = false, aspect_ratio = 1);
    savefig("Grounding_Grid_Design.pdf");
    return CONDX, CONDY, CONDZ, Node, A, S
end
#Function that generate square grid data
function MC_Malha(r::Float64,L::Float64,nMalha::Int64,h::Float64,segR::Float64)
    n = floor(Int, L/(segR)); #Número de segmentações de cada "eletrodo"
    n_aux = nMalha*(nMalha + 1);
    Pos1 = Array{Float64}(undef,n_aux,10); #inicializando matriz auxiliar 1
    Pos2 = Array{Float64}(undef,n_aux,10); #inicializando matriz auxiliar 2
    col = 1; #inicilizando variáives dos loops (linha)
    lin = 1; #inicilizando variáives dos loops (coluna)
    ll = 1; #inicilizando variáives dos loops (busca de nós duplicados)
    kk =1; #inicilizando variáives dos loops (busca de nós duplicados)
    Last_node = 1; #inicializando variável que computa o último nó do Loop
    CONDX = Array{Float64}(undef,n*(2*n_aux),2); #inicialização posicao do vetor x
    CONDY = Array{Float64}(undef,n*(2*n_aux),2); #inicialização posicao do vetor y
    CONDZ = Array{Float64}(undef,n*(2*n_aux),2); #inicialização posicao do vetor z
    Nos = Array{Int64}(undef,n*(2*n_aux),2); #inicialização posicao dos nós que conectam cada segmento
    Node = Array{Int64}(undef,2*n_aux,2);
    Condutores = Array{Float64}(undef,n+1,);
    #Posicionando informações da malha na matriz
    #Pos1 = Matriz contendo informações dos eletrodos paralelos
    #Pos2 = Matriz contendo informações dos eletroso ortognais
    for col = 1:(nMalha + 1)
        for lin = 1:nMalha
            Pos1[lin+nMalha*(col-1),1] = (lin-1)*L;
            Pos2[lin+nMalha*(col-1),5] = (lin-1)*L;

            Pos1[lin+nMalha*(col-1),2] = (lin)*L;
            Pos2[lin+nMalha*(col-1),6] = (lin)*L;

            Pos1[lin+nMalha*(col-1),5] = (col-1)*L;
            Pos2[lin+nMalha*(col-1),1] = (col-1)*L;

            Pos1[lin+nMalha*(col-1),6] = (col-1)*L;
            Pos2[lin+nMalha*(col-1),2] = (col-1)*L;

            Pos1[lin+nMalha*(col-1),3] = -h;
            Pos2[lin+nMalha*(col-1),3] = -h;

            Pos1[lin+nMalha*(col-1),4] = -h;
            Pos2[lin+nMalha*(col-1),4] = -h;
        end
    end
    #Início dos eletrodos
    for lin = 1:n_aux
        for col = 1:n_aux
            if Pos1[lin,1] == Pos2[col,1] && Pos1[lin,5] == Pos2[col,5]
                Pos1[lin,7] = col + n_aux;
                Pos1[lin,8] = 1;
            end
        end
    end
    for lin = 1:n_aux
        for col = 1:n_aux
            if Pos1[lin,1] == Pos2[col,2] && Pos1[lin,5] == Pos2[col,6]
                Pos1[lin,7] = col + n_aux;
                Pos1[lin,8] = 2;
            end
        end
    end
    #Fim dos Eletrodos
    for lin = 1:n_aux
        for col = 1:n_aux
            if Pos1[lin,2] == Pos2[col,1] && Pos1[lin,6] == Pos2[col,5]
                Pos1[lin,9] = col + n_aux;
                Pos1[lin,10] = 1;
            end
        end
    end
    for lin = 1:n_aux
        for col = 1:n_aux
            if Pos1[lin,2] == Pos2[col,2] && Pos1[lin,6] == Pos2[col,6]
                Pos1[lin,9] = col + n_aux;
                Pos1[lin,10] = 2;
            end
        end
    end
    #Corrigindo nós "duplicados"
    for lin = 1:n_aux
        for col = 1:n_aux
            if Pos1[lin,2] == Pos1[col,1] && Pos1[lin,6] == Pos1[col,5]
                Pos1[lin,9] = col;
                Pos1[lin,10] = 1;
            end
        end
    end
    for lin = 1:n_aux
        for col = 1:n_aux
            if Pos1[lin,1] == Pos1[col,2] && Pos1[lin,5] == Pos1[col,6]
                Pos1[lin,7] = col;
                Pos1[lin,8] = 2;
            end
        end
    end
    for col = 1:nMalha + 1
        for lin = 1:nMalha
            for ll = 1:n_aux
                for kk = 1:n_aux
                    if Pos1[ll,1] == Pos2[kk,1] && Pos1[ll,5] == Pos2[kk,5]
                        Pos2[kk,7] = ll;
                        Pos2[lin+nMalha*(col-1),8] = 1;
                    end
                    if Pos1[ll,2] == Pos2[kk,2] && Pos1[ll,6] == Pos2[kk,6]
                        Pos2[kk,9] = ll;
                        Pos2[lin+nMalha*(col-1),10] = 2;
                    end
                end
            end
            for ll = 1:n_aux
                for kk = 1:n_aux
                    if Pos1[ll,2] == Pos2[kk,1] && Pos1[ll,6] == Pos2[kk,5]
                        Pos2[kk,7] = ll;
                        Pos2[kk,8] = 2;
                    end
                    if Pos1[ll,1] == Pos2[kk,2] && Pos1[ll,5] == Pos2[kk,6]
                        Pos2[kk,9] = ll;
                        Pos2[kk,10] = 1;
                    end
                end
            end
        end
    end
    Pos = [Pos1' Pos2']'; size(Pos)
    Last_node = 1;
    for lin = 1:2*n_aux
        #global Last_node;
        #Direção x
        Condutores = linspace(Pos[lin,1],Pos[lin,2],n+1);
        CONDX[((lin-1)*n)+1:(lin)*n,:] = [Condutores[1:end-1] Condutores[2:end]];
        #Direção y
        Condutores = linspace(Pos[lin,3],Pos[lin,4],n+1);
        CONDY[((lin-1)*n)+1:(lin)*n,:] = [Condutores[1:end-1] Condutores[2:end]];
        #Direção z
        Condutores = linspace(Pos[lin,5],Pos[lin,6],n+1);
        CONDZ[((lin-1)*n)+1:(lin)*n,:] = [Condutores[1:end-1] Condutores[2:end]];

        #salvando os nós de cada segmentação
        if Pos[lin,7] > lin && Pos[lin,9] > lin #PRIMEIRO ELETRODO
            No = collect(Last_node:1:(Last_node+n));
        end
        if Pos[lin,7] < lin && Pos[lin,9] > lin
            No = collect(Last_node:1:(Last_node+n-1));
            No = [Node[Int(Pos[lin,7]),Int(Pos[lin,8])] No']';
        end
        if Pos[lin,7] > lin && Pos[lin,9] < lin
            No = collect(Last_node:1:(Last_node+n-1));
            No = [No' Node[Int(Pos[lin,9]),Int(Pos[lin,10])]]';
        end
        if Pos[lin,7] < lin && Pos[lin,9] < lin
            No = Last_node:1:(Last_node+n-2);
            No = [Node[Int(Pos[lin,7]),Int(Pos[lin,8])] No' Node[Int(Pos[lin,9]),Int(Pos[lin,10])]];
        end
        Nos[((lin-1)*n)+1:(lin)*n,:] = [No[1:end-1] No[2:end]];
        Node[lin,:] = [No[1] No[end]];
        Last_node = max(No) + 1;
    end
    #Inicializando matrizes de conectividade
    A = zeros(n*(2*n_aux),Last_node-1);
    S = zeros(n*(2*n_aux),Last_node-1);
    for lin = 1:(n*(2*n_aux)) #Numero de eletrodos
        A[lin,Nos[lin,1]] = 0.5;
        S[lin,Nos[lin,1]] = 1.0;
        A[lin,Nos[lin,2]] = 0.5;
        S[lin,Nos[lin,2]] = -1.0;
    end

    #=writedlm("CONDX.csv",CONDX,',');
    writedlm("CONDY.csv",CONDY,',');
    writedlm("CONDZ.csv",CONDZ,',');
    writedlm("Node.csv",Node,',');
    writedlm("A.csv",A,',');
    writedlm("S.csv",S,',');=#
    plot([CONDX[:,1]' CONDX[end,2]]',[CONDZ[:,1]' CONDZ[end,2]]',linestyle =
       :dot,seriestype=:scatter,title="Malha de Aterramento",legend = false, aspect_ratio = 1);
    savefig("Malha_Aterramento.pdf");
    return CONDX, CONDY, CONDZ, Node, A, S
end
function MC_Horizontal(r::Float64,L::Float64,h::Float64,segR::Float64)
    n = floor(Int, L/(segR)); #Número de segmentações de cada "eletrodo"
    n_aux = 1; #número de hastes horizontais (cabos contrapeso)
    Pos = Array{Float64}(undef,n_aux,10); #inicializando matriz auxiliar
    CONDX = Array{Float64}(undef,n*(n_aux),2); #inicialização posicao do vetor x
    CONDY = Array{Float64}(undef,n*(n_aux),2); #inicialização posicao do vetor y
    CONDZ = Array{Float64}(undef,n*(n_aux),2); #inicialização posicao do vetor z
    Nos = Array{Int64}(undef,n*(n_aux),2); #inicialização posicao dos nós que conectam cada segmento
    Node = Array{Int64}(undef,n_aux,2);
    Nos = Array{Int64}(undef,n*(n_aux),2);
    Condutores = Array{Float64}(undef,n+1,);
	Δx = L/n;
	mm = 1; posx = 0.0;
	for linha = 1:n
		CONDY[mm,1] = -h; CONDY[mm,2] = -h;
		CONDZ[mm,1] = 0.0; CONDZ[mm,2] = 0.0;
		CONDX[mm,1] = posx;
		posx = posx + Δx;
		CONDX[mm,2] = posx;
		Nos[mm,1] = mm;
		Nos[mm,2] = mm + 1;
		if linha == 1
			Node[1,1] = Nos[mm,1]
		end
		if linha == n
			Node[1,2] = Nos[mm,2]
		end
		mm = mm + 1;
	end
	Last_node = Nos[end,2] + 1;
    #MATRIZES DE CONECTIVIDADE
    A = zeros(n*(n_aux),Last_node - 1);
    S = zeros(n*(n_aux),Last_node - 1);
    #S = A;
    for linha = 1:length(CONDX[:,1]) #NUMERO DE ELETRODOS
        A[linha,Nos[linha,1]] = 0.5;
        S[linha,Nos[linha,1]] = 1.0;
        A[linha,Nos[linha,2]] = 0.5;
        S[linha,Nos[linha,2]] = -1.0;
    end
    #=writedlm("CONDX.csv",CONDX,',');
    writedlm("CONDY.csv",CONDY,',');
    writedlm("CONDZ.csv",CONDZ,',');
    writedlm("Node.csv",Node,',');
    writedlm("A.csv",A,',');
    writedlm("S.csv",S,',');=#

    plot([CONDX[:,1]' CONDX[end,2]]',[CONDZ[:,1]' CONDZ[end,2]]',linestyle =
        :dot,seriestype=:scatter,title="Grounding Grid",legend = false, aspect_ratio = 1);
    savefig("Grounding_Grid_Design.pdf");
    return CONDX, CONDY, CONDZ, Node, A, S
end
#Function that generate grounding composed of vertical rods only
function MC_Vertical_Rods(r::Float64,haste::Float64,n_haste::Int64,distX::Float64,h::Float64,segR::Float64)
    ny = floor(Int, haste/(segR)); #Número de segmentações de cada "eletrodo"
	n = n_haste*ny;
	CONDX = Array{Float64}(undef,n,2); #inicialização posicao do vetor x
    CONDY = Array{Float64}(undef,n,2); #inicialização posicao do vetor y
    CONDZ = Array{Float64}(undef,n,2); #inicialização posicao do vetor z
    Nos = Array{Int64}(undef,n,2); #inicialização posicao dos nós que conectam cada segmento
    Node = Array{Int64}(undef,n_haste,2);
    Nos = Array{Int64}(undef,n,2);
    Last_node = 1; #inicializando variável que computa o último nó do Loop
	Δy = haste/ny;
	mm = 1;
	for col = 1:n_haste
		posx = (col-1)*distX;
		posy = h;
		posz = 0.0;
		for linha = 1:ny
			CONDY[mm,1] = posy;
			posy = posy + Δy
			CONDY[mm,2] = posy;
			CONDX[mm,1] = posx; CONDX[mm,2] = posx;
			CONDZ[mm,1] = posz;	CONDZ[mm,2] = posz;
			Nos[mm,1] = mm + col - 1;
			Nos[mm,2] = mm + col;
			if linha == 1
				Node[col,1] = Nos[mm,1]
			end
			if linha == ny
				Node[col,2] = Nos[mm,2]
			end
			mm = mm + 1;
		end
	end
	Last_node = Nos[end,2] + 1;
    #MATRIZES DE CONECTIVIDADE
    A = zeros(n,Last_node - 1);
    S = zeros(n,Last_node - 1);

    for linha = 1:length(CONDX[:,1]) #NUMERO DE ELETRODOS
        A[linha,Nos[linha,1]] = 0.5;
        S[linha,Nos[linha,1]] = 1.0;
        A[linha,Nos[linha,2]] = 0.5;
        S[linha,Nos[linha,2]] = -1.0;
    end
    #=writedlm("CONDX.csv",CONDX,',');
    writedlm("CONDY.csv",CONDY,',');
    writedlm("CONDZ.csv",CONDZ,',');
    writedlm("Node.csv",Node,',');
    writedlm("A.csv",A,',');
    writedlm("S.csv",S,',');=#

    plot([CONDX[:,1]' CONDX[end,2]]',-[CONDY[:,1]' CONDY[end,2]]',linestyle =
        :dot,seriestype=:scatter,title="Grounding Grid",legend = false, aspect_ratio = 1);
    savefig("Grounding_Grid_Design.pdf");
    return CONDX, CONDY, CONDZ, Node, A, S
end
#Function that generate the Distribution Line Grounding (CEMIG)
function MC_Dist_Line(r::Float64,n_haste::Int64,h::Float64,segR::Float64)
	haste = 2.5; #Vertical Rod Length
	distX = 3.0; #Horizontal Rod Length
	Lx = (n_haste-1)*distX;
    ny = floor(Int, haste/(segR)); #Número de segmentações de cada "eletrodo"
	nx = floor(Int, Lx/(segR));
	n = nx + n_haste*ny;
	CONDX = Array{Float64}(undef,n,2); #inicialização posicao do vetor x
    CONDY = Array{Float64}(undef,n,2); #inicialização posicao do vetor y
    CONDZ = Array{Float64}(undef,n,2); #inicialização posicao do vetor z
    Nos = Array{Int64}(undef,n,2); #inicialização posicao dos nós que conectam cada segmento
    Node = Array{Int64}(undef,n_haste,2);
    Nos = Array{Int64}(undef,n,2);
    Last_node = 1; #inicializando variável que computa o último nó do Loop
	Δx = Lx/nx;
	Δy = haste/ny;
	mm = 1; posz = 0.0;
	#Vertical Rods
	for col = 1:n_haste
		posx = (col-1)*distX;
		posy = h;
		for linha = 1:ny
			CONDY[mm,1] = posy;
			posy = posy + Δy
			CONDY[mm,2] = posy;
			CONDX[mm,1] = posx; CONDX[mm,2] = posx;
			CONDZ[mm,1] = posz;	CONDZ[mm,2] = posz;
			Nos[mm,1] = mm + col - 1;
			Nos[mm,2] = mm + col;
			if linha == 1
				Node[col,1] = Nos[mm,1]
			end
			if linha == ny
				Node[col,2] = Nos[mm,2]
			end
			mm = mm + 1;
		end
	end
	Last_node = Nos[mm-1,2] + 1;
	#Horizontal Rods
	posx = 0.0;
	Comp = 0.0; mmm = 1;
	for linha = 1:nx
		CONDY[mm,1] = h; CONDY[mm,2] = h;
		CONDX[mm,1] = posx; posx = posx + Δx;
		CONDX[mm,2] = posx;
		CONDZ[mm,1] = posz;	CONDZ[mm,2] = posz;
		if abs(Comp - posx + Δx) < Δx
			Nos[mm,1] = Node[mmm,1]; mmm = mmm + 1;
			Comp = Comp + distX;
			if linha > 1
				Nos[mm - 1, 2] = Nos[mm,1];
			end
		else
			Nos[mm,1] = Last_node; Last_node = Last_node + 1;
		end
		if linha == nx
			Nos[mm,2] = Node[mmm,1];
		else
			Nos[mm,2] = Last_node;
		end
		mm = mm + 1;
	end
	# Last_node = Nos[end,2] + 1;
    #MATRIZES DE CONECTIVIDADE
    A = zeros(n,Last_node - 1);
    S = zeros(n,Last_node - 1);

    for linha = 1:length(CONDX[:,1]) #NUMERO DE ELETRODOS
        A[linha,Nos[linha,1]] = 0.5;
        S[linha,Nos[linha,1]] = 1.0;
        A[linha,Nos[linha,2]] = 0.5;
        S[linha,Nos[linha,2]] = -1.0;
    end
    #=writedlm("CONDX.csv",CONDX,',');
    writedlm("CONDY.csv",CONDY,',');
    writedlm("CONDZ.csv",CONDZ,',');
    writedlm("Node.csv",Node,',');
    writedlm("A.csv",A,',');
    writedlm("S.csv",S,',');=#

    plot([CONDX[:,1]' CONDX[end,2]]',-[CONDY[:,1]' CONDY[end,2]]',linestyle =
        :dot,seriestype=:scatter,title="Grounding Grid",legend = false, aspect_ratio = 1);
    savefig("Grounding_Grid_Design.pdf");
    return CONDX, CONDY, CONDZ, Node, A, S, Nos
end
