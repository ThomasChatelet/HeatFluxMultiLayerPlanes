module HeatFluxMultiLayerPlanes

using Plots
using NumericalIntegration
using MyPhysicalConstants
using OpticalProperties
using QuadGK
using MultiLayerNFRHT
using Revise
using Dates
using DrWatson

#includet("C:\\Users\\tchatelet\\Documents\\Julia\\PlanplanFlux.jl")
srcdir()
#using GR #A installer si pb échelle
plotlyjs()
gr()
#C:\\Users\\tchatelet\\Documents\\Julia\\HeatFluxInfinitePlanesBulk
hbar = 1.0545718e-34
k_b = 1.38064852e-23
Corps_noir = Cst(1.00 + im*1e-5)

@quickactivate "HeatFluxMultiLayerPlanes"

function difftab(lengg,Homemade,Merchfunc,Errorcompare,abscisse3,Errorrelative,abscisse2)
    for j in 1:lengg
                Errorcompare[j] = Homemade[j] - Merchfunc[j]
                Errorrelative[j] = Errorcompare[j]/Merchfunc[j]
                abscisse3[j] = abscisse2[j]
                #println(io, round(abscisse2[j],digits=6),  round(Homemade[j],digits=6),  round(Merchfunc[j],digits=6));
    end
end

function BEdistribdotEnergy(w,T)
    return hbar* w* 1.0/(exp((hbar*w)/(k_b*T))-1)
end

#Definition de kperp
function Kz(omega, permit, kpar)
    return sqrt( (omega*omega*permit)/(c0*c0) - kpar*kpar +0.0im)
end

"""
Définition des coefficients de réflexion et de transmission selon transverse electrique et magnétique
"""
function reflectcoefTE(omega,kpar,Eps1,Eps2) #Patch4mat
    return (( Kz(omega,Eps1,kpar) - Kz(omega,Eps2,kpar)))/( Kz(omega,Eps1,kpar) + Kz(omega,Eps2,kpar))
end

function transmiscoefTE(omega,kpar,Eps1,Eps2) #Patch4mat
    return  (2.0*Kz(omega,Eps1,kpar))/( Kz(omega,Eps1,kpar) + Kz(omega,Eps2,kpar))
end

function transmiscoefTM(omega,kpar,Eps1,Eps2) #Patch4mat
    return ( 2.0*Kz(omega,Eps1,kpar) * sqrt(Eps1+0im) * sqrt(Eps2+0im) +0.0im ) /( Kz(omega,Eps1,kpar)*Eps2 + Kz(omega,Eps2,kpar)*Eps1 )
end

function reflectcoefTM(omega,kpar,Eps1,Eps2)
    return ( Kz(omega,Eps1,kpar)*Eps2 - Kz(omega,Eps2,kpar)*Eps1      )/( Kz(omega,Eps1,kpar)*Eps2 + Kz(omega,Eps2,kpar)*Eps1 )
end


#=
function gencomplextab(leng)
    l = [0.0 + 0im]
    for i in 1:leng-1
        append!(l,0.0 + 0im)
    end
    return l
end

function checkcoefTE(Mat1,Mat2,omega)
    leng = 1250
    kparallele = 10 .^range(3,9, length = leng)
    Eps1 = permittivity(Mat1,omega)
    Eps2 = permittivity(Mat2,omega)
    R_TE = gencomplextab(leng)
    R_TM = gencomplextab(leng)
    T_TE = gencomplextab(leng)
    T_TM = gencomplextab(leng)
    r_TE = zeros(leng)
    r_TM = zeros(leng)
    t_TE = zeros(leng)
    t_TM = zeros(leng)
    consistancyTE = zeros(leng)
    consistancyTM = zeros(leng)
    for i in 1:leng
        R_TE[i] = reflectcoefTE(omega,kparallele[i],Eps1,Eps2)
        R_TM[i] = reflectcoefTM(omega,kparallele[i],Eps1,Eps2)
        T_TE[i] = transmiscoefTE(omega,kparallele[i],Eps1,Eps2)
        T_TM[i] = transmiscoefTM(omega,kparallele[i],Eps1,Eps2)
        r_TE[i] =  abs2(R_TE[i])
        r_TM[i] =  abs2(R_TM[i])
        t_TE[i] =  abs2(T_TE[i])
        t_TM[i] =  abs2(T_TM[i])
        #t_TE[i] =  (real(Kz(omega,kparallele[i],Eps2))*abs2(T_TE[i]))/(real(Kz(omega,kparallele[i],Eps1)))
        #t_TM[i] =  real(Kz(omega,kparallele[i],Eps2)/Eps2)*abs2(T_TM[i])*abs(Eps2/Eps1)/(real(Kz(omega,kparallele[i],Eps1)/Eps1))
        consistancyTE[i] = real(r_TE[i]) + real(t_TE[i])
        consistancyTM[i] = real(r_TM[i]) + real(t_TM[i])
        # real(k2z)/real(k0z)*abs(t)^2
    end
    plot(kparallele,consistancyTE,xlabel = "kparallele",  label = "consistancyTE", xaxis=:log, ylim=(-0.01, 1.25))
    plot!(kparallele,r_TE,xlabel = "kparallele",  label = "R_TE", xaxis=:log)
    plot!(kparallele,t_TE,xlabel = "kparallele",  label = "T_TE", xaxis=:log, ylim=(-0.01, 1.35) )
    #=plot(kparallele,consistancyTM,xlabel = "kparallele",  label = "consistancyTM", xaxis=:log)
    plot!(kparallele,R_TM,xlabel = "kparallele",  label = "R_TM", xaxis=:log)
    plot!(kparallele,T_TM,xlabel = "kparallele",  label = "T_TM", xaxis=:log)=#
    #print(T_TE)
end

omega = 1.00e12
kpar = omega/(2*c0)

checkcoefTE(Layer(Au()),Layer(Al()),omega)

tabMat1 = [Corps_noir, Layer(Au()),Layer(Sic()),Layer(Vacuum),Layer((Ti)),Layer(TiW_v2),Layer(Si_n_doped(nSi_masetti,1e17)),Layer(Si_p_doped(pSi_masetti,1e17))]
tabMat2 = [Corps_noir, Layer(Au()),Layer(Sic()),Layer(Vacuum),Layer((Ti)),Layer(TiW_v2),Layer(Si_n_doped(nSi_masetti,1e17)),Layer(Si_p_doped(pSi_masetti,1e17))]

function masscheckcoefTE(tabMat1,tabMat2,omega)
    size1, = size(tabMat1)
    size2, = size(tabMat2)
    cd(datadir())
    cd("exp_raw")
    datetoday = formatdat()
    gendir = datetoday * " reflectcoef mass check "
    mkdir(gendir)
    cd(gendir)
    for i in 1:size1
        Mat1 = tabMat1[i]
        for j in 1:size2
            Mat2 = tabMat2[j]
            Mat1str = identifymat(Mat1)
            Mat2str = identifymat(Mat2)
            checkcoefTE(Mat1,Mat2,omega)
            #rm("Hflux" * Mat1str * " " * Mat3str)
            savefig("Coefcheck" * Mat1str * " " * Mat2str *" .png")
        end
    end
end

#tt = transmiscoefTE(omega,kpar,permittivity(Layer(TiW_v2),omega),1)
#rr = reflectcoefTE(omega,kpar,permittivity(Layer(TiW_v2),omega),1)

#masscheckcoefTE(tabMat1,tabMat2,omega)
=#
"""
Intégrande des differentes parties de k parallèle
"""

#=
function integpropagTEMultilayer(omega,kpar,d,Eps2,tablayer,tablayer2,tablayer2,tablaysize,tabepaisseur)
    #Eps3 = permittivity(tabLayer[tablaysize-1],omega)
    (r_TE1,b,r_TM1,d,e,f,sE,sM) = calculreflexionmulticouches(omega,kpar,tablayer,tablaysize,tabepaisseur)
    tempsize, = size(tablayer2)
    (r_TE3,b,r_TM3,d,e,f,sE,sM) = calculreflexionmulticouches(omega,kpar,tablayer2,tempsize,tabepaisseur)
    return kpar*((((1 - abs2(r_TE1)) * (1 - abs2(r_TE3))))/(abs2(1 - (r_TE1 * r_TE3 *exp(2im*d*Kz(omega, Eps2, kpar))))))
end

function integpropagTMMultilayer(omega,kpar,d,Eps2,tablayer,tablayer2,tablaysize,tabepaisseur)
    tempsize, = size(tablayer2)
    (r_TE3,b,r_TM3,d,e,f,sE,sM) = calculreflexionmulticouches(omega,kpar,tablayer2,tempsize,tabepaisseur)
    (r_TE1,b,r_TM1,d,e,f,sE,sM) = calculreflexionmulticouches(omega,kpar,tablayer,tablaysize,tabepaisseur)
    return kpar*(((1 - abs2(r_TM1)) * (1 - abs2(r_TM3)))/(abs2(1 - (r_TM1 * r_TM3 * exp(2im*d*Kz(omega, Eps2, kpar))))))
end

function propagshape(r1,r3,d,kpar,w)
    return kpar*((((1 - abs2(r1)) * (1 - abs2(r3))))/(abs2(1 - (r1 * r3 * exp(2im*d*Kz(omega, Eps2, kpar))))))
end
=#
function integpropagTEMultilayer(omega,kpar,d,Eps2,tablayer,tablayer2)
    #Eps3 = permittivity(tabLayer[tablaysize-1],omega)
    (r_TE1,r_TM1,sE,sM) = calculreflexionmulticouches(omega,kpar,tablayer)
    (r_TE3,r_TM3,sE,sM) = calculreflexionmulticouches(omega,kpar,tablayer2)
    return kpar*((((1 - abs2(r_TE1)) * (1 - abs2(r_TE3))))/(abs2(1 - (r_TE1 * r_TE3 *exp(2im*d*Kz(omega, Eps2, kpar))))))
end

function integpropagTMMultilayer(omega,kpar,d,Eps2,tablayer,tablayer2)
    (r_TE1,r_TM1,sE,sM) = calculreflexionmulticouches(omega,kpar,tablayer)
    (r_TE3,r_TM3,sE,sM) = calculreflexionmulticouches(omega,kpar,tablayer2)
    return kpar*(((1 - abs2(r_TM1)) * (1 - abs2(r_TM3)))/(abs2(1 - (r_TM1 * r_TM3 * exp(2im*d*Kz(omega, Eps2, kpar))))))
end



function integevaTEMultilayer(omega,kpar,d,Eps2,tablayer,tablayer2)
    (r_TE1,r_TM1,sE,sM) = calculreflexionmulticouches(omega,kpar,tablayer)
    (r_TE3,r_TM3,sE,sM) = calculreflexionmulticouches(omega,kpar,tablayer2)
    k2perp = Kz(omega, Eps2, kpar)
    return 4*kpar*exp(-2*d*imag(k2perp))*((imag(r_TE1)) * (imag(r_TE3)))/(abs2(1-(r_TE1 * r_TE3 * exp(-2*d*imag(k2perp)))))
end

function integevaTMMultilayer(omega,kpar,d,Eps2,tablayer,tablayer2)
    (r_TE1,r_TM1,sE,sM) = calculreflexionmulticouches(omega,kpar,tablayer)
    (r_TE3,r_TM3,sE,sM) = calculreflexionmulticouches(omega,kpar,tablayer2)
    k2perp = Kz(omega, Eps2, kpar)
    return (4*kpar*exp(-2*d*imag(k2perp)))*(((imag(r_TM1)) * (imag(r_TM3)))/(abs2(1-(r_TM1 * r_TM3 * exp(-2*d*imag(k2perp))))))
end


"""
Intégration sur Kparallele
"""
function FpropagTE(omega,d,T1,T3,tablayer,tablayer2)
    bornedecoup = omega/c0
    Eps2 = 1.0
    prefact = (1.0/(4*pi*pi))*(BEdistribdotEnergy(omega,T1) - BEdistribdotEnergy(omega,T3))
    valTE, err3 = quadgk( kpar -> integpropagTEMultilayer(omega,kpar,d,Eps2,tablayer,tablayer2), 0, bornedecoup, rtol=1e-6)
    return prefact*real(valTE) #catch
end

function FpropagTM(omega,d,T1,T3,tablayer,tablayer2)
    bornedecoup = omega/c0
    Eps2 = 1.0 + 0.0im
    prefact = (1.0/(4*pi*pi))*(BEdistribdotEnergy(omega,T1) - BEdistribdotEnergy(omega,T3))
    bornesup = Inf
    valTM, err5 = quadgk( kpar -> integpropagTMMultilayer(omega,kpar,d,Eps2,tablayer,tablayer2), 0, bornedecoup, rtol=1e-6)
    return prefact*(real(valTM))
end

function FevaTE(omega,d,T1,T3,tablayer,tablayer2)
    bornedecoup = omega/c0
    Eps2 = 1.0
    prefact = (1.0/(4*pi*pi))*(BEdistribdotEnergy(omega,T1) - BEdistribdotEnergy(omega,T3))
    valTE2, err4 = quadgk( kpar -> integevaTEMultilayer(omega,kpar,d,Eps2,tablayer,tablayer2), bornedecoup, Inf, rtol=1e-6)
    return prefact*( real(valTE2) )
end

function FevaTM(omega,d,T1,T3,tablayer,tablayer2)
    bornedecoup = omega/c0
    Eps2 = 1.0+ 0.0im
    prefact = (1.0/(4*pi*pi))*(BEdistribdotEnergy(omega,T1) - BEdistribdotEnergy(omega,T3))
    valTM2, err6 = quadgk( kpar -> integevaTMMultilayer(omega,kpar,d,Eps2,tablayer,tablayer2), bornedecoup, Inf, rtol=1e-6)
    return prefact*( real(valTM2) )
end

function generates0()
    s0 = [0.0 + 0im,0.0 + 0im,0.0 + 0im,0.0 + 0im]
    s0[1] = 1.0 + 0im
    s0[4] = 1.0 + 0im
    #s0 = Complex(s0)
    return s0
end

s0 = generates0()

function evolves0TE(S_TE,Mat1,Mat2,omega,kparallele,epaisseur)
    #dl = epaisseur de Mat1
    S11 = S_TE[1]
    S12 = S_TE[2]
    S21 = S_TE[3]
    S22 = S_TE[4]
    Eps1 = permittivity(Mat1,omega)
    Eps2 = permittivity(Mat2,omega)
    exp1kz = exp(1im*epaisseur*Kz(omega,Eps1,kparallele))
    exp2kz = exp(2im*epaisseur*Kz(omega,Eps1,kparallele))
    r_TE = reflectcoefTE(omega,kparallele,Eps1,Eps2)
    t_TE = transmiscoefTE(omega,kparallele,Eps1,Eps2)
    S_TE[1] = (S11 * t_TE * exp1kz)   /(1.0-S12*r_TE*exp2kz)
    S_TE[2] = (S12*exp2kz - r_TE)     /(1.0-S12*r_TE*exp2kz)
    S_TE[3] = S21 + (1.0/t_TE)*(S_TE[1]*S22*r_TE*exp1kz)
    S_TE[4] = (1.0/t_TE)*(S22*(r_TE*S_TE[2]+1.0+0im)*exp1kz)
end

function evolves0TM(S_TM,Mat1,Mat2,omega,kparallele,epaisseur)
    #dl = epaisseur de Mat1
    S11 = S_TM[1]
    S12 = S_TM[2]
    S21 = S_TM[3]
    S22 = S_TM[4]
    Eps1 = permittivity(Mat1,omega)
    Eps2 = permittivity(Mat2,omega)
    exp1kz = exp(1im*epaisseur*Kz(omega,Eps1,kparallele))
    exp2kz = exp(2im*epaisseur*Kz(omega,Eps1,kparallele))
    r_TM = reflectcoefTM(omega,kparallele,Eps1,Eps2)
    t_TM = transmiscoefTM(omega,kparallele,Eps1,Eps2)
    S_TM[1] = (S11 * t_TM * exp1kz)   /(1.0-S12*r_TM*exp2kz)
    S_TM[2] = (S12*exp2kz - r_TM)     /(1.0-S12*r_TM*exp2kz)
    S_TM[3] = S21 + (1.0/t_TM)*(S_TM[1]*S22*r_TM*exp1kz)
    S_TM[4] = (1.0/t_TM)*(S22*(r_TM*S_TM[2]+1.0+0im)*exp1kz)
end


s0 = generates0()
tabLayer = generatelayertab(1)

#=
omega = 1e14
kpar = omega/(2.0*c0)
s0 = generates0()
Current_mat = Layer(Al())
evolves0TE(s0,Current_mat,Current_mat,omega,kpar,0)
s0
tt = transmiscoefTE(omega,kpar,permittivity(Current_mat,omega),permittivity(Current_mat,omega))
rr = reflectcoefTE(omega,kpar,permittivity(Current_mat,omega),permittivity(Current_mat,omega))
abs(rr) - abs(s0[3])

abs(tt) + abs(rr)
real(tt) + real(rr)
=#

function generatethicknesstab(leng) #Generate random thicknesses for Layers  10-8 => 10-4
    tabthick = zeros(leng)
    for i in 1:leng
        tabthick[i]= 1*10^(-1*(rand(5:7) + rand(Float64)))
    end
    return tabthick
end

function generatethicknesstab2(leng) #Generate random thicknesses for Layers  10-8 => 10-4
    tabthick = zeros(leng)
    for i in 1:leng
        tabthick[i]= 1*10^(-7)
    end
    return tabthick
end

function generatelayertab(leng)
    tablay = [Layer(Vacuum)]
    pop!(tablay) #Type issue
    tabMat1 = [(Corps_noir),(Al()),(Au()),(Sic()),((Ti)),(TiW_v2),(Si_n_doped(nSi_masetti,1e17)),(Si_p_doped(pSi_masetti,1e17)),(Vacuum)]
    for i in 1:leng
        draw = rand(1:1)
        objet = tabMat1[draw]
        push!(tablay,Layer(objet,1*10^(-1*(rand(5:7) + rand(Float64)))))
    end
    return tablay
end

function generatelayertab2(leng,key,dimthick)
    tablay = [Layer(Vacuum)]
    #pop!(tablay)
    tabMat1 = [(Corps_noir),(Al()),(Au()),(Sic()),((Ti)),(TiW_v2),(Si_n_doped(nSi_masetti,1e17)),(Si_p_doped(pSi_masetti,1e17)),(Vacuum)]
    for i in 1:leng
        draw = rand(key:key)
        objet = tabMat1[draw]
        #push!(tablay,Layer(objet,1*10^(-1.0*dimthick)))
        push!(tablay,Layer(objet,dimthick))
    end
    #push!(tablay,Layer(Vacuum))
    return tablay
end

function calculreflexionmulticouches(omega,kparallele,tabllayer)
    S_TransverseTE = generates0()
    S_TransverseTM = generates0()
    nbcouches, = size(tabllayer)
    rTE_end = 0
    rTM_end = 0
    #=
    if (nbcouches != size(tabLayer))
        println("Nombre d'épaisseurs != Nombre de layers")
        return
    end
    =#
    #inner = nbcouches - 1
    for i in 2:nbcouches
        evolves0TE(S_TransverseTE,tabllayer[i-1],tabllayer[i],omega,kparallele,tabllayer[i-1].thickness)
        evolves0TM(S_TransverseTM,tabllayer[i-1],tabllayer[i],omega,kparallele,tabllayer[i-1].thickness)
    end
    rTE_end = S_TransverseTE[3]
    rTM_end = S_TransverseTM[3]
    tTE_end = S_TransverseTE[1]
    tTM_end = S_TransverseTM[1]
    return (rTE_end, rTM_end, S_TransverseTE,S_TransverseTM)
end

end
#=
rt(laytab,te(),kpar,omega)
omega = 1e14
kpar = omega/(2.0*c0)
(a,b,c,d,e,f,sE,sM) = calculreflexionmulticouches(omega,kpar,tabLayer)
(a,b,sE,sM) = calculreflexionmulticouches(omega,kpar,tabLayer)
Current_mat =  Layer(Si_n_doped(nSi_masetti,1e17))
tt = transmiscoefTE(omega,kpar,permittivity(Current_mat,omega),1.0 + 0im)
rr = reflectcoefTE(omega,kpar,permittivity(Current_mat,omega),1.0 + 0im)
rrM = reflectcoefTM(omega,kpar,permittivity(Current_mat,omega),1.0 + 0im)
ttM = transmiscoefTM(omega,kpar,permittivity(Current_mat,omega),1.0 + 0im)

abs(rr) - abs(sE[3])
=#
"""
Formattage des output
"""

function formatdat()
    Momentexec = string(Dates.now())
    res = ""
    for i in 1:10
        res*=Momentexec[i]
    end
    res*=" "
    for i in 12:19
        if (Momentexec[i] != ':')
            res*=Momentexec[i]
        end
    end
    return res
end

function identifymat(Mat)
    if (Mat == Layer(Al())) return "Al" end
    if (Mat == Layer(Sic())) return "CarbureSilicium" end
    if (Mat == Layer(Corps_noir)) return "Corps noir" end
    if (Mat == Layer(Vacuum)) return "Vacuum" end
    if (Mat == Layer(Au())) return "Au" end
    if (Mat == Layer(Ti)) return "Ti" end
    if (Mat == Layer(TiW_v2)) return "TiW_v2" end
    if (Mat == Layer(Si_n_doped(nSi_masetti,1e17))) return "Si_n_doped(nSi_masetti,1e17)" end
    if (Mat == Layer(Si_p_doped(pSi_masetti,1e17))) return "Si_p_doped(pSi_masetti,1e17)" end
    return "Inc"
end

"""
Intégration sur omega
"""

function FluxnetechangeMultilayer(d,T1,T3,tablayer,tablayer2; w1 = 1e10, w2 = 1e16) #UndefVarError: FluxnetechangeMultilayer not defined
    #return quadgk(integrandtot,0,Inf,rtol=1e-6)
    valTE,err20 = quadgk( omega -> FpropagTE(omega,d,T1,T3,tablayer,tablayer2), w1, w2; rtol=1e-6)
    valTE2,err21 = quadgk( omega -> FevaTE(omega,d,T1,T3,tablayer,tablayer2), w1, w2; rtol=1e-6)
    valTM,err22 = quadgk( omega -> FpropagTM(omega,d,T1,T3,tablayer,tablayer2), w1, w2; rtol=1e-6)
    valTM2,err23 = quadgk( omega -> FevaTM(omega,d,T1,T3,tablayer,tablayer2), w1, w2; rtol=1e-6)
    return (valTE, valTE2, valTM, valTM2, valTE + valTE2, valTM + valTM2, valTE + valTE2 + valTM + valTM2 )
end



#FluxnetechangeMultilayer(1e-3,1,0,tablayer)

function multilayertest(nbcouches)
    tablayer = generatelayertab(nbcouches)
    pushfirst!(tablayer,tablayer2, Layer(Vacuum))
    tabepaisseur = generatethicknesstab(nbcouches)
    pushfirst!(tabepaisseur, 0.99999999)
    FluxnetechangeMultilayer(1e-3,1,0,tablayer,tablayer2,nbcouches+1,tabepaisseur)
end

#multilayertest(20)

function CalcFlux(tablayer,tablayer2,T1,T3,wbot,wtop,message)
    lengg = 25
    println("Starting " * message)
    datetoday = formatdat()
    abscisse2 = 10 .^range(-9,-5, length = lengg)
    Homemade = zeros(lengg)
    Merchfunc = zeros(lengg)
    TE_contrib = zeros(lengg)
    TM_contrib = zeros(lengg)
    valTE = zeros(lengg)
    valTM = zeros(lengg)
    PropagTEtab = zeros(lengg)
    PropagTMtab = zeros(lengg)
    EvaTEtab = zeros(lengg)
    EvaTMtab = zeros(lengg)
    Mat1str = message #identifymat(Mat1)
    Mat3str = " " #identifymat(Mat3)
    cd(datadir())
    cd("exp_raw")
    #rm("Hflux" * Mat1str * " " * Mat3str)
    gendir = datetoday * " Hflux" * Mat1str * " " * Mat3str * "Heat transfer @ " * string(T1) * "K " * string(T3) * "K "
    mkdir(gendir)
    cd(gendir)
    ioA = open("Abscisse .txt", "w+");
    ioH = open("Calcul HFlux Homemade " * Mat1str * " " * Mat3str * ".txt", "w+");
    ioM = open("Calcul HFlux Merch " * Mat1str * " " * Mat3str * " .txt", "w+");
    ioDiff = open("Calcul Flux Diff Homemade-Merch " * Mat1str * " " * Mat3str* " .txt", "w+");
    for i in 1:lengg #Abscisse2 itération sur d
                #Merchfunc[i] = total_heat_transfer(tablayer, tablayer2, Layer(Vacuum, abscisse2[i]), T1 , T3, 1e10, 1e16 ;tolkx=1e-6,tolw=1e-6)[5]
                if (i%3 == 0) println(repr(i) * " Precalcul Homemade" * Mat1str * " " * Mat3str) end
                (PropagTEtab[i], EvaTEtab[i], PropagTMtab[i],EvaTMtab[i], valTE[i], valTM[i], Homemade[i] )= FluxnetechangeMultilayer(abscisse2[i],T1,T3,tablayer,tablayer2)
                #return (valTE, valTE2, valTM, valTM2, valTE + valTE2, valTM + valTM2, valTE + valTE2 + valTM + valTM2 )
                if (i%3 == 0) println(i) end
                println(ioA, abscisse2[i])
                println(ioH, Homemade[i])
                println(ioM, Merchfunc[i])
    end
    close(ioA)
    close(ioH)
    close(ioM)
    for i in 1:lengg
        (PropagTEtab[i], EvaTEtab[i], PropagTMtab[i],EvaTMtab[i], valTE[i], valTM[i], Homemade[i] ) = (abs(PropagTEtab[i]), abs(EvaTEtab[i]),abs( PropagTMtab[i]),abs(EvaTMtab[i]), abs(valTE[i]), abs(valTM[i]), abs(Homemade[i]) )
    end
    plot(abscisse2,Homemade,xlabel = "Plane separation distance",  label = "Homemade logscale", xaxis=:log, yaxis=:log)
    #plot!(abscisse2,Merchfunc,xlabel = "Plane separation distance", label = "Merchfunc logscale", xaxis=:log, yaxis=:log)
    savefig("HeatFluxSemilogged" * Mat1str * " " * Mat3str * " " * string(T1) * "K " * string(T3) * "K " * " .png")
    plot(abscisse2,Homemade,xlabel = "Plane separation distance",  label = "Homemade logscale", xaxis=:log, yaxis=:log,ylim=(1e-1, 1e8))
    plot!(abscisse2,PropagTEtab,xlabel = "Plane separation distance", label = "PropagTEtab contribution", xaxis=:log, yaxis=:log,ylim=(1e-1, 1e8))
    plot!(abscisse2,EvaTEtab,xlabel = "Plane separation distance", label = "EvaTEtab contribution", xaxis=:log, yaxis=:log,ylim=(1e-1, 1e8))
    plot!(abscisse2,PropagTMtab,xlabel = "Plane separation distance", label = "PropagTMtab contribution", xaxis=:log, yaxis=:log,ylim=(1e-1, 1e8))
    plot!(abscisse2,EvaTMtab,xlabel = "Plane separation distance", label = "EvaTMtab contribution", xaxis=:log, yaxis=:log,ylim=(1e-1, 1e8))
    savefig("HeatFluxContributionslogged" * Mat1str * " " * Mat3str * " " * string(T1) * "K " * string(T3) * "K " * " .png")
    #savefig("C:\\Users\\tchatelet\\Documents\\Julia\\Results-plot", "HeatFluxContributionslogged" * Mat1str * " " * Mat3str * " " * string(T1) * "K " * string(T3) * "K " * " .png")
    lengg = 10
    #=
    Errorcompare = zeros(lengg)
    Errorrelative = zeros(lengg)
    abscisse3 = zeros(lengg)
    difftab(lengg,Homemade,Merchfunc,Errorcompare,abscisse3,Errorrelative,abscisse2)
    plot(abscisse3,Errorrelative, xlabel = "Plane separation distance",label ="Relative error Hmade _ TotalHeatTransferFlux",Title = "Heat transfer @ 300K 600K")
    savefig("Relative_error " * Mat1str * " " * Mat3str * ".png")
    close(ioDiff)
    =#
    println("endof " * Mat1str * " " * Mat3str)
end



function massgenerateplots(Tc,Tf)
    #Tc = 100
    #Tf = 700
    tablayerAl = generatelayertab2(1,2,9e-1)
    tablayerAl2 = generatelayertab2(1,2,9e-1)
    tablayerSic = generatelayertab2(1,4,9e-1)
    tablayerSic2 = generatelayertab2(1,4,9e-1)
    tablayerAu = generatelayertab2(1,3,9e-1)
    tablayerAu2 = generatelayertab2(1,3,9e-1)
    while (Tc <= Tf)
        CalcFlux(tablayerAl,tablayerAl2, Tc ,Tf, 1e8, 1e15, "1 Layer of Al - 1 Layer of Al thick = 1m")
        CalcFlux(tablayerAl, tablayerSic, Tc ,Tf, 1e8, 1e15, "1 Layer of Al - 1 Layer of SiC thick = 1m")
        CalcFlux(tablayerSic,tablayerSic2, Tc ,Tf , 1e8, 1e15, "1 Layer of SiC - 1 Layer of SiC thick = 1m")
        CalcFlux(tablayerAu,tablayerAu2, Tc ,Tf , 1e8, 1e15, "1 Layer of Au - 1 Layer of Au thick = 1m")
        Tf +=50
    end
end

massgenerateplots(100,700)

#https://sci-hub.se/10.1364/OE.20.001903

CalcFlux(Layer(Al()),Layer(Al()),500,300, 1e10, 1e16)
CalcFlux(Layer(Al()),Layer(Al()),1,0, 1e10, 1e16)

CalcFlux(500,300, 1e10, 1e16)

CalcFlux(500,300, 1e10, 1e16)
CalcFlux(500,300, 1e10, 1e16)
CalcFlux(500,300, 1e10, 1e16)


ioA = open("Abscisse .txt", "w+");
ioH = open("Calcul HFlux Homemade .txt", "w+");
ioM = open("Calcul HFlux Merch .txt", "w+");
ioDiff = open("Calcul Flux Diff Homemade-Merch .txt", "w+");

#Test

total_heat_transfer(Layer(Al()), Layer(Al()), Layer(Vacuum, 1e-6), 500 , 300, 5e10, 1e16 ; tolkx=1e-6,tolw=1e-6)[5]
Fluxnetechange(1e-6,500,300,Layer(Al()),Layer(Al()))

total_heat_transfer(Layer(Corps_noir),Layer(Corps_noir) , Layer(Vacuum, 1e-9), 1 , 0, 1e8, 1e15 ; tolkx=1e-6,tolw=1e-6)[5]
FluxnetechangeMultilayer(1e-6,1,0,Corps_noir,Corps_noir)

FluxnetechangeMultilayer(1e-7,1,0,tablayer)

tablayer

plot(abscisse2,Homemade,xlabel = "Plane separation distance",  label = "Homemade logscale", xaxis=:log)
plot!(abscisse2,Merchfunc,xlabel = "Plane separation distance", label = "Merchfunc logscale", xaxis=:log)
savefig("C:\\Users\\tchatelet\\Documents\\Julia\\HeatFluxSemilogged.png")

end
tablayer = generatelayertab2(1,1)
tablayer2 =  generatelayertab2(1,1)
FluxnetechangeMultilayer(1e-9,1,0,tablayer,tablayer2)


Tc = 100
Tf = 700
tablayerAl = generatelayertab2(5,2)
tablayerAl2 = generatelayertab2(5,2)
tablayerSic = generatelayertab2(5,4)
tablayerSic2 = generatelayertab2(5,4)
tablayerAu = generatelayertab2(5,3)
tablayerAu2 = generatelayertab2(5,3)
while (Tc <= Tf)
    #CalcFlux(tablayerAl, tablayerSic, Tc ,Tf, 1e8, 1e15, "Layer of Al-Sic")
    CalcFlux(tablayerSic,tablayerSic2, Tc ,Tf , 1e8, 1e15, "Layer of Sic-Sic")
    CalcFlux(tablayerAu,tablayerAu2, Tc ,Tf , 1e8, 1e15, "Layer of Au-Au")
    #CalcFlux(tablayerAl,tablayerAl2, Tc ,Tf, 1e8, 1e15, "Layer of Al-Al")
    Tc +=50
end
#=
function evolve_amplitude_TE(Al,Bl,S_TE)
    BlMultilayer1= 1.0/(S_TE[4]) * (Bl - S_TE[3])
    AlMultilayer1= S_TE[1] + S_TE[2]*Bl
    return (AlMultilayer1,BlMultilayer1)
end

function evolve_amplitude_TM(Al,Bl,S_TM)
    BlMultilayer1= 1.0/(S_TM[4]) * (Bl - S_TM[3])
    AlMultilayer1= S_TM[1] + S_TM[2]*Bl
    return (AlMultilayer1,BlMultilayer1)
end

tabtest = ones(8)
typeof(size(tabtest))
a, = size(tabtest)

=#
