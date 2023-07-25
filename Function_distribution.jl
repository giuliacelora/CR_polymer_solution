module Module_distribution
export compute_homogeneous_equilibria
using BifurcationKit, LinearAlgebra, Plots, Parameters, Setfield

################## implicit definition homogeneous equilibrium ###########################
function Fun_equilibrium(X,p)
    @unpack α,η,Nmono,ϕP,ϕCL,χ,Nz= p
    
    ϕps=10^X[1] # estimate value of ϕ+/ϕs
    A= 0 # normalising factor
    Q=0 # mean charge
    for z in 0:Nz    
        temp=ϕps.^z.*binomial(Nz,z).*exp(-(α+χ .*ϕP).*z-η.*z.^2 ./(2.0.*Nz)) 
        Q+=temp.*z
        A+=temp
    end  
    Q*=(1/A) 
    ϕs= (1-ϕP-ϕCL)/(1+ϕps)
    electro_neutrality=ϕps*ϕs-ϕCL.+1. /Nmono.*Q*ϕP
    return [electro_neutrality]
end

################## function to extract moments distribution ###########################
function mean_charge(X,ϕP,p)
    """
    return the mean and standard deviation of the charge distribution
    """
    @unpack η,Nmono,χ,Nz,α,ϕCL= p
    ϕps=10^X[1]
    A = 0
    Q =0 
    SECOND_MOMENT =0
    for z in 0:Nz    
        temp=ϕps.^z.*binomial(Nz,z).*exp(-(α+χ .*ϕP).*z-η.*z.^2 ./(2.0.*Nz)) 
        Q+=temp.*z
        SECOND_MOMENT+=temp.*z^2
        A+=temp
    end  
    Q*=(1/A)
    SECOND_MOMENT*=(1/A)
    S=sqrt(SECOND_MOMENT-Q^2)
    return Q,S
end

function compute_homogeneous_equilibria(par_mod)
    sol=Dict()
    sol0=[0.25]
    
    ############ use continuation to solve equilibrium problem ################
    print("solving non-linear problem\n\n")
    for x0 in range(0.005,0.495,step=0.005)
        par_mod =@set par_mod.ϕCL=x0

        Plim1=1-x0 
        Plim2=(2*x0-1)/(par_mod.Nz/par_mod.Nmono-1)

        Plim=minimum([Plim1,Plim2]) # delimit range of existence solution

        opt_newton = NewtonPar(tol = 1e-12, verbose =false,maxIter=50)
        prob = BifurcationProblem(Fun_equilibrium,sol0,par_mod,(@lens _.ϕP),recordFromSolution= (x, p) -> mean_charge(x,p,par_mod)[1];);
        out= newton(prob, opt_newton);
        prob = BifurcationKit.reMake(prob, u0 =out.u);

        
        opt_newton = NewtonPar(tol = 1e-12, verbose =false, maxIter = 3);
        opts_br = ContinuationPar(dsmin=0.0001,dsmax=0.01,ds=0.0075, pMin=0.0001,pMax=Plim,detectBifurcation=1,
           newtonOptions = opt_newton,saveSolEveryStep=1,maxSteps=500);
        branch=BifurcationKit.continuation(prob,PALC(tangent=Bordered()),opts_br;plot=false,verbosity=0,bothside = false);
        sol[x0]=branch;
    end
    ############ processing output in the right format ################
    print("preparing output\n\n")

    xϕP=Vector{Float64}() # volume fraction polymer
    xϕL=Vector{Float64}() # volume fraction counter-ions
    xH=Vector{Float64}()  # volume fraction co-ions (H+)
    meanQ=Vector{Float64}()  # polymer mean charge, Q
    stdQ=Vector{Float64}()   # polymer charge standard deviation, S
    for x0 in range(0.005,0.495,step=0.005)
        for el in sol[x0].sol
            Z=10^el.x[1]
            D=1-el.p-x0

            append!(xH,Z*D/(1+Z))
            append!(xϕP,el.p)
            append!(xϕL,x0)
            temp=mean_charge(el.x,el.p,par_mod)
            append!(meanQ,temp[1])
            append!(stdQ,temp[2])
        end
    end
    M=hcat(xϕP,xϕL,xH,meanQ,stdQ) 
    print("done\n\n")
    ################################################
    return M
    end
end