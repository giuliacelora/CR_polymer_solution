module Module_phase_diagram_FC
export compute_binodals
using BifurcationKit, LinearAlgebra, Plots, Parameters, Setfield,DelimitedFiles,Revise
includet("Functions_Binodal_FC.jl")
import ..FC_model
################## coexisting conditions FC model ###########################
function Tangent_construction(X,p)
    @unpack P = p

    f=similar(X) # output vector
    
    vars=FC_model.transform_variable(X,P,p) # transform into volume fractions

    # phase I: dilute phase variables
    var_I=vars[1,:]
    ϕP_I=var_I[1]
    ϕCL_I=var_I[2]       
    chemical_potentials_I=FC_model.gradient_psi(var_I,p)
    ψ_I=FC_model.psi(var_I,p) # effective free energy phase I

    # phase II: condensed phase variables
    var_II=vars[2,:]
    ϕP_II=var_II[1]
    ϕCL_II=var_II[2]
    chemical_potentials_II=FC_model.gradient_psi(var_II,p)
    ψ_II=FC_model.psi(var_II,p) # effective free energy phase II

    μP=chemical_potentials_I[1]
    μCL=chemical_potentials_I[2]
    
    f[1]=μP-chemical_potentials_II[1] # continuity effective chemical potential polymer
    f[2]=μCL-chemical_potentials_II[2] # continuity effective chemical potential counter-ions

    f[3]=(ψ_II-ψ_I-μP*(ϕP_II-ϕP_I)-μCL*(ϕCL_II-ϕCL_I))/(ϕP_II-ϕP_I) # continuity effective osmotic pressure

    return f
end


function compute_binodal(fixpar)
    par_mod=(Nz=fixpar.Nz,χ=fixpar.χ,P=log10(0.1),λ=fixpar.λ,Nmono=fixpar.Nmono);
    
    ############# load the initial guess #############
    name="initial_condition_FC_chi_0.95"

    M=readdlm(name*".txt", ',', Float64, '\n');
    indeces=sortperm((M[1:end,1].-par_mod.Nz).^2);
    idx=[indeces[1]]
    n=1
    for el in indeces[2:20]
        X=minimum(abs.(el.-indeces[1:n]))
        if X>1 
            append!(idx,el)
        end
        n+=1
    end
   
    ############ 
    branches=Dict()
    sol0=zeros(3,) 
    n=0
    ############ use continuation to solve equilibrium problem ################
    print("computing binodal\n\n")
    for x0 in idx[1]
        sol0=M[x0,2:end]
        opt_newton = NewtonPar(tol = 1e-11, verbose =false,maxIter=50)
        prob = BifurcationProblem(Tangent_construction,sol0,par_mod,(@lens _.P),recordFromSolution= (x, p) -> x[1];);
        out= newton(prob, opt_newton);
        prob = BifurcationKit.reMake(prob, u0 =out.u);
        opt_newton = NewtonPar(tol = 1e-8, verbose =false, maxIter = 3)
        opts_br = ContinuationPar(dsmin=0.00001,dsmax=0.01,ds=-0.005,pMin=-10.0,pMax=0.,detectBifurcation=0,
        newtonOptions = opt_newton,saveSolEveryStep=1,maxSteps=2000)
        branches[n]=BifurcationKit.continuation(prob,PALC(tangent=Bordered()),opts_br;plot=false,bothside = true);
        n+=1
    end
    ############ processing output in the right format ################
    print("preparing output\n\n")

    output=[]
    for n in 0:length(branches)-1

        for el in branches[n].sol

            vars=FC_model.transform_variable(el.x,el.p,par_mod) 
            var_I=vars[1,:]
            ϕp_I=FC_model.compute_ϕp(var_I,par_mod)
            
            var_II=vars[2,:]
            ϕp_II=FC_model.compute_ϕp(var_II,par_mod)
            if length(output)<1
                output=vcat(var_I,ϕp_I,var_II,ϕp_II)'
            else
                new_point=vcat(var_I,ϕp_I,var_II,ϕp_II)'
                distance=sqrt(minimum(sum((1 .-output./new_point).^2,dims=2)))
                temp=sum((new_point.-output).^2,dims=2)
                if distance>1e-1
                    output=[output; var_I' ϕp_I' var_II' ϕp_II']
                end
            end
        end
    end
    print("output computed correctly\n\n")

    ################################################
    return output
    end
end

