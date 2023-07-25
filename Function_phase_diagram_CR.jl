module Module_phase_diagram_CR
export compute_binodals
using BifurcationKit, LinearAlgebra, Plots, Parameters, Setfield,DelimitedFiles,Revise
includet("Functions_Binodal_CR.jl")
import ..CR_model
################## coexisting conditions CR model ###########################
function Tangent_construction(X,p)
    @unpack P= p

    f=similar(X) # output vector
    
    vars=CR_model.transform_variable(X,P) # transform into volume fractions

    # phase I: dilute phase variables
    var_I=vars[1,:]
    ϕP_I=var_I[1]
    ϕCL_I=var_I[2]       

    f[1]=CR_model.electroneutral_condition(var_I,p)    

    chemical_potentials_I=CR_model.gradient_psi(var_I,p)
    ψ_I=CR_model.psi(var_I,p) # effective free energy phase I

    # phase II: condensed phase variables
    var_II=vars[2,:]
    ϕP_II=var_II[1]
    ϕCL_II=var_II[2]
    f[2]=CR_model.electroneutral_condition(var_II,p)
    chemical_potentials_II=CR_model.gradient_psi(var_II,p)
    ψ_II=CR_model.psi(var_II,p) # effective free energy phase II

    μP=chemical_potentials_I[1]
    μCL=chemical_potentials_I[2]
    
    f[3]=μP-chemical_potentials_II[1] # continuity effective chemical potential polymer
    f[4]=μCL-chemical_potentials_II[2] # continuity effective chemical potential counter-ions

    f[5]=(ψ_II-ψ_I-μP*(ϕP_II-ϕP_I)-μCL*(ϕCL_II-ϕCL_I)) # continuity effective osmotic pressure

    return f
end


function compute_binodal(fixpar)
    par_mod=(η=fixpar.η,α=fixpar.α,Nz=fixpar.Nz,χ=fixpar.χ,P=log10(0.05),λ=fixpar.λ,Nmono=fixpar.Nmono);
    
    ############# load the initial guess #############
    name="initial_condition_alpha_-6.5chi_0.95"

    M=readdlm(name*".txt", ',', Float64, '\n');
    indeces=sortperm((M[1:end,1].-par_mod.η).^2);
    idx=[indeces[1]]
    n=1
    for el in indeces[2:20]
        X=minimum(abs.(el.-indeces[1:n]))
        if X>1 
            append!(idx,el)
        end
        n+=1
    end
   
    x0_vec=sortperm(M[idx,7]);
    ############ 
    branches=Dict()
    sol0=zeros(5,) 
    n=0
    ############ use continuation to solve equilibrium problem ################
    print("computing binodal\n\n")
    for x0 in [idx[x0_vec[1]],idx[x0_vec[end]]]
        vars=CR_model.inv_transform_variable(M[x0,2:4],M[x0,7:9])
        var_I,var_II=vars[1,:],vars[2,:]

        par_mod=@set par_mod.P=var_I[1]
        sol0[1:2].=var_I[2:3]
        sol0[3:5].=var_II
        opt_newton = NewtonPar(tol = 1e-10, verbose =false,maxIter=50)
        prob = BifurcationProblem(Tangent_construction,sol0,par_mod,(@lens _.P),recordFromSolution= (x, p) -> x[1];);
        out= newton(prob, opt_newton);
        prob = BifurcationKit.reMake(prob, u0 =out.u);
        opt_newton = NewtonPar(tol = 1e-8, verbose =false, maxIter = 3)
        opts_br = ContinuationPar(dsmin=0.00001,dsmax=0.01,ds=0.005, pMin=-10.0,pMax=0.,detectBifurcation=0,
        newtonOptions = opt_newton,saveSolEveryStep=1,maxSteps=2000)
        
        branches[n]=BifurcationKit.continuation(prob,PALC(tangent=Bordered()),opts_br;plot=false,bothside = true);
        n+=1
    end
    ############ processing output in the right format ################
    print("preparing output\n\n")

    output=[]
    for n in 0:length(branches)-1

        for el in branches[n].sol

            vars=CR_model.transform_variable(el.x,el.p) 
            var_I=vars[1,:]
            charge_I=CR_model.charge_moments(var_I,par_mod)
            
            var_II=vars[2,:]
            charge_II=CR_model.charge_moments(var_II,par_mod)
            if length(output)<1
                output=vcat(var_I,charge_I,var_II,charge_II)'
            else
                new_point=vcat(var_I,charge_I,var_II,charge_II)'
                distance=sqrt(minimum(sum((1 .-output./new_point).^2,dims=2)))
                temp=sum((new_point.-output).^2,dims=2)
                if distance>1e-1
                    output=[output; var_I' charge_I' var_II' charge_II']
                end
            end
        end
    end
    print("output computed correctly\n\n")

    ################################################
    return output
    end
end

