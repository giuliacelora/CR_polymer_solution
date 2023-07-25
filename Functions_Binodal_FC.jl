module FC_model
export psi,gradient_psi,transform_variable

using LinearAlgebra, SparseArrays,Parameters, Setfield


function transform_variable(X,P,par)
    @unpack Nz,Nmono= par
    ξ=Nz/Nmono
    ϕP_I=10^P    
    ϕCL_I=compute_ϕCL(X[1],ϕP_I,ξ)

    ϕP_II=10^X[2]/(1+10^X[2])
    ϕCL_II=compute_ϕCL(X[3],ϕP_II,ξ)

    return [ϕP_I ϕCL_I; ϕP_II ϕCL_II]
end
function compute_ϕCL(X,P,xi)
    temp=10^(X[1])/(10+10^(X[1]))
    top=(1+P*(xi-1))/2
    v1=P*xi+temp*top
    return v1
end
function compute_ϕp(vars,par)
    @unpack Nz,χ,λ,Nmono= par

    ϕP,ϕCL=vars
   
    return ϕCL-Nz/Nmono*ϕP
end
function psi(vars,par)
    @unpack Nz,χ,λ,Nmono= par
    ϕP,ϕCL=vars
    κ= sqrt(2*λ*ϕCL)
    ϕp=ϕCL-Nz/Nmono*ϕP
    ϕs=1-ϕP-ϕCL-ϕp
    term1= ϕP*log(ϕP)/Nmono+ϕs*log(ϕs)+ϕCL*log(ϕCL)+ϕp*log(ϕp)
    term2= ϕs*χ*ϕP
    term3= 1/4/pi*(κ/2*(κ+2-κ^2)/(κ+1)-log(1+κ))
    return term1+term2+term3
end
function dpsidϕP(vars,par)
    @unpack Nz,χ,λ,Nmono= par
    ϕP,ϕCL=vars

    κ= sqrt(2*λ*ϕCL)
    ϕp=ϕCL-Nz/Nmono*ϕP
    ϕs=1-ϕP-ϕCL-ϕp
    term1= (1+log(ϕP))/Nmono-Nz/Nmono*(1+log(ϕp))+(Nz/Nmono-1)*(1+log(ϕs))
    term2= (1-2*ϕP*(1-Nz/Nmono)-2*ϕCL)*χ
    term3= 0
    return term1+term2+term3
end
function dpsidϕCL(vars,par)
    @unpack Nz,χ,λ,Nmono= par
    ϕP,ϕCL=vars

    κ= sqrt(2*λ*ϕCL)
    dκdCL = λ/κ

    ϕp=ϕCL-Nz/Nmono*ϕP
    ϕs=1-ϕP-ϕCL-ϕp
    term1= log(ϕCL)+log(ϕp)-2*log(ϕs)
    term2= -2*ϕP*χ
    term3= -1/4/pi*(κ^2/(1+κ))*dκdCL
    return term1+term2+term3
end
function gradient_psi(vars,par)
    v1=dpsidϕP(vars,par)
    v2=dpsidϕCL(vars,par)
    return [v1,v2]
end

end