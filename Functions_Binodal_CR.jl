module CR_model
export psi,gradient_psi,transform_variable,inv_transform_variable

using LinearAlgebra, SparseArrays,Parameters, Setfield


#### redifination of variables to help satisfying conditions ######
function transform_variable(X,P)
    var_I=transform_variable_I(X,P)
    var_II=transform_variable_II(X[4:end],10^X[3]/(1+10^X[3]))
    return [var_I';var_II']
end
function transform_variable_I(X,a)
    v1=10^a
    temp=(1-v1)*10^X[1]/(10+10^(X[1]))
    v2=temp*(10/(10+10^(X[2])))
    v3=temp*10^(X[2])/(10+10^(X[2]))
    return [v1,v2,v3]
end

function transform_variable_II(X,a)
    v1=a*10^(X[1])/(10+10^(X[1]))
    v2=(10/(10+10^(X[1])))*a
    v3=(1-a)*10^(X[2])/(10+10^(X[2]))
    return [v1,v2,v3]
end
function inv_transform_variable(var_I,var_II)
    inv_var_I=inv_transform_variable_I(var_I)
    inv_var_II=inv_transform_variable_II(var_II)

    return [inv_var_I';inv_var_II']
end
function inv_transform_variable_I(var)
    ϕP,ϕCL,ϕp=var
    a=log10(ϕP)
    temp=ϕCL+ϕp   
    b=1+log10(temp/(1-ϕP-temp))
    c=1+log10(ϕp/ϕCL)

    return [a,b,c]
end
function inv_transform_variable_II(var)
    ϕP,ϕCL,ϕp=var
    tot=ϕP+ϕCL
    a=log10(tot/(1-tot))
    temp=ϕP
    b=1+log10(ϕP/ϕCL)
    temp=ϕp/(1-tot)
    c=1+log10(temp/(1-temp))
    
    return [a,b,c]
end

function charge_moments(var,par)
    @unpack η,α,Nz,χ,λ,Nmono= par
    ϕP,ϕCL,ϕp=var
    κ= sqrt(2*λ*ϕCL)
    
    D=1-ϕP-ϕCL
    ϕps=ϕp/(D-ϕp)    

    A= 0
    Q = 0
    S=0
    
    for z in 0:Nz    
        temp=ϕps.^z.*binomial(Nz,z).*exp(-(α+χ .*ϕP).*z-η.*z.^2 ./(2.0.*Nz)) 
        A+=temp
        Q+=temp.*z
        S+=temp.*z^2
    end  
    Q*=(1/A)
    S*=(1/A)
    S=sqrt(S-Q^2)
    return [Q,S]
end

function psi(vars,par)
    @unpack η,α,Nz,χ,λ,Nmono= par
    ϕP,ϕCL,ϕp=vars

    
    κ= sqrt(2*λ*ϕCL)
   
    ϕs=1-sum(vars)
    ϕps=ϕp/ϕs

    A= 0
    Q=0
    for z in 0:Nz  
        temp=ϕps^z*binomial(Nz,z)*exp(-(α+χ*ϕP)*z-η*z^2/(2.0*Nz))
        Q+=temp.*z
        A+=temp
    end
    Q*=(1/A)
    D=1-ϕP-ϕCL

    term1= ϕP*(log(ϕP)-log(A))/Nmono
    term2= (1-ϕP-ϕCL)*(log(D)-log(1+ϕps))+(1-ϕP-2*ϕCL)*χ*ϕP
    term3= 1/4/pi*(κ/2*(κ+2-κ^2)/(κ+1)-log(1+κ))
    term4= ϕCL*log(ϕCL*ϕps)
    return term1+term2+term3+term4
end


function electroneutral_condition(vars,par)
    @unpack η,α,Nz,χ,λ,Nmono= par
    ϕP,ϕCL,ϕp=vars
    κ = sqrt(2*λ*ϕCL)
    ϕs=1-ϕP-ϕCL-ϕp
    A= 0
    Q=0
    for z in 0:Nz    
        temp=(ϕp/ϕs).^z.*binomial(Nz,z).*exp(-(α+χ .*ϕP).*z-η.*z.^2 ./(2.0.*Nz)) 
        Q+=temp.*z
        A+=temp
    end  

    Q*=(1/A)
    electroneutrality=ϕp-ϕCL.+1/Nmono.*Q*ϕP        

    return electroneutrality
end
function dpsidϕP(vars,par)
    @unpack η,α,Nz,χ,λ,Nmono= par
    ϕP,ϕCL,ϕp=vars

    κ = sqrt(2*λ*ϕCL)
    ϕs=1-sum(vars)

    A= 0
    Q=0
    VARIANCE=0

    for z in 0:Nz    
        temp=(ϕp/ϕs).^z.*binomial(Nz,z).*exp(-(α+χ .*ϕP).*z-η.*z.^2 ./(2.0.*Nz)) 
        Q+=temp.*z
        VARIANCE+=temp.*z^2
        A+=temp
    end  
    Q*=(1/A)
    VARIANCE*=(1/A)    
    dQdϕps = (VARIANCE-Q^2)

    ϕps=ϕp/ϕs
    D=1-ϕP-ϕCL
    ########### electroneutrality: E(φP,φCL,φps(φP,φCL))=0 => ∂φps/∂φP=-(∂E/∂φP)/(∂E/∂φps) ########
    dEdϕps=D/(1+ϕps)^2+1/Nmono*dQdϕps/ϕps*ϕP
    dEdϕP=-ϕps/(1+ϕps)+(Q-ϕP*χ*dQdϕps)/Nmono
    dϕpsdϕP=-dEdϕP/dEdϕps
   
    
    term1= (log(ϕP)-log(A)+1-Q*(dϕpsdϕP/ϕps-χ)*ϕP)/Nmono
    term2= (1-2*ϕP-2*ϕCL)*χ+(-log(D)+log(1+ϕps))-(1-ϕP-ϕCL)*(1/D+1/(1+ϕps)*dϕpsdϕP)
    term3= 0
    term4= ϕCL/ϕps*dϕpsdϕP
    return term1+term2+term3+term4
end


function dpsidϕCL(vars,par)
    @unpack η,α,Nz,χ,λ,Nmono= par
    ϕP,ϕCL,ϕp=vars
    ϕs=1-sum(vars)

    κ = sqrt(2*λ*ϕCL)
    A= 0
    Q=0
    VARIANCE=0

    for z in 0:Nz    
        temp=(ϕp/ϕs).^z.*binomial(Nz,z).*exp(-(α+χ .*ϕP).*z-η.*z.^2 ./(2.0.*Nz))  
        Q+=temp.*z
        VARIANCE+=temp.*z^2
        A+=temp
    end  
    Q*=(1/A)
    VARIANCE*=(1/A)
    
    ϕps=ϕp/ϕs
    dQdϕps = (VARIANCE-Q^2)/ϕps
    D=1-ϕP-ϕCL
    ########### electroneutrality: E(φP,φCL,φps(φP,φCL))=0 => ∂φps/∂φCL=-(∂E/∂φCL)/(∂E/∂φps) ########

    dEdϕps=D/(1+ϕps)^2+dQdϕps*ϕP/Nmono
    dEdϕCL=-(ϕps/(1+ϕps)+1)

    dϕpsdϕCL=-dEdϕCL/dEdϕps
    dκdϕCL = λ/κ
    
    term1= -Q*dϕpsdϕCL/ϕps*ϕP/Nmono
    term2= (-log(D)+log(1+ϕps))-(1-ϕP-ϕCL)*(1/D+1/(1+ϕps)*dϕpsdϕCL)-2*χ*ϕP
    term3= -1/4/pi*(κ^2/(1+κ))*dκdϕCL
    term4= (log(ϕCL*ϕps)+1+ϕCL*dϕpsdϕCL/ϕps)
    return term1+term2+term3+term4
end
function gradient_psi(vars,par)
    v1=dpsidϕP(vars,par)
    v2=dpsidϕCL(vars,par)
    return [v1,v2]
end

end