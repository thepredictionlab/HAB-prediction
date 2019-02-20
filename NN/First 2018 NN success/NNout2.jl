function NNout2(x0::Array{Any,1},μ::Array{Float64,2},L::Array{Int64,2})
    x=x0[1];
    z=zeros(size(x,1),L[end]);
    N=length(L); #number of layers, so N-2 middle layers
    n=sum(L[2:end])+sum(L[2:end].*L[1:end-1]); #total number of NN parameters
    w=μ;
    W=Any[];
    b=Any[];
    push!(W,[]);
    push!(b,[]); #so they are 1xN

    LC=vcat(0,L[2:end].*L[1:end-1]+L[2:end]);
    LC=cumsum(vec(LC));

    for k=2:N
        push!(W, zeros(L[k],L[k-1]) );
        push!(b, zeros(L[k],1) );
        for i=1:L[k-1]
            W[k][:,i]=w[LC[k-1] + (i-1)*L[k] + 1 : LC[k-1] + i*L[k] ];
        end
            b[k]=w[LC[k-1] + L[k-1]*L[k] + 1 : LC[k-1] + (L[k-1]+1)*L[k] ];
    end

    W=W[2:end];
    b=b[2:end];

    for kd=1:size(x,1)

            M=length(W);
            yout=x[kd,:];
            G(x)=tanh(x); #activation function
            for k=1:M
                yout=G.(W[k]*yout+b[k]);
            end

        z[kd,:]=yout[:];
    end

    return(z)
end
