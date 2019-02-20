function BProp2(data::Array{Any,2},W::Array{Any,1},b::Array{Any,1},L::Array{Int64,2})
    N=length(L); #number of layers, so N-2 middle layers
    n=sum(L[2:end])+sum(L[2:end].*L[1:end-1]);#total number of NN parameters
    LC=vcat(0,L[2:end].*L[1:end-1]+L[2:end]);
    LC=cumsum(vec(LC));

    xb=data[1];
    yb=data[2];
    a=Any[];
    z=Any[];
    S=Any[];
    # display(yb)
    #forward iteration

    push!(a,xb);
    push!(z,[]);
    push!(S,[]);

    G(x)=(exp.(2*x) .- 1.0) ./ (exp.(2*x) .+ 1.0);
    GP(x)=( 4.0.*exp.(2.0.*x) )./( (1.0 .+ exp.(2x)).^2 );

    for k=2:N
        push!(z,W[k]*a[k-1]+b[k]);
        push!(a,G.(z[k]));
        push!(S,Diagonal(GP.(vec(z[k]))));
    end
    # display(size(z))
    # display(size(a))
    # display(N-1)
    # push!(z,W[N]*a[N-1]+b[N]);
    # push!(a,z[N]); #final activation function is f(x)=x;
    # push!(S,Matrix{Float64}(I, L[end],L[end])); #same

    #initialize gradients
    ∂Qb=copy(b);
    ∂QW=copy(W);

    #backward iteration
    # Here we are using special cost function F(x)=5x^2 for x<0 and x@ x>=0
    # function DCostCYABV(x)
    #     if x>=0
    #         y=2.0*x
    #     else
    #         y=10.0*x
    #     end
    #     return y
    # end

    # ∂Qb[N]=GP(z[N])*DCostCYABV(a[N] - [yb]);
    ∂Qb[N]=GP(z[N]).*(a[N] - [yb]);

    for k=1:N-2
        ∂Qb[N-k]=S[N-k] * W[N-k+1]' * ∂Qb[N-k+1];
    end

    for k=2:N
        ∂QW[k]=∂Qb[k] * a[k-1]';
    end

    #reshape gradients into single array to match w

    ∂Qw=zeros(n,1);

    for k=2:N
        for i=1:L[k-1]
            ∂Qw[ LC[k-1] + (i-1)*L[k] + 1 : LC[k-1] + i*L[k]  ] = ∂QW[k][:,i];
        end
            ∂Qw[ LC[k-1] + L[k-1]*L[k] + 1 : LC[k-1] + (L[k-1]+1)*L[k] ] = ∂Qb[k];
    end

    # function CostCYABV(x)
    #     if x>=0
    #         y=x^2
    #     else
    #         y=5*x^2
    #     end
    #     return y
    # end
    # Q=CostCYABV(a[N] - [yb]);
    Q=.5*norm(a[N] - [yb]);
    ∂Q=∂Qw;

    return (∂Q,Q)

end
