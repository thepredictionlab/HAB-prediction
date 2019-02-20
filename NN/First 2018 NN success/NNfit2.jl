function NNfit2(Data::Array{Any,1},L::Array{Int64,2},T::Int64,α::Float64)
    N=length(L); #number of layers, so N-2 middle layers
    n=sum(L[2:end])+sum(L[2:end].*L[1:end-1]); #total number of NN parameters
    LC=vcat(0,L[2:end].*L[1:end-1]+L[2:end]);
    LC=cumsum(vec(LC));
    ND=length(Data[1]);
    # initialize posterior parameters
    μ=randn(n,1);
    μs=Any[];
    dQ=Any[];

    # cost function values;
    f = [];
    # global f;
    iter=1;
    stop=0;
    while iter<T && stop==0;
        function SGD1(Data,μ,L)
            N=length(L); #number of layers, so N-2 middle layers
            n=sum(L[2:end])+sum(L[2:end].*L[1:end-1]); #total number of NN parameters
            LC=vcat(0,L[2:end].*L[1:end-1]+L[2:end]);
            LC=cumsum(vec(LC));
            ND=length(Data);
            sgd=shuffle!(collect(1:ND));
            μin=μ;
            for kd=1:ND

                # display(iter)
                Δμ=zeros(n,1);
                Δρ=zeros(n,1);
                fs=0; # cost function values
                Qtot=0;

                w = copy(μ);# .+ 0*.00001*(log.(1.0.+exp.(ρ))).*ϵ; # mc sample weights

                #first build W and b
                W=Any[];
                b=Any[];
                push!(W,[]);
                push!(b,[]); #so they are 1xN

                for k=2:N
                    push!(W, zeros(L[k],L[k-1]) );
                    push!(b, zeros(L[k],1) );
                    # wc = w[LC[k-1]+1 : LC[k]]
                    for i=1:L[k-1]
                        W[k][:,i]=w[LC[k-1] + (i-1)*L[k] + 1 : LC[k-1] + i*L[k] ];
                    end
                        b[k]=w[LC[k-1] + L[k-1]*L[k] + 1 : LC[k-1] + (L[k-1]+1)*L[k] ];
                end

                #compute stochastic gradient of quadratic cost via backprop

                (∂Q,Q)=BProp2([[Data[kd][1]] [Data[kd][2]]] ,W,b,L);

                #Likelihood cost
                Lcost = Q;

                #Complexity cost
                Ccost =  .5 * sum( (w.^2));
                #Likelihood cost
                ∂Lcost = ∂Q;

                #Gradient
                Δμ= Δμ + ∂Q;#
                Δρ= Δρ;

                #Cost function
                fs = fs + Lcost;

                #update parameters
                μ = μ - α*Δμ;

            end
            return(μ,fs);
        end
        (μ,fs)=SGD1(Data,μ,L);

        push!(μs,μ);
        #store cost function value
        push!(f,fs);
        #evaluate stopping criteria
        # if iter > 2 && abs((f[end]-f[end-1])/f[end-1]) < 1e-5;
        #     stop=1;
        # end

        iter=iter+1;
    end
    return(μs,iter,f);
end
