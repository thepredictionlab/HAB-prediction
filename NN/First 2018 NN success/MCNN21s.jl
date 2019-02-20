μn=[];
score=[];
MC=100;
for NN=1:5

    # sample network with N node middle layers
    L=[2 NN 1]
    μt=[];
    scoret=zeros(MC,2);
    for k=1:MC
        (μ,iter,f)=NNfit2(Dtcya,L,3000,.15);
        push!(μt,μ);

        # score on test set
        DX=[];
        push!(DX,Xcya[indcyatrain,:]);
        zt=NNout(DX,μ[end],L);
        scoret[k,1]=CostCYABV(zt-Ycya[indcyatrain]);

        # score on validation set
        DX=[];
        push!(DX,Xcya[indcyaval,:]);
        zv=NNout(DX,μ[end],L);
        scoret[k,2]=CostCYABV(zv-Ycya[indcyaval]);

    end
    push!(score,scoret);


end

DXV=[];
push!(DXV,Xcya[indcyaval,:]);
DX=[];
push!(DX,Xcya[indcyatrain,:]);
L=[4 4 1];
n=sum(L[2:end])+sum(L[2:end].*L[1:end-1])

(μ,iter,f)=NNfit2(Dtcya,L,20000,.35);
(f[end]-f[end-1])/f[end-1]
plot(f)

z=NNout2(DX,μ[end],L);
plot(Ycya[indcyatrain,:])
plot!(sign.(z))

z=NNout2(DXV,μ[end],L);
plot(Ycya[indcyaval,:])
plot!(sign.(z))
