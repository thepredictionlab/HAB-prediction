using LinearAlgebra
using Random
using Plots



site=5;

datesm=(convert(Array{Float64,2},Matrix(dates)));
datesy=round.(datesm,RoundDown);

cyanm=Matrix(cyano)[:,site];
tox4m=Matrix(tox4)[:,site];
tempm=Matrix(temp);
rainm=Matrix(rain);
tbv14=Matrix(tbv1)[:,site];
tbv24=Matrix(tbv2)[:,site];
tbv34=Matrix(tbv3)[:,site];
tbv44=Matrix(tbv4)[:,site];
satdatesm=Matrix(sattime);
colorm=Matrix(COL);

#interpolate cyanobacteria
cya0=copy(cyanm);
indn = findall(x-> isnan(x), cyanm);
cya0[indn]=zeros(length(indn),1);


#color Data
colori=convert(Array{Float64},colorm);
satdatesi=convert(Array{Float64},satdatesm);
col=zeros(length(datesm),7);
dcol=round.(satdatesi,digits=2);
testsat=zeros(length(satdatesi),1);
for k=1:length(satdatesi)
    ind=findall(x-> x==dcol[k], round.(datesm,digits=2)[:] )
    for j=1:length(ind)
        col[ind[j],:]=colori[k,:];
    end
    if length(ind)>0
        testsat[k]=1
    end
end

#interpolate color
indcolor=findall(x-> x>0.0, col[:,1]);
colint=copy(col);
for k=1:7
    colint[1:(indcolor[1]-1),k]=(collect(1:(indcolor[1]-1))).*(col[indcolor[1],k])/indcolor[1];
end

for k=1:length(indcolor)-1
    for j=1:7
        colint[indcolor[k]+1: indcolor[k+1]-1,j]=(collect(1:indcolor[k+1]-1-indcolor[k] )).*(col[indcolor[k+1],j])/(indcolor[k+1]-indcolor[k]);
    end
end


#interpolate temp and rain
raini=(convert(Array{Float64,2},rainm));
indnr=findall(x-> isnan(x), raini[:]);
for k=1:length(indnr)
    raini[indnr[k]]=0;
end

tempi=(convert(Array{Float64,2},tempm));
indnt=findall(x-> isnan(x), tempi[:]);
for k=1:length(indnt)
    tempi[indnt[k]]=0;
end
indnt=findall(x-> x==1.0, tempi[:]);
for k=1:length(indnt)
    tempi[indnt[k]]=0;
end

for k=1:length(tempi)
    if tempi[k]==0.0
        if k<length(tempi)
            if tempi[k+1]>0.0 && tempi[k-1]>0.0
                tempi[k]=.5*(tempi[k-1]+tempi[k+1]);
            end
            if tempi[k-1]>0.0
                tempi[k]=tempi[k-1]
            end
        end
        tempi[k]=tempi[k-1];
    end

end

# moving average
tempiw=zeros(size(tempi));
rainiw=zeros(size(raini));
colw=zeros(size(colint));
cyaw=zeros(size(cya0));

for k=8:length(tempi)
    tempiw[k]=maximum(tempi[k-7:k]);
    rainiw[k]=maximum(raini[k-7:k]);
    colw[k,:]=maximum(colint[k-7:k,:],dims=1);
    cyaw[k]=maximum(cya0[k-7:k]);
end

for k=1:7
    tempiw[k]=maximum(tempi[1:k]);
    rainiw[k]=maximum(raini[1:k]);
    colw[k,:]=maximum(col[1:k,:],dims=1);
    cyaw[k]=maximum(cya0[1:k]);
end

tempim=zeros(size(tempi));
rainim=zeros(size(raini));

for k=31:length(tempi)
    tempim[k]=maximum(tempi[k-30:k]);
    rainim[k]=maximum(raini[k-30:k]);
end

for k=1:30
    tempim[k]=maximum(tempi[1:k]);
    rainim[k]=maximum(raini[1:k]);
end

# day of year
doy=zeros(size(datesm));
for k=1:length(datesm)
    doy[k]=round(modf(datesm[k])[1]*365);
end

# X0=[tempiw rainiw doy];
X0=[tempiw rainiw colw[:,4] cyaw];
# X0=[tempiw rainiw];
# X0=doy;

# ## microcystin
#
# itox=findall(x-> !isnan(x), tox4m)
# datestox=datesm[itox];
#
# Ytox=tox4m[itox];
# Xtox=X0[itox,:];
#
# #shift
# Xtox=Xtox[7:end,:];
# Ytox=Ytox[1:end-7];
#
# mxtox=[sum(Xtox[:,1].^2)^.5 sum(Xtox[:,2].^2)^.5];
# mytox=[sum(Ytox.^2)];
#
# Xtox=Xtox./mxtox;
# Ytox=Ytox./mytox;

## cyanobacteria

#shift
X1=X0[7:end,:];
Y1=cyanm[1:end-7];
dates0=datesm[7:end-7];

icya=findall(x-> !isnan(x), cyanm)
datescya=datesm[icya];

Ycya=Y1[icya];
Xcya=X1[icya,:]

#normalize Y;
Ycya=Ycya/1e6;
# Ycya=log.(1.0 .+ Ycya);

# mXcya=[sum(Xcya[:,1].^2)^.5 sum(Xcya[:,2].^2)^.5];
mXcya=sum(Xcya.^2, dims=1).^.5;
mYcya=[sum(Ycya.^2)];

Xcya=Xcya./mXcya;
# Ycya=Ycya./mYcya;

# Make Ycya categorical
indlow=findall(x-> x<.25, Ycya);
Ycya=ones(length(Ycya),1);
Ycya[indlow]=-1.0.*ones(length(indlow),1);

ind=findall(x-> x<2018, datescya);
indt=collect(1:10:length(ind))
indcyatrain=setdiff(ind,indt)
indcyaval=indt;
Dtcya=[];
for k=1:length(indcyatrain)
    push!(Dtcya,[[Xcya[indcyatrain[k],:]] [Ycya[indcyatrain[k]]]]);
end
# Dtcya=[[Xcya[indcyatrain,:]] [Ycya[indcyatrain]]];


indcya2018=findall(x-> x>=2018, datescya);
D2018cya=[[Xcya[indcya2018]] [Ycya[indcya2018]]];
