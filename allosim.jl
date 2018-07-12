using(DataFrames)
using(CSV)
using(RCall)
using(Distributions)
include("$(homedir())/Dropbox/Postdoc/2018_eigenvec/src/laplacian.jl")
include("$(homedir())/Dropbox/Postdoc/2018_eigenvec/src/eigencluster.jl")
# include("$(homedir())/2018_eigenvec/src/laplacian.jl")
# include("$(homedir())/2018_eigenvec/src/eigencluster.jl")


nsp = 100;
sp = collect(1:nsp);
nmeas = 3;

b0error = 0.000001;
b0dist = Normal(0,b0error);
etaerror = 0.01;
etadist = Normal(0,etaerror);

b0 = rand(collect(0.01:0.001:10.0),nmeas);
eta = repeat([0.75],outer=nmeas);
M = rand(collect(1:0.1:1000),nsp);

data = Array{Float64}(nsp,nmeas);
for i=1:nsp
    for j=1:nmeas
        data[i,j] = abs(b0[j]*(1+rand(b0dist))*M[i]^(eta[j]*(1+rand(etadist))));
    end
end


PC = SharedArray{Float64}(nsp,nsp);
# measmeans = mean(pcdatatr[!isnothing(pcdatatr)],1);
#Build similarity matrix
@time @sync @parallel for i = 0:(nsp^2 - 1)
    a = mod(i,nsp) + 1;
    b = Int64(floor(i/nsp)) + 1;
    if a == b
        PC[a,b] = 0.0;
        continue
    end
    ct = 0;
    ctones = 0;
    for j = 1:nmeas
        ct += log(minimum([data[a,j],data[b,j]])/maximum([data[a,j],data[b,j]]));
        ctones += 1;
    end
    ctscaled = exp(ct/ctones);
    PC[a,b] = Float64(ctscaled); #/Float64(ctones);s
end

S = laplacian(PC,10);
ev = eigs(S; nev=10,which=:SR);
eval = ev[1];
evecs = ev[2];


namespace = string("$(homedir())/2018_eigenvec/figures/allosim.pdf");
R"""
pdf($namespace, height = 12, width = 15)
par(mfrow=c(2,2))

library(RColorBrewer)
pal = colorRampPalette(brewer.pal(11,"Spectral"))($nsp)
plot($(evecs[:,2]),$(evecs[:,3]),col=pal[$(sortperm(sortperm(M)))],pch=16)

library(RColorBrewer)
pal = colorRampPalette(brewer.pal(11,"Spectral"))($nsp)
plot($(evecs[:,3]),$(evecs[:,4]),col=pal[$(sortperm(sortperm(M)))],pch=16)
# text($(evecs[:,3]),$(evecs[:,4]),$sp,cex=0.5)

# library(RColorBrewer)
# pal = colorRampPalette(brewer.pal(11,"Spectral"))($nsp)
# plot($(evecs[:,3]),$(evecs[:,5]),col=pal[$(sortperm(sortperm(M)))],pch=16)
# # text($(evecs[:,3]),$(evecs[:,5]),$sp,cex=0.5)

library(scatterplot3d) 
library(RColorBrewer)
pal = colorRampPalette(brewer.pal(11,"Spectral"))($nsp);
s3d = scatterplot3d(x=cbind($(evecs[:,2]),$(evecs[:,3]),$(evecs[:,4])),color=pal[$(sortperm(sortperm(M)))],pch=16,xlab='ev2',ylab='ev3',zlab='ev4',scale.y=1,angle=80,type='h');
# text(s3d$xyz.convert(cbind($(evecs[:,2]),$(evecs[:,3]),$(evecs[:,4]))),labels=$sp,cex=0.5);


pal = colorRampPalette(brewer.pal(11,"Spectral"))($nsp);
s3d = scatterplot3d(x=cbind($(data[:,1]),$(data[:,2]),$(data[:,3])),color=pal[$(sortperm(sortperm(M)))],pch=16,xlab='ev2',ylab='ev3',zlab='ev4',scale.y=1,angle=80,type='h');

dev.off()
"""


#If nmeas is low-D, we can visualize the data
# R"""
# par(mfrow=c(1,1))
# library(scatterplot3d) 
# library(RColorBrewer)
# pal = colorRampPalette(brewer.pal(11,"Spectral"))($nsp);
# s3d = scatterplot3d(x=cbind($(data[:,1]),$(data[:,2]),$(data[:,3])),color=pal[$(sortperm(sortperm(M)))],pch=16,xlab='ev2',ylab='ev3',zlab='ev4',scale.y=1,angle=80,type='h');
# # dev.off()
# """
