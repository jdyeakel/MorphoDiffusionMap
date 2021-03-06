using(DataFrames)
using(CSV)
using(RCall)
include("$(homedir())/Dropbox/Postdoc/2018_eigenvec/src/laplacian.jl")
include("$(homedir())/Dropbox/Postdoc/2018_eigenvec/src/eigencluster.jl")

pcdata = CSV.read("$(homedir())/Dropbox/Postdoc/2018_eigenvec/fieldbirds.csv",header=true);
fpcdata = CSV.read("$(homedir())/Dropbox/Postdoc/2018_eigenvec/fieldbirds.csv",header=true);
# sp = pcdata[:common_name];
genus = pcdata[:Genus];
species = pcdata[:species];
genspinds = mapslices(join, [Array{String}(genus) Array{String}(species)], 2);
gensp = unique(genspinds);
nsp = length(gensp);
pcdatatr = pcdata[:,7:size(pcdata)[2]];
nmeas = size(pcdatatr)[2];
meas = Array{String}(names(pcdatatr));

data = zeros(Float64,nsp,nmeas);
orderlist = Array{Int64}(nsp);
#Find averages for values across species
for i=1:nsp
    inds = find(x->x==gensp[i],genspinds);
    for j=1:nmeas
        data[i,j] = mean(pcdatatr[inds,j]);
    end
    orderlist[i] = find(x->x==String(pcdata[:Order][inds[1]]),unique(pcdata[:Order]))[1];
end
bodymass = data[:,1];
#remove body size data from analysis
data = data[:,2:size(data)[2]];
nmeas = size(data)[2];

PC = Array{Float64}(nsp,nsp);
# measmeans = mean(pcdatatr[!isnothing(pcdatatr)],1);
#Build similarity matrix
for i = 0:(nsp^2 - 1)
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

R"par(mfrow=c(2,2))"

R"""
library(RColorBrewer)
pal = colorRampPalette(brewer.pal(11,"Spectral"))($nsp)
plot($(evecs[:,2]),$(evecs[:,3]),col=pal[$(sortperm(sortperm(bodymass)))],pch=16)
text($(evecs[:,2]),$(evecs[:,3]),$orderlist,cex=0.5)
"""
# 
# R"""
# library(RColorBrewer)
# pal = colorRampPalette(brewer.pal(11,"Spectral"))(max($orderlist))
# plot($(evecs[:,2]),$(evecs[:,3]),col=pal[$(orderlist)],pch=16)
# """


R"""
library(RColorBrewer)
pal = colorRampPalette(brewer.pal(11,"Spectral"))($nsp)
plot($(evecs[:,3]),$(evecs[:,4]),col=pal[$(sortperm(sortperm(bodymass)))],pch=16)
# text($(evecs[:,3]),$(evecs[:,4]),$sp,cex=0.5)
"""


R"""
library(RColorBrewer)
pal = colorRampPalette(brewer.pal(11,"Spectral"))($nsp)
plot($(evecs[:,3]),$(evecs[:,5]),col=pal[$(sortperm(sortperm(bodymass)))],pch=16)
# text($(evecs[:,3]),$(evecs[:,5]),$sp,cex=0.5)
"""

# 
# R"""
# library(RColorBrewer)
# pal = colorRampPalette(brewer.pal(11,"Spectral"))($nsp)
# plot($(evecs[:,3]),$(evecs[:,5]),col=pal[$(sortperm(sortperm(bodymass)))],pch=16)
# # text($(evecs[:,3]),$(evecs[:,5]),$sp,cex=0.5)
# """
# 


R"""
library(scatterplot3d) 
library(RColorBrewer)
pal = colorRampPalette(brewer.pal(11,"Spectral"))($nsp)
s3d = scatterplot3d(x=cbind($(evecs[:,2]),$(evecs[:,3]),$(evecs[:,4])),color=pal[$(sortperm(sortperm(data[:,1])))],pch=16,xlab='ev2',ylab='ev3',zlab='ev4',scale.y=1,angle=80,type='h',xlim=c(-0.1,0.1))
# text(s3d$xyz.convert(cbind($(evecs[:,2])*-1,$(evecs[:,3])*-1,$(evecs[:,4])*-1)),labels=$sp,cex=0.5)
"""




#PCA
d = Array{Float64}(data);
R"""
library(ggfortify)
data = $d;
meas = $(meas)[-1];
rownames(data) = $gensp;
colnames(data) = meas;
autoplot(prcomp(data),label=T,label.size=2,loadings=T,loadings.label=T,loadings.label.size  = 2)
"""

