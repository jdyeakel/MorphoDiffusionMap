using(DataFrames)
using(CSV)
using(RCall)
include("$(homedir())/Dropbox/Postdoc/2018_eigenvec/src/laplacian.jl")
include("$(homedir())/Dropbox/Postdoc/2018_eigenvec/src/eigencluster.jl")

pcdata = CSV.read("$(homedir())/Dropbox/Postdoc/2018_eigenvec/serranomodernbirds.csv",header=true);
fpcdata = CSV.read("$(homedir())/Dropbox/Postdoc/2018_eigenvec/serranopaleobirds.csv",header=true);

alldata = vcat(pcdata,fpcdata);

genspinds = alldata[:SPECIES];
gensp = unique(genspinds);
nsp = length(gensp);
pcdatatr = alldata[:,5:size(alldata)[2]];
nmeas = size(pcdatatr)[2];
meas = Array{String}(names(pcdatatr));
yearsbp = Array(alldata[:BP]);

# data = DataFrame();
data = zeros(Float64,nsp,nmeas);
orderlist = Array{Int64}(nsp);
yearsbp2 = Array{Int64}(nsp);
#Find averages for values across species
for i=1:nsp
    inds = find(x->x==gensp[i],genspinds);
    for j=1:nmeas
        meandata = find(!ismissing,pcdatatr[inds,j]);
        if length(meandata) > 0
            data[i,j] = mean(pcdatatr[inds[find(!ismissing,pcdatatr[inds,j])],j]);
        else
            data[i,j] = NaN;
        end
    end
    orderlist[i] = find(x->x==String(alldata[:TAXA][inds[1]]),unique(alldata[:TAXA]))[1];
    yearsbp2[i] = alldata[:BP][inds[1]];
end

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
        if !isnan(data[a,j]) && !isnan(data[b,j])
            #NOTE: there are some negative numbers in the dataset. Right now, I'm assuming these should be positive numbers (taking abs)
            ct += log(minimum([abs(data[a,j]),abs(data[b,j])])/maximum([abs(data[a,j]),abs(data[b,j])]));
            ctones += 1;
        end
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
palbp = numeric($nsp)h
for (i in 1:$nsp) {
    if ($(yearsbp2)[i] == 0) {
        palbp[i] = 'black'
        } else {
        palbp[i] = 'gray'
        }
    palbp[]
    }
plot($(evecs[:,2]),$(evecs[:,3]),col=pal[$(sortperm(sortperm(yearsbp2)))],pch=16)
text($(evecs[:,2]),$(evecs[:,3]),$yearsbp2,cex=0.5)
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
plot($(evecs[:,3]),$(evecs[:,4]),col=pal[$(sortperm(sortperm(yearsbp2)))],pch=16)
# text($(evecs[:,3]),$(evecs[:,4]),$sp,cex=0.5)
"""


R"""
library(RColorBrewer)
pal = colorRampPalette(brewer.pal(11,"Spectral"))($nsp)
plot($(evecs[:,3]),$(evecs[:,5]),col=pal[$(sortperm(sortperm(yearsbp2)))],pch=16)
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
s3d = scatterplot3d(x=cbind($(evecs[:,2]),$(evecs[:,3]),$(evecs[:,4])),color=pal[$(sortperm(sortperm(yearsbp2)))],pch=16,xlab='ev2',ylab='ev3',zlab='ev4',scale.y=1,angle=80,type='h',xlim=c(-0.1,0.1))
# text(s3d$xyz.convert(cbind($(evecs[:,2])*-1,$(evecs[:,3])*-1,$(evecs[:,4])*-1)),labels=$sp,cex=0.5)
"""
