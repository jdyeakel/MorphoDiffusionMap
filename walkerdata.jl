using(DataFrames)
using(CSV)
using(RCall)
include("$(homedir())/Dropbox/Postdoc/2018_eigenvec/src/laplacian.jl")
include("$(homedir())/Dropbox/Postdoc/2018_eigenvec/src/eigencluster.jl")

##############################
# HERBIVORE POSTCRANIAL DATA
##############################

pcdata = CSV.read("$(homedir())/Dropbox/Postdoc/2018_eigenvec/pc_bones.csv",header=true);
sp = pcdata[:common_name];
genus = pcdata[:Genus];
species = pcdata[:Species];
_array = Array{Char}(length(genus)); _array[1:length(genus)] = ' ';
gensp = mapslices(join, [Array{String}(genus) _array Array{String}(species)], 2);
nsp = size(pcdata)[1];
nmeas = size(pcdata)[2];

#take out body mass data info
# remove = [find(x->x==:mass_max_kg,names(pcdata));find(x->x==:mass_min_kg,names(pcdata))]
# keep = deleteat!(collect(1:nmeas),remove);
# pcdata = pcdata[:,keep];
# nmeas = size(pcdata)[2];
bodymass = Array(pcdata[:mass_max_kg]);
pcdatatr = copy(pcdata[:,6:nmeas]);
#Delete measurements that are missing for ALL species
nmeas = size(pcdatatr)[2];
todelete = Array{Int64}(0);
for i=1:nmeas
    if all(ismissing.(pcdatatr[:,i]))
        push!(todelete,i);
    end
end
pcdatatr = pcdatatr[:,setdiff(collect(1:nmeas),todelete)];
meas = Array{String}(names(pcdatatr));
nmeas = size(pcdatatr)[2];

# pcdatatr = Array(pcdatatr) ./ Array(pcdatatr[:femur_length]);


#scale to 'mass'
fl = copy(Array(pcdatatr[:femur_length]));
for i=1:nmeas
    pcdatatr[:,i] = Array(pcdatatr[:,i]) ./ fl;
end


nomeas = Array{Int64}(nsp);
for i=1:nsp
    nomeas[i] = length(find(ismissing,Array(pcdatatr[i,:])));
end
bydataamt = sortperm(nomeas);

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
        if !ismissing(pcdatatr[a,j]) && !ismissing(pcdatatr[b,j])
            ct += log(minimum([pcdatatr[a,j],pcdatatr[b,j]])/maximum([pcdatatr[a,j],pcdatatr[b,j]]));
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

ranked = sortperm(evecs[:,2]);


namespace = string("$(homedir())/Dropbox/Postdoc/2018_eigenvec/figures/walker.pdf");
R"""
pdf($namespace, height = 12, width = 15)
par(mfrow=c(2,2))

library(RColorBrewer)
pal = colorRampPalette(brewer.pal(11,"Spectral"))($nsp)
plot($(evecs[:,2]),$(evecs[:,3]),col=pal[$(sortperm(sortperm(bodymass)))],pch=16)

library(RColorBrewer)
pal = colorRampPalette(brewer.pal(11,"Spectral"))($nsp)
plot($(evecs[:,3]),$(evecs[:,4]),col=pal[$(sortperm(sortperm(bodymass)))],pch=16)
# text($(evecs[:,3]),$(evecs[:,4]),$sp,cex=0.5)

library(RColorBrewer)
pal = colorRampPalette(brewer.pal(11,"Spectral"))($nsp)
plot($(evecs[:,3]),$(evecs[:,5]),col=pal[$(sortperm(sortperm(bodymass)))],pch=16)
# text($(evecs[:,3]),$(evecs[:,5]),$sp,cex=0.5)

library(scatterplot3d) 
library(RColorBrewer)
pal = colorRampPalette(brewer.pal(11,"Spectral"))($nsp);
s3d = scatterplot3d(x=cbind($(evecs[:,2]),$(evecs[:,3]),$(evecs[:,4])),color=pal[$(sortperm(sortperm(bodymass)))],pch=16,xlab='ev2',ylab='ev3',zlab='ev4',scale.y=1,angle=80,type='h');
text(s3d$xyz.convert(cbind($(evecs[:,2]),$(evecs[:,3]),$(evecs[:,4]))),labels=$sp,cex=0.5);
dev.off()
"""



R"""
image($(PC[ranked,ranked]))
"""

R"""
library(RColorBrewer)
pal = colorRampPalette(brewer.pal(11,"Spectral"))($nsp)
plot($(evecs[:,2]),$(evecs[:,3]),col=pal[$(sortperm(sortperm(bodymass)))],pch=16)
text($(evecs[:,2]),$(evecs[:,3]),$sp,cex=0.5)
"""

R"""
library(RColorBrewer)
pal = colorRampPalette(brewer.pal(11,"Spectral"))($nsp)
plot($(evecs[:,2])*-1,$(evecs[:,4]),col=pal[$(sortperm(sortperm(bodymass)))],pch=16)
text($(evecs[:,2])*-1,$(evecs[:,4]),$sp,cex=0.5)
"""


R"""
library(scatterplot3d) 
library(RColorBrewer)
pal = colorRampPalette(brewer.pal(11,"Spectral"))($nsp)
s3d = scatterplot3d(x=cbind($(evecs[:,2])*-1,$(evecs[:,3])*-1,$(evecs[:,4])*-1),color=pal[$(sortperm(sortperm(nomeas)))],pch=16,xlab='ev2',ylab='ev3',zlab='ev4',scale.y=1,angle=80,type='h')
text(s3d$xyz.convert(cbind($(evecs[:,2])*-1,$(evecs[:,3])*-1,$(evecs[:,4])*-1)),labels=$sp,cex=0.5)
"""

#3D plot with plotly
R"""
library(plotly)
t <- list(
  family = "sans serif",
  size = 14,
  color = toRGB("grey50"))
species = $sp;
df = data.frame(species,$(evecs[:,2])*-1,$(evecs[:,3])*-1,$(evecs[:,4])*-1);
colnames(df) = c('sp','ev2','ev3','ev4');
p <- plot_ly(df, x = ~ev2, y = ~ev3, z = ~ev4,
        mode = 'text',
        text = ~species,
        textposition = 'middle right',
        marker = list(color = ~ev2, colorscale = c('#FFE1A1', '#683531'), showscale = TRUE)) %>%
        add_markers() %>%
        add_text(textfont = t, textposition = "top right") %>%
  layout(scene = list(xaxis = list(title = 'ev2'),
                     yaxis = list(title = 'ev3'),
                     zaxis = list(title = 'ev4')),
         annotations = F)
"""



baskspcodes = CSV.read("$(homedir())/Dropbox/Postdoc/2018_eigenvec/baskervilledata/Table_S1.csv",header=true);
baskcode = Array(baskspcodes[:code]);
baskgensp = Array(baskspcodes[:species]);


basktrophic = CSV.read("$(homedir())/Dropbox/Postdoc/2018_eigenvec/baskervilledata/Table_S2.csv",header=true);
preds = basktrophic[:pred_code];
preys = basktrophic[:prey_code];
npreds = length(preds);

trophic = Array{String}(length(preds),2);
predsunique = unique(preds);
preysunique = unique(preys);
for i = 1:length(predsunique)
    id = find(x->x==predsunique[i],preds);
    codetosp = find(x->x==predsunique[i],baskcode);
    trophic[id,1] = repeat(baskgensp[codetosp],outer=length(id));
end
for i = 1:length(preysunique)
    id = find(x->x==preysunique[i],preys);
    codetosp = find(x->x==preysunique[i],baskcode);
    trophic[id,2] = repeat(baskgensp[codetosp],outer=length(id));
end

R"""
library(RColorBrewer)
pal = colorRampPalette(brewer.pal(11,"Spectral"))($nsp)
plot($(evecs[:,2]),$(evecs[:,4]),col=pal[$(sortperm(sortperm(nomeas)))],pch=16)
text($(evecs[:,2]),$(evecs[:,4]),$sp,cex=0.5)
"""
for i=1:length(trophic[:,1])
    if in(trophic[i,1],gensp) && in(trophic[i,2],gensp)
        # println(trophic[i,1], " -> ",trophic[i,2])
        #draw a link
        eigencoord1 = [evecs[:,2][find(x->x==trophic[i,1],gensp)];evecs[:,4][find(x->x==trophic[i,1],gensp)]];
        eigencoord2 = [evecs[:,2][find(x->x==trophic[i,2],gensp)];evecs[:,4][find(x->x==trophic[i,2],gensp)]];
        R"""
        coords = rbind($(eigencoord1),$(eigencoord2))
        lines(x=coords[,1],y=coords[,2],col=alpha('gray',0.5),lwd=0.5)
        """
    end
end

namespace = "$(homedir())/Dropbox/PostDoc/2018_eigenvec/walkertrophic.pdf";
R"""
library(scatterplot3d) 
library(RColorBrewer)
pal = colorRampPalette(brewer.pal(11,"Spectral"))($nsp)
s3d = scatterplot3d(x=cbind($(evecs[:,2])*-1,$(evecs[:,3])*-1,$(evecs[:,4])*-1),color=pal[$(sortperm(sortperm(nomeas)))],pch=16,xlab='ev2',ylab='ev3',zlab='ev4',scale.y=0.5,angle=88,type='h')
text(s3d$xyz.convert(cbind($(evecs[:,2])*-1,$(evecs[:,3])*-1,$(evecs[:,4])*-1)),labels=$sp,cex=0.5)
"""
for i=1:length(trophic[:,1])
    if in(trophic[i,1],gensp) && in(trophic[i,2],gensp)
        # println(trophic[i,1], " -> ",trophic[i,2])
        #draw a link
        eigencoord1 = [evecs[:,2][find(x->x==trophic[i,1],gensp)]*-1;evecs[:,3][find(x->x==trophic[i,1],gensp)]*-1;evecs[:,4][find(x->x==trophic[i,1],gensp)]*-1];
        eigencoord2 = [evecs[:,2][find(x->x==trophic[i,2],gensp)]*-1;evecs[:,3][find(x->x==trophic[i,2],gensp)]*-1;evecs[:,4][find(x->x==trophic[i,2],gensp)]*-1];
        R"""
        coords = rbind($(eigencoord1),$(eigencoord2))
        lines(s3d$xyz.convert(cbind(coords[,1],coords[,2],coords[,3])),col=alpha('gray',0.5),lwd=0.5)
        """
    end
end


R"""
# scatter3D($(evecs[:,2]),$(evecs[:,3]),$(evecs[:,4]),pch=16,colkey=F,add=T)
# text3D($(evecs[:,2]),$(evecs[:,3]),$(evecs[:,4]),labels=$sp,cex=0.5,add=T,colkey=F)
dev.off()
"""

zdata = CSV.read("$(homedir())/Dropbox/Postdoc/2018_eigenvec/Zliobaitedata/Zliobaitedata.csv",header=true);
padata = Array(zdata[:,10:22]);
zsite = CSV.read("$(homedir())/Dropbox/Postdoc/2018_eigenvec/Zliobaitedata/Zliobaitesite.csv",header=true);

site = Array(zsite[:site]);
abbrev = Array(zsite[:abbrev]);
zsp = Array(zdata[:species]);
padata = Array(zdata[:,10:22]);
nsite = length(site);

spmissing = Array{Int64}(nsite);
for ii = 1:nsite
    zspsite = zsp[find(isodd,padata[:,ii])];
    localsp = intersect(zspsite,gensp);
    spmissing[ii] = length(setdiff(zspsite,gensp));
end



namespace = "$(homedir())/Dropbox/PostDoc/2018_eigenvec/walkerbysite2.pdf";
R"""
library(plot3D)
library(scales)
pdf($namespace,width=15,height=10)
par(mfrow=c(3,5),
    oma = c(5,4,0,0) + 0.1, 
    mar = c(0,0,1,1) + 0.1)
"""

for ii = sortperm(spmissing)
    zspsite = zsp[find(isodd,padata[:,ii])];
    localsp = intersect(zspsite,gensp);
    # spmissing[ii] = length(setdiff(zspsite,gensp));
    keep = indexin(localsp,gensp)
    
    R"""
    scatter3D($(evecs[keep,2]),$(evecs[keep,3]),$(evecs[keep,4]),pch=16,colkey=F,theta=0,phi=0,ticktype="detailed",xlab='EV2',ylab='EV3',zlab='EV4',main=$(site[ii]))
    text3D($(evecs[keep,2]),$(evecs[keep,3]),$(evecs[keep,4]),labels=$(sp[keep]),cex=0.5,add=T,colkey=F)
    """
    
end
R"dev.off()"








#PCA of the same dataset
using MultivariateStats
dataset = pcdatatr;
nmeas = size(dataset)[2];
#turn every missing datapoint into the mean for that measurement
d = Array(dataset);
for i=1:nmeas
    d[find(ismissing,d[:,i]),i] = mean(d[find(!ismissing,d[:,i]),i]);
end
d[find(ismissing,d)] = NaN;
d = Array{Float64}(d);
fit(PCA,d)
R"""
library(ggfortify)
data = $d;
rownames(data) = $sp;
colnames(data) = $meas;
autoplot(prcomp(data),label=T,label.size=2,loadings=T,loadings.label=T,loadings.label.size  = 2)
"""

