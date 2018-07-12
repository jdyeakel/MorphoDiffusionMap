using(DataFrames)
using(CSV)
using(RCall)
include("$(homedir())/Dropbox/Postdoc/2018_eigenvec/src/laplacian.jl")
include("$(homedir())/Dropbox/Postdoc/2018_eigenvec/src/eigencluster.jl")

# df = CSV.read("$(homedir())/robertson_1929_matr.csv",header=false);
# pols = CSV.read("$(homedir())/robertson_1929_pol.csv",header=false);

mdata = CSV.read("$(homedir())/Dropbox/Postdoc/2018_eigenvec/herbdata/Cerlingtrim.csv",header=true);


#Construct presence absence matrix
sp = unique(mdata[:common]);
nsp = length(sp);

loc = sort(unique(mdata[:ecosystem]));
nloc = length(loc);

PA = zeros(Int64,nloc,nsp);
for i = 1:nsp
    cn = sp[i];
    cnrows = find(x->x==cn,mdata[:common]);
    cnloc = unique(mdata[cnrows,:ecosystem]);
    #turn location into positions of loc
    PAones = findin(loc,cnloc)
    
    PA[PAones,i] = 1;
end

cerlingclass = ["F","M","M","M","AF","M","M","M","F","F","F","M","F","M","M","F","M","M","M","M","F","G","M","M","M","AF","M","G","M","M"];

R"""
# pdf("~/2018_eigenvec/herblocmatrix.pdf",height=5,width=6)
image($PA)
# dev.off()
"""
# 
# PA = PA'
# nloc = size(PA)[1];
# nsp = size(PA)[2];

A = Array{Float64}(nloc,nloc);
#Build similarity matrix
for i = 0:(nloc^2 - 1)
    a = mod(i,nloc) + 1;
    b = Int64(floor(i/nloc)) + 1;
    
    if a == b
        A[a,b] = 0.0;
        continue
    end
    # if a==1
    #     println(b)
    # end
    ct = 0;
    ctones = 0;
    for j = 1:nsp
        ct += (PA[a,j]*PA[b,j]) + (1 - PA[a,j])*(1 - PA[b,j]);
        # ctones += PA[a,j] + PA[b,j];
    end
    A[a,b] = Float64(ct); #/Float64(ctones);
    
end


S = laplacian(A,10);

ev = eigs(S; nev=10,which=:SR);
eval = ev[1];
evecs = ev[2];

set21 = find(x->x>0,evecs[:,2]);
set22 = find(x->x<0,evecs[:,2]);
forest = loc[set21]
nonforest = loc[set22]

ranking = sortperm(evecs[:,2]);
loc[ranking];
rankingalph = sortperm(loc[ranking])
# rankingscaled = rankingalph/maximum(rankingalph);

set31 = find(x->x>0,evecs[:,3]);
set32 = find(x->x<0,evecs[:,3]);


#Forests
loc[intersect(set21,set31)]

afroalpine = loc[intersect(set22,set32)]



mapmat = CSV.read("$(homedir())/Dropbox/Postdoc/2018_eigenvec/herbdata/mapmat.csv",header=true);
map = mapmat[:MAP];
mat = mapmat[:MAT];

R"""
library(RColorBrewer)
rankingalph = $rankingalph;
pal = colorRampPalette(brewer.pal(9,"Spectral"))(length($rankingalph))
plot($mat,$map,pch=16,col=pal[rankingalph],xlab='MAT',ylab='MAP')
text($mat+0.75,$map+33,$(loc),cex=0.8)
"""

R"""
image($(PA[ranking,:]))
"""



R"""
pdf("~/matriximage.pdf",height=5,width=6)
image($S)
dev.off()
"""
R"""
pdf("~/adjmatrix.pdf",height=5,width=6)
image($A)
dev.off()
"""


# HERBIVORE ANALYSIS

PA = PA'
nloc = size(PA)[1];
nsp = size(PA)[2];

A = Array{Float64}(nloc,nloc);
#Build similarity matrix
for i = 0:(nloc^2 - 1)
    a = mod(i,nloc) + 1;
    b = Int64(floor(i/nloc)) + 1;
    
    if a == b
        A[a,b] = 0.0;
        continue
    end
    # if a==1
    #     println(b)
    # end
    ct = 0;
    ctones = 0;
    for j = 1:nsp
        ct += (PA[a,j]*PA[b,j]) + (1 - PA[a,j])*(1 - PA[b,j]);
        # ctones += PA[a,j] + PA[b,j];
    end
    A[a,b] = Float64(ct); #/Float64(ctones);
    
end

S = laplacian(A,10);

ev = eigs(S; nev=10,which=:SR);
eval = ev[1];
evecs = ev[2];

ranking2 = sortperm(evecs[:,2]);
sp[ranking2];
rankingscaled = ranking/maximum(ranking);

R"image(x=seq(1,30),y=seq(1,60),z=$(PA[ranking,ranking2]))"

R"plot($(sum(PA[ranking,ranking2],2)))"


set21 = find(x->x>0,evecs[:,2]);
set22 = find(x->x<0,evecs[:,2]);
open = loc[set21]
nonopen = loc[set22]

set31 = find(x->x>0,evecs[:,3]);
set32 = find(x->x<0,evecs[:,3]);


nonrhino = sp[intersect(set31,set22)]
blackrhino = sp[intersect(set32,set22)]

set41 = find(x->x>0,evecs[:,4]);
set42 = find(x->x<0,evecs[:,4]);

savanna = sp[intersect(set41,intersect(set31,set22))]
mixed = sp[intersect(set42,intersect(set31,set22))]

set51 = find(x->x>0,evecs[:,5]);
set52 = find(x->x<0,evecs[:,5]);

sp[intersect(set52,intersect(set42,intersect(set31,set22)))]

set61 = find(x->x>0,evecs[:,6]);
set62 = find(x->x<0,evecs[:,6]);

sp[intersect(set61,intersect(set52,intersect(set42,intersect(set31,set22))))]

#Forests
loc[intersect(set21,set31)]

afroalpine = loc[intersect(set22,set32)]

set41 = find(x->x>0,evecs[:,4]);
set42 = find(x->x<0,evecs[:,4]);

set51 = find(x->x>0,evecs[:,5]);
set52 = find(x->x<0,evecs[:,5]);

##############################
# HERBIVORE POSTCRANIAL DATA
##############################

pcdata = CSV.read("$(homedir())/Dropbox/Postdoc/2018_eigenvec/pc_bones.csv",header=true);
sp = pcdata[:common_name];
genusname = pcdata[:Genus];
speciesname = pcdata[:Species];



nsp = size(pcdata)[1];
nmeas = size(pcdata)[2];

#take out body mass data info
remove = [find(x->x==:mass_max_kg,names(pcdata));find(x->x==:mass_min_kg,names(pcdata))]
keep = deleteat!(collect(1:nmeas),remove);
pcdata = pcdata[:,keep];
nmeas = size(pcdata)[2];
pcdatatr = pcdata[:,6:nmeas];

#scale to 'mass'
pcdatatr = Array(pcdatatr) ./ pcdatatr[:femur_length];

nmeas = size(pcdatatr)[2];
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
    # if a==1
    #     println(b)
    # end
    ct = 0;
    # ct = Array{Float64}(nmeas).*0 + 1;
    ctones = 0;
    for j = 1:nmeas
        if !ismissing(pcdatatr[a,j]) && !ismissing(pcdatatr[b,j])
            # ct += (sqrt((pcdatatr[a,j]-pcdatatr[b,j])^2)); #/mean(pcdatatr[!ismissing.(pcdatatr[:,j]),j]);
            ct += log(minimum([pcdatatr[a,j],pcdatatr[b,j]])/maximum([pcdatatr[a,j],pcdatatr[b,j]]));
            ctones += 1;
        end
        # ctones += PA[a,j] + PA[b,j];
    end
    ctscaled = exp(ct/ctones);
    PC[a,b] = Float64(ctscaled); #/Float64(ctones);
    
end


S = zeros(Float64,nsp,nsp);
# S = copy(-PC);


for i = 1:nsp

    val = zeros(Float64,10);
    lok = zeros(Int64,10);

    for j = 1:nsp
        if PC[i,j] > val[10]
            val[10] = PC[i,j];
            lok[10] = j;
        end
        for k=1:9
            if val[11-k] > val[10-k]
                v = val[11-k];
                val[11-k] = val[10-k];
                val[10-k] = v;
                l = lok[11-k];
                lok[11-k] = lok[10-k];
                lok[10-k] = l;
            end
        end
    end
    S[i,lok] = -PC[i,lok];
    S[lok,i] = -PC[lok,i];
end
rowsums = sum(S,2);
S[diagind(S)] = -rowsums;

ev = eigs(S; nev=10,which=:SR);
eval = ev[1];
evecs = ev[2];

ranked = sortperm(evecs[:,2]);
sp[sortperm(evecs[:,2])]

R"""
image($(PC[ranked,ranked]))
"""

R"""
plot($(evecs[:,2]),$(evecs[:,3]))
text($(evecs[:,2]),$(evecs[:,3]),$sp,cex=0.5)
"""

R"""
library(plot3D)
scatter3D($(evecs[:,2]),$(evecs[:,3]),$(evecs[:,4]))
text3D($(evecs[:,2]),$(evecs[:,3]),$(evecs[:,4]),labels=$sp,cex=0.5,add=T)
"""
    
set21 = find(x->x>0,evecs[:,2]);
intersect(set21,sortperm(evecs[:,3])
v =evecs[:,3]


set22 = find(x->x<0,evecs[:,2]);
sp[set21]
largemammals = sp[set22]

set31 = find(x->x>0,evecs[:,3]);
set32 = find(x->x<0,evecs[:,3]);
sp[set31]
sp[set32]
sp[intersect(set31,set22)]

#Groups carnivores and herbivores by size-similarities (predict trophic relationships?)
smallmammals = sp[intersect(set32,set21)]
mesomammals = sp[intersect(set31,set21)]
megamammals = sp[intersect(set31,set22)]
elephant = sp[intersect(set32,set22)]

set41 = find(x->x>0,evecs[:,4]);
set42 = find(x->x<0,evecs[:,4]);
smallmammals2 = sp[intersect(set42,intersect(set32,set21))]
mediummammals2 = sp[intersect(set41,intersect(set31,set21))]
mesomammals2 = sp[intersect(set41,intersect(set32,set21))]
megamammals2 = sp[intersect(set41,intersect(set31,set22))]
elephant2 = sp[intersect(set42,intersect(set32,set22))]

set51 = find(x->x>0,evecs[:,5]);
set52 = find(x->x<0,evecs[:,5]);


giraffe = sp[intersect(set51,intersect(set41,intersect(set31,set22)))]

smallmammals3 = sp[intersect(set51,intersect(set42,intersect(set32,set21)))]
sp[intersect(set52,intersect(set41,intersect(set31,set21)))]
sp[intersect(set52,intersect(set42,intersect(set32,set21)))]
sp[intersect(set52,intersect(set41,intersect(set31,set22)))]



#Import food web data from Baskerville


