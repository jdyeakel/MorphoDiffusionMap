using(DataFrames)
using(CSV)
using(RCall)
eigencluster = function(sp,evecs,n)
    carray = Array{Array}();
    carray[1] = sortperm(evecs[:,2]);
    println("Ordered List ","=",sp[carray[1]])
    for i=2:n
        eigenvec = evecs[:,i];
        rarray = Array{Array}(0);
        for j = 1:length(carray)
            set21 = find(x->x>0,eigenvec);
            set22 = find(x->x<0,eigenvec);
    
            set21 = intersect(carray[j],set21);
            set22 = intersect(carray[j],set22);
            if length(set21) > 0
                push!(rarray,set21);
            end
            if length(set22) > 0
                push!(rarray,set22)
            end
        end
        carray = copy(rarray);
    end
    for i=1:length(carray)
        println()
        println("Cluster ",i,"=",sp[carray[i]])
    end
    # println(showall(carray))
end


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
plot($(evecs[:,2]),$(evecs[:,4]))
text($(evecs[:,2]),$(evecs[:,4]),$sp,cex=0.5)
"""

R"""
library(plot3D)
scatter3D($(evecs[:,2]),$(evecs[:,3]),$(evecs[:,4]))
text3D($(evecs[:,2]),$(evecs[:,3]),$(evecs[:,4]),labels=$sp,cex=0.5,add=T)
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

namespace = "$(homedir())/Dropbox/PostDoc/2018_eigenvec/walkertrophic.pdf";
R"""
library(plot3D)
library(scales)
pdf($namespace,width=10,height=8)
scatter3D($(evecs[:,2]),$(evecs[:,3]),$(evecs[:,4]),pch=16,colkey=F,theta=0,phi=0,ticktype="detailed",xlab='EV2',ylab='EV3',zlab='EV4')
text3D($(evecs[:,2]),$(evecs[:,3]),$(evecs[:,4]),labels=$sp,cex=0.5,add=T,colkey=F)
"""
for i=1:length(trophic[:,1])
    if in(trophic[i,1],gensp) && in(trophic[i,2],gensp)
        # println(trophic[i,1], " -> ",trophic[i,2])
        #draw a link
        eigencoord1 = [evecs[:,2][find(x->x==trophic[i,1],gensp)];evecs[:,3][find(x->x==trophic[i,1],gensp)];evecs[:,4][find(x->x==trophic[i,1],gensp)]];
        eigencoord2 = [evecs[:,2][find(x->x==trophic[i,2],gensp)];evecs[:,3][find(x->x==trophic[i,2],gensp)];evecs[:,4][find(x->x==trophic[i,2],gensp)]];
        R"""
        coords = rbind($(eigencoord1),$(eigencoord2))
        lines3D(x=coords[,1],y=coords[,2],z=coords[,3],add=T,colkey=F,col=alpha('gray',0.5),lwd=1,lwd=0.5)
        """
    end
end
R"""
# scatter3D($(evecs[:,2]),$(evecs[:,3]),$(evecs[:,4]),pch=16,colkey=F,add=T)
# text3D($(evecs[:,2]),$(evecs[:,3]),$(evecs[:,4]),labels=$sp,cex=0.5,add=T,colkey=F)
dev.off()
"""



#Build trophic adjacency matrix
#Species in both datasets
intpred = intersect(trophic[:,1],gensp);
intprey = intersect(trophic[:,1],gensp);
ncombsp = length(combsp)
    
        


for i =