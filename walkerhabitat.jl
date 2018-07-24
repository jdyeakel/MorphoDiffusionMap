using(DataFrames)
using(CSV)
using(RCall)
include("$(homedir())/Dropbox/Postdoc/2018_eigenvec/src/laplacian.jl")
include("$(homedir())/Dropbox/Postdoc/2018_eigenvec/src/eigencluster.jl")
include("$(homedir())/Dropbox/Postdoc/2018_eigenvec/src/walkeranalysis.jl")
include("$(homedir())/Dropbox/Postdoc/2018_eigenvec/src/eigendistance.jl")


pcdata,
pcdatatr,
wgensp,
eval,
evecs = walkeranalysis();

#subselect just the herbivores
wgensp = Array{String}(wgensp);
worder = pcdata[:Order];
wherbs = find(x->x!="Carnivora" && x!="Primates",worder);


#Co-occurance matrix
zdata = CSV.read("$(homedir())/Dropbox/Postdoc/2018_eigenvec/Zliobaitedata/Zliobaitedata.csv",header=true);
zsitenames = CSV.read("$(homedir())/Dropbox/Postdoc/2018_eigenvec/Zliobaitedata/Zliobaitesite.csv",header=true);
zsite = zsitenames[:abbrev];

zgensp = Array{String}(zdata[:species]);
nzsp = length(zgensp);
nzsite = length(zsite);
zorder = zdata[:order];
zherbs = find(x->x!="Primates",zorder);

#Intersection of species in walker and zil datasets
gensp = Array{String}(intersect(wgensp[wherbs],zgensp[zherbs]));

#Co-occurance matrix of species across sites in zilobate dataset
com = Array{Float64}(nzsp,nzsp);
for i=1:nzsp
    for j=1:nzsp
        occ = Array{Int64}(2,13);
        occ[1,:] = Array{Int64}(zdata[i,11:23]);
        occ[2,:] = Array{Int64}(zdata[j,11:23]);
        com[i,j] =sum(prod(occ,1))/sum(occ[1,:]);
    end
end

#QUESTION: for each site, what is the eigenvector interval distance for species in that site?

# veci = 2;

mo = Array{Float64}(nzsite);
rmo = Array{Float64}(nzsite);
eiginterval = Array{Float64}(nzsite);
reiginterval = Array{Float64}(nzsite);
for i=1:nzsite
    #find the species that are in the site, but also in the walker database
    site = find(x->x==zsite[i],Array{String}(names(zdata)));
    occurance = Array(zdata[:,site]);
    spsite = intersect(zgensp[find(isodd,occurance)],gensp);
    loc = Array{Int64}(length(spsite));
    wloc = Array{Int64}(length(spsite));
    for j=1:length(spsite)
        loc[j] = find(x->x==spsite[j],zgensp)[1];
        wloc[j] = find(x->x==spsite[j],wgensp)[1];
    end
    occvalue = com[loc,loc];
    mo[i] = mean(occvalue[setdiff(collect(1:length(loc)^2),diagind(occvalue))]);
    
    eigmatrix = Array{Float64}(length(spsite),length(spsite));
    for j=1:length(spsite)
        for k=1:length(spsite)
            # eigmatrix[j,k] = sqrt((evecs[wloc,veci][j]-evecs[wloc,veci][k])^2);
            eigmatrix[j,k] = eigendistance(j,k,evecs,wloc);
        end
    end
    # eiginterval[i] = std(evecs[wloc,veci]);
    eiginterval[i] = mean(eigmatrix);
    
    its = 1000;
    rrmo = SharedArray{Float64}(its);
    rreiginterval = SharedArray{Float64}(its);
    @sync @parallel for j=1:its
        
        rloc = rand(zherbs,length(spsite));
        roccvalue = com[rloc,rloc];
        rrmo[j] = mean(roccvalue[setdiff(collect(1:length(rloc)^2),diagind(roccvalue))]);
        
        rloc2 = rand(wherbs,length(spsite));
        # rreiginterval[j] = std(evecs[rloc2,veci]);
        
        rreigmatrix = Array{Float64}(length(spsite),length(spsite));
        for l=1:length(spsite)
            for k=1:length(spsite)
                # rreigmatrix[l,k] = sqrt((evecs[rloc2,veci][l]-evecs[rloc2,veci][k])^2);
                rreigmatrix[l,k] = eigendistance(l,k,evecs,rloc2);
            end
        end
        # eiginterval[i] = std(evecs[wloc,veci]);
        rreiginterval[j] = mean(rreigmatrix);
        
        
    end
    rmo[i] = mean(rrmo)
    reiginterval[i] = mean(rreiginterval);
end
R"""
par(mfrow=c(1,1))
boxplot($([eiginterval,reiginterval]),ylim=c(0.2,0.4))
"""
        
zsite[sortperm(eiginterval)]





#Question: what is the average similarity of the species in the eigenclusters relative to the average similarity?
MO = Array{Array}(10);
for ii=1:10
eclust = eigencluster(sp,evecs,ii);
meanoccurance = Array{Float64}(length(eclust))
for i=1:length(eclust)
    clustersp = Array{String}(gensp[eclust[i]]);
    speciesinz = intersect(clustersp,zsp);
    zsppos = sort(indexin(speciesinz,zsp));
    covalue = com[zsppos,zsppos];
    meanoccurance[i] = mean(UpperTriangular(covalue));
end
MO[ii] = meanoccurance;
end
