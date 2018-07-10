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


pcdata = CSV.read("$(homedir())/Dropbox/Postdoc/2018_eigenvec/fieldbirds.csv",header=true);
# sp = pcdata[:common_name];
genus = pcdata[:Genus];
species = pcdata[:species];
gensp = mapslices(join, [Array{String}(genus) Array{String}(species)], 2);
nsp = length(unique(gensp));
nmeas = 
#Find averages for values across species
#Across genera
for i=1:nsp
    inds = find(x->x==gensp[i],gensp);
    
