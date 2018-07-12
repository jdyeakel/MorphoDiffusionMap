using(DataFrames)
using(CSV)
using(RCall)
include("$(homedir())/Dropbox/Postdoc/2018_eigenvec/src/laplacian.jl")
include("$(homedir())/Dropbox/Postdoc/2018_eigenvec/src/eigencluster.jl")


zdata = CSV.read("$(homedir())/Dropbox/Postdoc/2018_eigenvec/Zliobaitedata/Zliobaitedata.csv",header=true);
zsite = CSV.read("$(homedir())/Dropbox/Postdoc/2018_eigenvec/Zliobaitedata/Zliobaitesite.csv",header=true);

sp = Array(zdata[:species]);
nsp = length(sp);
morphdata = Array(zdata[:,3:9]);
nmeas = size(morphdata)[2];
padata = Array(zdata[:,10:22]);

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
    ctones = 0;
    for j = 1:nmeas
        if morphdata[a,j] == morphdata[b,j]
            # ct += (sqrt((pcdatatr[a,j]-pcdatatr[b,j])^2)); #/mean(pcdatatr[!ismissing.(pcdatatr[:,j]),j]);
            ct += 1;
        end
        ctones += 1;
        # ctones += PA[a,j] + PA[b,j];
    end
    PC[a,b] = Float64(ct/ctones); #/Float64(ctones);
end

S = laplacian(PC,10);

ev = eigs(S; nev=10,which=:SR);
eval = ev[1];
evecs = ev[2];

ranked = sortperm(evecs[:,2]);



R"""
plot($(evecs[:,3]),$(evecs[:,4]))
text($(evecs[:,3]),$(evecs[:,4]),$sp,cex=0.5)
"""

R"""
library(plot3D)
scatter3D($(evecs[:,3]),$(evecs[:,4]),$(evecs[:,5]))
text3D($(evecs[:,3]),$(evecs[:,4]),$(evecs[:,5]),labels=$sp,cex=0.5,add=T)
"""
    
