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
        if CM[a,j] == CM[b,j]
            # ct += (sqrt((pcdatatr[a,j]-pcdatatr[b,j])^2)); #/mean(pcdatatr[!ismissing.(pcdatatr[:,j]),j]);
            ct += 1;
        end
        ctones += 1;
        # ctones += PA[a,j] + PA[b,j];
    end
    PC[a,b] = Float64(ct/ctones); #/Float64(ctones);
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