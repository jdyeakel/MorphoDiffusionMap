
using(DataFrames)
using(CSV)

# df = CSV.read("$(homedir())/Dropbox/PostDoc/2018_foodwebs/robertson_1929_matr.csv",header=false);

df = CSV.read("$(homedir())/robertson_1929_matr.csv",header=false);
pols = CSV.read("$(homedir())/robertson_1929_pol.csv",header=false);

#preprocessing
ints = find(!iszero,sum(Array(df),2));
df = df[ints,:];
pols = pols[ints,:];
npol = size(df)[1];
nplant = size(df)[2];

A = SharedArray{Float64}(npol,npol);
#Build similarity matrix
@time @sync @parallel for i = 0:(npol^2 - 1)
    a = mod(i,npol) + 1;
    b = Int64(floor(i/npol)) + 1;
    
    if a == b
        A[a,b] = 0.0;
        continue
    end
    # if a==1
    #     println(b)
    # end
    ct = 0;
    ctones = 0;
    for j = 1:nplant
        ct += (df[a,j]*df[b,j]); #+ (1 - df[a,j])*(1 - df[b,j])
        ctones += df[a,j] + df[b,j];
    end
    A[a,b] = Float64(ct)/Float64(ctones);
    
end

S = zeros(Float64,npol,npol);
for i = 1:npol
    
    val = zeros(Float64,10);
    loc = zeros(Int64,10);
    
    for j = 1:npol
        if A[i,j] > val[10]
            val[10] = A[i,j];
            loc[10] = j;
        end
        for k=1:9
            if val[11-k] > val[10-k]
                v = val[11-k];
                val[11-k] = val[10-k];
                val[10-k] = v;
                l = loc[11-k];
                loc[11-k] = loc[10-k];
                loc[10-k] = l;
            end
        end
    end
    S[i,loc] = -A[i,loc];
    S[loc,i] = -A[loc,i];
end
rowsums = sum(S,2);
S[diagind(S)] = -rowsums;

ev = eigs(S; nev=10,which=:SR);
eval = ev[1];
evecs = ev[2];

set1 = find(x->x>0,evecs[:,2]);
set2 = find(x->x<0,evecs[:,2]);
pols[set1,:]
pols[set2,:]

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
