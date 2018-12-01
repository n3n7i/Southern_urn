
## Swarm Tree imp.

### KNN adjacency matrix  (dead param rx)
function ws_adj(items, rx = 0.1, px = 10, xf = xL5_dist)
  if(px > size(items,1))
    px = size(items,1);
    end;
  n = size(items,1);
  adj = falses(n, n);
  adj2 = zeros(n);
  for iter = 1:n
    adj2 = xf(items[iter,:], items);
    rp = sortperm(adj2[:]);
#    rm = minimum(adj2);
  #  println("$(rp[1:5])");
#    adj[iter,:] = adj2 .<= (rm + (maximum(adj2) - rm) * rx);
    adj[iter, rp[1:px]] = true; 
    adj[iter, iter] = false;
    end;
  return ((adj .+ adj') .> 0);
  end;


### Swarm Core

function unigrap2(wdata, adj, iid, iw = 10, lim = 990)
n = size(wdata,1);
cvec = trues(n);
c2vec = falses(n);
p = zeros(n);
w = ones(n) .* 1e7;
s = ones(n) .* 1e12;
t = zeros(n);
r = zeros(n);
iter = 1;
p[iid] = iter;
w[iid] = wdata[iid] .+ iw;
s[iid] = w[iid];
r[iid] = 1:length(iid); ## Ok, iid is a vector
nv = 1:n;
subiter = 1;
c2vec[1] = true
while (sum(c2vec) > 0)
  iter += 1;
  c2vec[:] = false;
  subiter = 0;
  for itern = iid
    subiter +=1;
    xvec = adj[itern, :][:];
  #  print(size(xvec));
    qval = wdata[xvec] .+ w[itern];
    q2val = qval .+ s[itern];
    nvec = (s[xvec] .> q2val) .+ (p[xvec] .!= 1) .+ (w[xvec] .> qval) .== 3;
    xvec[xvec] = nvec;
  #  print("nv-", sum(nvec), " ");
    s[xvec] = q2val[nvec];
    w[xvec] = qval[nvec]; ##weight integ
    p[xvec] = iter; ##depth 
    t[xvec] = itern; ##prev
    r[xvec] = r[itern]; ##root
    c2vec[xvec] = true;
   # print("c2-",sum(c2vec), " ");
    end;
##  print(subiter, " ");
  cvec[iid] = false;
  iid = nv[c2vec];
  if(sum(c2vec) == 0)
##    println("\n:", iter);
    iter = lim +1;
    end;
#  print(iid);
  end;
return([s w p t r]);
end;


## Swarm data Extraction
## really needs a cleanup

function grapext(gstat)
n = size(gstat,1);
gpath = round(Int64, maximum(gstat[:,5]));
gdepth = round(Int64, maximum(gstat[:,3]));
gtrace = falses(gpath, gdepth);
gmax = zeros(Int64, gpath);
nvec = 1:n;
for iter = 1:gdepth
  v = gstat[:,3] .== iter;
  gtrace[gstat[v,5], iter] = true;
  gmax[gstat[v,5]] = iter;
  end;
pmat = zeros(Int64, gpath, gdepth);
wvec = zeros(gpath);
for iter = 1:gpath
  v = ((gstat[:,3] .== gmax[iter]) .+ (gstat[:,5] .== iter)) .== 2;
  v[v] = gstat[v,1] .== minimum(gstat[v,1]);
  pmat[iter, gmax[iter]] = (nvec[v])[1];
  wvec[iter] = gstat[ (nvec[v])[1], 2];
  end;
for iter = gdepth:-1:2
  v = pmat[:, iter] .!= 0;
  pmat[v, iter-1] = gstat[pmat[v, iter], 4];
  end;
pcheck = pmat[:][pmat[:] .!=0];
pc = 0;
for iter = 1:size(pcheck,1)-1
  pc += sum(pcheck .== pcheck[iter]) -1;
  end;
if(pc > 0)
  println("Grapext Error!!");
  #sleep(15);
  end;
return (pmat, gmax, wvec);
end;


## Randomized start

function initwnodes(xdat, ii)
  nvec = 1:size(xdat,1);
  ivec = zeros(Int64, ii);
  ik = minimum([size(xdat,1), ii])
  (a3,b3) = xS_cluster(xdat[:,2:3], ik, 20, xL5_dist);
  for iter = 1:ik
    v = b3[:] .== iter;
    v[v] = xdat[v,4] .== minimum(xdat[v,4]);
    ivec[iter] = (nvec[v])[1];
    end;
  return(ivec);
  end;

## Subselection helper 

function submap(v, v2)
  n = size(v,1);
  k = 0;
  xv = zeros(Int64, n);
  for iter = 1:n
    k += v[iter];
    xv[iter] = k;
    end;
  return(xv[v2], xv[v]);
  end;


## Swarm Main

function ws_main(xdat, group1 = 10, wlim = 990, blim = 50, zlim = 1.1, zpred = 6, zlim2 = 10)
n = size(xdat,1);
hvec = 1:n;
n1 = round(Int64, ceil((sum(xdat[:,4]) / wlim) * zlim));
slim = round(Int64, ceil((size(xdat,1) / n1) * zlim2));
(a2, b2) = xS_cluster(xdat[:,4], group1, 20, xL5_dist);
chord = zeros(Int64, n1, slim);
chorn = zeros(Int64, n1);
whorn = zeros(Float64, n1) .+ 10.0;
invhchor = hvec[b2[:].==1];
chornow = invhchor[ initwnodes(xdat[b2[:].==1,:], n1) ];
sdat = [];
for iter = 1:group1
  v = b2[:] .== iter;
  if(sum(v) == 0)
    continue;
    end;
  v[chornow] = true;
  v2 = v[v];
  sdat = xdat[v, :];
  (hchornow, hchor) = submap(v, chornow);
  invhchor = hvec[v];
  iterc = 0;
  gvec = 1:size(v2,1);
  while (iterc <= blim)
    v2[hchornow] = true;

#    println("\n rem-", sum(v2));

    (gchorn, gchor) = submap(v2, hchornow);
#    sdat = xdat[v, :];
    sadj = ws_adj(sdat[v2,2:3], 0, zpred);
    gstat = unigrap2(sdat[v2,4], sadj, gchorn, whorn, blim);
    (chorx, cz, wz) = grapext(gstat);
    whorn = wz;
    invgchor = gvec[v2];
    for iterb = 1:n1
      chord[iterb, chorn[iterb]+1:chorn[iterb]+cz[iterb]-1] = invhchor[ invgchor[ chorx[iterb, 1:cz[iterb]-1] ] ];
      chorn[iterb] = chorn[iterb] + cz[iterb]-1;
      hchornow[iterb] = invgchor[ chorx[iterb, cz[iterb]] ];
      end;

#    println("max-chorn : ", maximum(chorn), ", slim : ", slim, ", sum ", sum(chorx .!=0));

    v2[invgchor[ chorx[:][chorx[:].!=0] ]] = false;
    iterc += 1;
#    if(iterc == blim)
#      println("Overrun!");
#      sleep(15);
#      #break;
#      end;
    if(sum(v2) == 0)
      iterc += blim;
      end;
    whorn = whorn .- sdat[hchornow, 4];
    end;
  chornow = invhchor[ hchornow ];
  #?
  end;
whorn = whorn .+ xdat[chornow, 4];
#println("zlim2 ");
for iter = 1:n1
  chord[iter, chorn[iter]+1] = chornow[iter];
  end;
chorn = chorn .+1;
print("x");
return(chord, chorn, whorn);
end;


## Swarm upper

function getsolve(data3, xc = 100, zc = 10, sparesol = 3.0, fanout = 6)
  n = size(data3,1);
  nv = 1:n;
  outx = zeros(Int64, n, 2);
  (a,b) = xS_cluster(data3[:, 2:3], xc, 10, xL5_dist);
  rc = 0;
  ix = 0;
  r = 1;
  for iter = 1:xc
  
    print("!", iter);
  
    v = b[:] .== iter;
    if(sum(v) == 0)
      continue;
      end;
    rv = nv[v];
    (r1,r2,r3) = ws_main(data3[b[:].==iter, :], zc, 1000, 25, sparesol, fanout);
##
pcheck = r1[:][r1[:] .!=0];
pc = 0;
for iter = 1:size(pcheck,1)
  pc += sum(pcheck .== pcheck[iter]) -1;
  end;
if(pc > 0)
  println("Santz Error!!");
  sleep(15);
  end;
##
    for iterb = 1:size(r1,1)
      if((ix+r2[iterb]) > 100000)
        print("ix !", iter);
        end;
      outx[ix+1:ix+r2[iterb],1] = rv[ r1[iterb, 1:r2[iterb]][:] ];
      outx[ix+1:ix+r2[iterb],2] = r;
      r = r+1;
      ix = ix + r2[iterb];
      end;
    rc = rc .+ sum(r3 .> 1010.0)
    end;
  println("Oversize count :: ", rc);
  return(outx);
end;
