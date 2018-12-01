##run(`ls ./`) #`../input`


function getstart()
  da = readcsv("../input/gifts.csv");
  data = map(Float64, da[2:end,:]);
  #data2 = mfb(data);
  println(da[1,:]," ", size(data,1));
  return (data); #,data2);
  end;


c_Earth = 6371;

function deg2rad(ang)
  return (ang ./ (180 / pi));
  end;

function havers(lat1, lng1, lat2, lng2)
  lat1 = deg2rad(lat1);
  lat2 = deg2rad(lat2);
  lng1 = deg2rad(lng1);
  lng2 = deg2rad(lng2);
  lat = lat2 .- lat1;
  lng = lng2 .- lng1;
  d = sin(lat ./ 2) .^ 2 .+ cos(lat1) .* cos(lat2) .* sin(lng ./ 2) .^ 2;
  #println("$([d,  sin(lat / 2), cos(lat1), cos(lat2), sin(lng / 2)])");
  h = 2 * c_Earth .* asin(sqrt(d));
  return (h);
  end;

function accsum(fvec, initw = 10)
  n = length(fvec);
  r = zeros(n+1);
  r[1] = initw;
  for iter = 1:n
    r[iter+1] = r[iter] + fvec[iter];
    end;
  return (r);
  end;

function autohav(dvec)
  home = [90  0];
  dv2 = [home; dvec; home];
  dis = havers(dv2[1:end-1,1], dv2[1:end-1,2], dv2[2:end,1], dv2[2:end,2]);
  return dis;
  end;

function getscore(res, da)
  n = maximum(res[:,2]);
  sc = 0;
  xx = 0;
  for iter = 1:n
    tvec = res[res[:,2] .== iter, 1];
    ow = accsum(da[tvec, 4]);
    od = autohav(da[tvec, 2:3]);
    sc += sum(ow .* od);
    xx += ow[end] > 1010;
    end;
  print(xx, ":oversize ");
  return(sc);
  end;


## cluster imp.

function xS_cluster(data, centers, iter, fxSd = xL5_dist)
n = size(data, 1);
m = size(data, 2);
x = randperm(n);

cents, dnow = nlnmodeB(data, centers, mean(data, dims=1));

println(size(dnow), size(data), size(cents));

##dnow = ones(n, 1) * 1e25;

dnow = dnow .* 2.0; 

inow = zeros(n,1);
oldcents = similar(cents);
centchange = 10;
for counti = 1:iter

  print("$counti ");
  if((centchange < 1e-4) & (counti > iter/2))
    continue;
    end;

  oldcents[:] .= cents[:];
  for countc = 1:centers

    dist = fxSd(cents[countc, :], data);
    ivec = dnow .> dist;
    xvec = inow .== countc;

    rvec = dotDollar(ivec, xvec);
    dnow[rvec] = dist[rvec];

##    print(typeof(rvec));
    inow[rvec] .= countc;

    #dnow[xvec] = dist[xvec];

    ##xindex = logic2index(inow .== countc);

##    print(size(xvec), size(data));

##    print(typeof(xvec), typeof(data));

    datax = data[xvec[:], :];  ##data[xindex,:];

    if(((counti == 1) + (counti == iter)) == 0) ##skip first & last update?

##      print(size(datax), size(cents), size(countc));

      cents[countc,:] = mean([datax; cents[countc,:]'], dims=1);

      end;
    ##print(".");
    end;

  centchange = sum(abs.(cents[:] .- oldcents[:]));

  print("$(centchange) $counti | ");
  end;
#if(mean(dnow) < 0)
  #print("\n simscore = $(mean(dnow) / centers), mean sim = $(mean(dnow))\n");
#else
  #print("\n distscore = $(mean(dnow) * centers), mean dist = $(mean(dnow))\n");
  #end;
return (cents, inow);
end;

function maxind(datavec)
i = 1;
for iter = 2:length(datavec)
  i = datavec[i] > datavec[iter] ? i : iter;
  end;
return i;
end;

function nlnmodeB(data, citer, cent1)
n = size(data,1);
cents = zeros(citer+1, size(data,2));
cents[1,:] = cent1;
dnow = ones(n, 1) * 1e25;
for iter = 2:citer+1
  dist = xL5_dist(cents[iter-1, :], data);
##  print("$(size(dnow)) $(size(dist)) dd\n")
  ivec = dnow .> dist;
  dnow[ivec] = dist[ivec];
  nlnid = maxind(dnow);
  cents[iter,:] = data[nlnid, :];
  end;

println(size(dnow), "?");

return (cents[2:end, :][sortperm(cents[2:end,1]),:], dnow);
end;

function dotDollar(vec1, vec2)
#n = length(vec1);
#v3 = falses(n);
#for(iter = 1:n)
  #v3 = vec1[iter] $ vec2[iter];
  #end;
return (vec1 .+ vec2) .> 0;
end;


function xL5_dist(cent, data)
##  print(size(cent), size(data));

  return sum(abs.(data .- cent'), dims=2);
  end;


function logic2index(logvec::BitArray)
n = sum(logvec .== true);
i = zeros(Int64, n,1);
j=0;
##println("length $n")
for iter = 1:size(logvec,1) 
  if(logvec[iter] == true)
    j = j+1;
    i[j] = iter;
    ##println("$iter $j");
    end;
  end;
return (i[:]);
end;

## eoimp.

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

###########

## All done! (Optimizer 1)

##########

## Optimizer 2 (Run North Heuristic)

function reopt(ox, d1, thres = 0.05)
  n = maximum(ox[:,2]);
  oxb = similar(ox);
  i = 0;
  for iter = 1:n
    cvec = ox[:,2] .== iter;
    cpath = ox[cvec, :];
    clist = d1[cpath[:,1][:], :];
    if(sum(cvec) == 0)
      continue;
      end;
    cwei = sum(clist[:,4]);
    cwrem = 1000 - cwei;
    dvec = ox[:,2] .> iter;
    dpath = ox[dvec, :];
    dlist = d1[dpath[:,1], :];
    evec = abs(dlist[:,3] .- clist[end,3]) .<= thres;
    elist = dlist[evec, :];
    elist = elist[sortperm(elist[:, 2]), :];
    xlist = Int[];
    j = 1;
    while( (cwrem > 1) + (j <= size(elist, 1)) == 2)
      if (elist[j, 4] <= cwrem)
        push!(xlist, elist[j, 1]);
        cwrem -= elist[j, 4];
        end;
      j += 1;
      end;
    zlist = similar(xlist);
    zlist[:] = iter;
    for iterb = 1:size(xlist,1)
      ox[ox[:,1] .== xlist[iterb], 2] = 0;
      end;
    npath = [cpath, [xlist zlist]];
    oxb[i+1:i+size(npath,1), :] = npath;
    i+= size(npath,1);
    print(" ", iter);
    end;
  return (oxb);
  end;

##########


##d1 = getstart();

##ox = getsolve(d1, 100, 10, 3.5, 6); ##using ok-ish(?) settings

##score = getscore(ox, d1);

##writecsv("subx.csv", [["GiftId" "TripId"]; ox[end:-1:1, :]]);

##print("\n\n Submission Scores ", score);

## Final checks

##k = zeros(100000);

##for iter = 1:100000

##  k[ox[iter, 1]] += 1;
  
##  end;
  
##kk = sum(k.==1) == 100000;

##print("  submission ok : ", kk);


##ox2 = similar(ox);

##ox2[:] = ox[:];

##ox3 = reopt(ox2, d1, 0.12);


##score = getscore(ox3, d1);

##writecsv("subxb.csv", [["GiftId" "TripId"]; ox3[end:-1:1, :]]);

##print("\n\n Submission Scores ", score);


