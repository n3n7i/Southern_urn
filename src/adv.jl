

using LinearAlgebra


function advDist(dataset, init=1)

  n = size(dataset,1)

  xdists = zeros(n,n);

  xcertain = zeros(n);

  xbase = dataset[init,:];

  xdists[:, init] = xL5_dist(dataset[init, :], dataset);

  for i in 1:n

    xtarg = dataset[i, :];

    zc = xdists[i, init]; ##xL5_dist(xbase, xtarg)

    for j in 1:n
 
      xdists[j, i] = xdists[j, init] + zc;

      end;

    xcertain[i] = zc;

    end;

  return (xdists, xcertain);

  end;


function triDist(dataset, init=1, xdef = 3, dfun = xL5_dist)

  n = size(dataset,1)

  xdists = zeros(n,xdef);

  xcertain = zeros(n);

  xbase = dataset[init,:];

  prexdist = dfun(dataset[init, :], dataset);

  x2 = maximum(prexdist);

  x3 = mean(prexdist);

  x4 = (x2+x3)/2;

  idvec = collect(1:n);

  for iter in 2:xdef

    idx = [idvec[vec( ( (x4*((iter-1)/xdef)) .< prexdist) .& (prexdist .< (x4*((iter)/xdef)) ) )]; rand(1:n)][1];

    xdists[:,iter] = dfun(dataset[idx, :], dataset);

    end;

  ##id2 = collect(1:n)[vec(prexdist .== x2)][1];

  ##id3 = collect(1:n)[vec( ((x3*(2/3)) .< prexdist) .& (prexdist .< (x3*(3/2))) )][1];

  xdists[:, 1] = prexdist;

##  xdists[:,2] = xL5_dist(dataset[id2, :], dataset);

##  xdists[:,3] = xL5_dist(dataset[id3, :], dataset);
  
  return xdists;

  end;


function getEst(td, id1, id2)

  return maximum(abs.(td[:, id1] .- td[:, id2]));

  end;


function getEstB(td, id1, id2)

  m = size(td,1);

  xm = 0;

  for iter=1:m

    tval = abs(td[iter,id1] - td[iter,id2]);

    (xm < tval) && (xm = tval);

    end;

  return xm;

  end;



function getEsts(tD, ids)

  n = size(ids,1);

  r = zeros(n,n);

  tDx = @view tD[:, ids];

  @inbounds @simd for i in 1:n

    for j in 1:n

      r[j,i] = getEst(tDx, i, j); ##ids[i], ids[j]);

      end;

    end;

  return r;

  end;


function getEstsB(tD, ids)

  n = size(ids,1);

  r = zeros(n,n);

  tDx = @view tD[:, ids];

  @inbounds @simd for i in 1:n

    for j in 1:n

      r[j,i] = getEstB(tDx, i, j); ##ids[i], ids[j]);

      end;

    end;

  return r;

  end;


function getEstC(tD, ids, targ, cutoff=1000)

  n = size(ids,1);

  r = zeros(n);

  b = BitArray(undef, n);

  b.= true;

  wid = size(tD,1);

  tDx = @view tD[:, ids];

  tGx = @view tD[:, targ];

  for j in 1:wid

    tvec = collect(1:n)[b];

    tDx2 = @view tDx[j, :];

    for i in tvec ##  @inbounds @simd

      tval = abs(tDx2[i] - tGx[wid]); ##tDx2[targ];

      (r[i] < tval) && (r[i] = tval); ##tDx2[i] - tDx2[targ] 

      b[i] = tval < cutoff;

      end;

    end;

  return (r, b);

  end;


function getEstsC(tD, ids, targ, cutoff=1000)

  m,n = length(ids), length(targ); 

  r = zeros(m,n);
 
  b = BitArray(undef, m,n);

  print(n, "* ")

  for iter=1:n

    r[:,iter], b[:,iter] = getEstC(tD, ids, targ[iter], cutoff);

    end;

  return (r, b);

  end;



function getEstsCh(tD, ids, targ, cutoff=1000, xhalt=100, rmax=xhalt)

  n = length(targ);

  if(n > xhalt) n = xhalt; end;

  r,b = getEstsC(tD, ids, targ[1:n], cutoff);

  b2 = any(b, dims=2);

  print(size(b2), " ", sum(b2), " | |  ");

  if(sum(b2) > rmax)

    r2 = minimum(r[vec(b2), :], dims=2);

    println(sum(r2.==0), " - nil");

    p2 = sortperm( r2[:] );

    ##b3 = @view b2[vec(b2)];

    @yulio_Do b2[vec(b2)][ p2[rmax+1:end] ] .= false;

    ##b3[ p2[rmax+1:end] ] .= false;

    println("\n $(size(b2)), $(size(r2)), $(size(p2)) XX");

    println(">>", sum(b2), "<<  ");

##  if()

##  b[b][]

    end;

  return (r,b,b2); ##getEstsC(tD, ids, targ[1:n], cutoff);

  end;





function getEsts2(tD, ids)

  n = size(ids,1);

  w = size(tD, 1);

  r = zeros(n,n);

  t = zeros(w, n);

  tDx = @view tD[:, ids];  

  @inbounds @simd for i in 1:n

    t[:,:] .= @view tDx[:, i];

    BLAS.axpy!(-1, tDx, t);

    r[:,i] = maximum(abs.(t),dims=1);

##    ##for j in 1:n

##      r[j,i] = getEst(tDx, i, j); ##ids[i], ids[j]);

##      end;

    end;

  return r;

  end;



