

module xcheck

  using DelimitedFiles

  if (!isdefined(Main, :xdata))

    @time d1 = readdlm("../cities.csv", ',', Any, '\n');

    d2 = Float64.(d1[2:end, 2:3]);

    end;

  end;


if (isdefined(xcheck, :d2))

    xdata = xcheck.d2;

    end;



function xL6_dist(cent, data)
##  print(size(cent), size(data));

  return sqrt.(sum((data .- cent').^2, dims=2));
  end;


function stepscan(data, init, activeperm)

  d0 = xL6_dist(data[init, :], data[activeperm, :])

  d0s = sortperm(d0[:]);

  return d0, d0s;

  end;


function stepsprime(data, init, maxgroup, minrat)

  n = size(data,1);

  xdist, xperm = stepscan(data, init, 1:n);

  ((n > maxgroup) && (xdist[xperm[2]]*minrat > xdist[xperm[maxgroup]])) && return xperm[1:maxgroup];

  j = maxgroup;

  k = xdist[xperm[2]]*minrat;

  while ((k < xdist[xperm[j]]) && (j>20))

    j-= 10;

    end;

  return xperm[1:j]  

  end;


function wx_adj(items, px = 10, xskip = 0, xf = xL6_dist)

  if(px > size(items,1))
    px = size(items,1);
    end;

  n = size(items,1);
  adj = zeros(n, n);
  adj2 = zeros(n);

  for iter = 1:n
    adj2 = xf(items[iter,:], items);
    rp = sortperm(adj2[:]);

    adj[iter, vec(rp[1+xskip:px+xskip])] .= adj2[vec(rp[1+xskip:px+xskip])]; 

    end;
  return adj;
  end;



sfx(x) = x .+ (x.==0) .* 10000;



function wx_mincost(xadj)

  n = size(xadj, 1);

  (minA, inds)  = findmin(sfx(xadj), dims=2); ## !

  ##for iter in 1:n
 
  end;
    


function zlogtunnel2(ladj, wadj, wmin)

  n = size(ladj,1);

  wmod = similar(wadj);

  wrate = similar(wadj);

  for iter in 1:n

    xvec = vec(ladj[iter, :]);
 
    xvals = wadj[iter, :]

    wmod[iter, xvec] = xvals[xvec] .- min.(wmin[xvec], wmin[iter]) .+ 1e-4; ## !

    wrate[iter, xvec] = xvals[xvec] ./ min.(wmin[xvec], wmin[iter]);

    end;

  return [wmod, wrate];

  end;


function cartmap(xcart, def=2)

  return map(x->x[def], xcart);

  end;


function histn(xhists, xrange=0)

  n = xrange;
  
  (n == 0) && (n = length(xhists));
 
  return [sum(xhists .== i) for i in 1:n]

  end;





function testrunz(data, init, maxgroup, minrat, maxwid)

  startv = stepsprime(data, init, maxgroup, minrat)

  n = length(startv);

  xmaxwid = min(maxwid, round(Int, n/2));

  startw = wx_adj(data[startv, :], xmaxwid, 1);

  startl = startw .> 0;

  startminw, sminCart = wx_mincost(startw);

  mid1, mid2 = zlogtunnel2(startl, startw, startminw); 

  v1, v2 = wx_mincost(mid1);

  v2b = cartmap(v2);

  v2c = histn(v2b);
  
  return startv, startw, startl, mid1, mid2, v1, v2b, v2c

  end;


function swarm_tunnels(vecA, vecB)

  jumpvec = vecA[vecB .== 1];

  oldvec = similar(jumpvec);

  xdepth = 1;

  while(length(jumpvec) > 0)  

    print(xdepth, " ", length(jumpvec), " || ");

    oldvec = copy(jumpvec);

    jumpvec = vecA[jumpvec][vecB[jumpvec].==1];

    xdepth += 1;

    sleep(0.5);

    end;

  return oldvec, xdepth;

  end;



function swarm_tunnelsB!(vecA, vecB)

  n = length(vecA);

  iden = collect(1:n);

  jumpvec = vecB .> 1;

  oldvec = similar(vecA, Bool);

  oldvec .= vecA .== -1;

  ##xdepth = 1;

  ##while(length(jumpvec) > 0)  

  for iter in iden[jumpvec]

    ##print(xdepth, " ", length(jumpvec), " || ");

    ##oldvec = copy(jumpvec);

    ##jumpvec = vecA[jumpvec][vecB[jumpvec].==1];

    if(vecB[iter] >= 2 && !oldvec[iter] && vecB[vecA[iter]] >= 2) 

      vecB[vecA[iter]] -= 1;

      oldvec[iter] = true;

      vecA[iter] = -1; ## !

      end;

    end;

  ##vecB .= vecB .+ (vec(oldvec) .==0);

  ##jumpvec = vecB .==1;
  
  println(histn(vecB.+1)[1:10]);

  ##return(oldvec, vecA, vecB);

  end;



function repDeep(vecA, item, maxdep)

  xitem = item;  

  for iter in 1:maxdep

    xitem = vecA[xitem];

    if(xitem == item || xitem < 1) 

      return iter;

      end;

    end;
 
  return maxdep;

  end;



function stunnel3(vecA, maxdep)

  n = length(vecA);

  vecC = similar(vecA, Int);

  vecC .= 0;
  
  for iter in 1:n

    vecC[iter] = repDeep(vecA, iter, maxdep);

    end;

  return vecC;

  end;



