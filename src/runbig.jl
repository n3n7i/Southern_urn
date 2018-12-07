

function initgroups(xdata, xp = (centers=100, iter = 20, submean=0.6))

  r1 = xS_cluster(xdata, xp.centers, xp.iter)

  r2 = histn(r1[2], xp.centers);

  r3 = mean(r2) .* xp.submean;

  r4 = r1[1][ (r2 .> r3) , :];

  return r4;

  end;

  
function getgroups(data)

  n = size(data,1);

  xscore = zeros(n) .+ 1e7;

  xid = zeros(Int, n);

  cents = initgroups(data)

  for iter in 1:size(cents,1)

    zdist = xL6_dist(cents[iter, :], data);

    zvec = zdist .< xscore;

    xid[vec(zvec)] .= iter;

##    println("\n\n", size(zvec), " ", size(xscore), " ", size(zdist));

    xscore[vec(zvec)] .= zdist[vec(zvec)];

    end;

  return xid;

  end;



function runpass(data, init, xgroup, xrate, xwid)

  tz2 = testrunz(data, init, xgroup, xrate, xwid)

  sv2 = logtune2(tz2[7], init, tz2[4]); 

  xa, xb, xc = 0,0,[];

  if (length(sv2) > 1)

    println("sv: $(length(sv2))");

    xa, xb, xc = xscoretune(sv2, tz2[2], tz2[4], init);

    end;

  vecstr = tz2[1][xc[xc .> 0]];

  return vecstr, xa, xb;

  end;



function rungroups(data)

  n = size(data,1);

  iden = collect(1:n);

  xid = getgroups(data);

  n2 = maximum(xid);

  xret = [];

  xscores = zeros(n2,2);

  for iter in 1:n2

    xnd = iden[xid .== iter];    

    xne, s1, s2 = runpass(data[xnd, :], 1, length(xnd), 150, 500);

    append!(xret, xnd[xne]);

    xscores[iter, :] = [s1, s2];

    end;    

  return xret, xscores;

  end;

