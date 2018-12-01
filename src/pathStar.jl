

function pathStar1(xdata, xids, xinit, exitdep=10, xsense = 250, xlim = 25, xret = 50, xprec = 4)

  n = size(xids, 1);

  b0 = BitArray(undef, n); ##search

  b1 = BitArray(undef, n); ##complete

  b0 .= false; b0[xinit] .= true;

  b1 .= false; #b0[xinit] = true;

  evec = copy(triDist(xdata[xids, :], xinit[1], xprec, xL6_dist)'); ##xL6_dist(xdata[xinit, :], xdata)');

  for iter in 1:exitdep    

    xvec = b0 .& (.!b1);

    print(iter, " vsum:", sum(xvec), "  ");

    (sum(xvec)==0) && continue;

    ##for wid in xids[xvec]

    zd, zb, zb2 = getEstsCh(evec, 1:n, xids[vec(xvec)], xsense, xlim, xret);

    ##xvec[xvec] = zb;

    b1 = b1 .| b0;

    b0 = zb2; ##any(zb, dims=2);     

    println("Loop, ", sum(b0), " ", sum(b1));
    
    end;

  return b1;

  end;



function pathStar2(subdata)

  

  end;
