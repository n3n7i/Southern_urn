

function revtarg(targvec, xtarg)

  iden = collect(1:length(targvec));

  return iden[vec(targvec .== xtarg)];

  end;


function scoretargs(targw, targen, xtarg)

  n = length(targen);

  xscore = 1e6;

  xid = 0;

  if( any(x-> x>0, targw[xtarg, .!targen]))

    for iter in 1:n

      if(!targen[iter] && targw[xtarg, iter] > 0)

        if(targw[xtarg, iter] < xscore)

          xscore = targw[xtarg, iter];

          xid = iter;

          end;

        end;

      end;

    end;

  return xid;

  end;



function logtune(targvecC, init, xscore)

  targvec = copy(targvecC);

  n = size(targvec, 1);

  retvec = zeros(Int, n);

  compvec = falses(n);

  xstop = false;

  xstep = init;

  xskip = 1;

  xiter = 1

  ziter = 1


  while (!xstop)

##    sleep(0.01);

    print("\b\b\b\b\b\b\b\b\b\b\b\b $xiter $ziter");

    xiter += 1;

    if(!compvec[targvec[xstep]])

      retvec[xstep] = targvec[xstep];

      compvec[xstep] = true;

      xstep = targvec[xstep];

      xskip = 1;

      ziter += 1;

      continue;

      end;

    zstep = revtarg(targvec, xstep);

    zstep = zstep[compvec[zstep].==false];

    if(length(zstep) == 0) 

      targvec[xstep] = scoretargs(xscore, compvec, xstep);

      if(targvec[xstep] == 0)

        xstop = true;

        end;

      end;

    if(length(zstep) == 1) 

      targvec[xstep] = zstep[1];

      end;

    if(length(zstep) > 1)

      println("skipopt!?")

      targvec[xstep] = zstep[xskip];

      if(xskip > length(zstep))

        xstop = true;

        end;

      xskip += 1; 

      end;

    end;

  return retvec;

  end;  


function xscoretune(jumpvec, s1, s2, init)

  if(jumpvec[init] == 0)

    println("Bad init!");

    end;

  zx = repDeep(jumpvec, init, length(jumpvec));

  

  res = zeros(Int, zx);

  ##compvec = falses(n);

  sc1 = 0;

  sc2 = 0;

  xstop = false;

  xstep = init;

  xiter = 1

  backs = "\b"^20

  nxstep = 0;

  while (!xstop)

    ##sleep(0.01);

    print("$backs $xiter $(round(sc1, sigdigits=4)) $(round(sc2, sigdigits=4))");

    res[xiter] = xstep;

    xiter += 1;

    ##if()

    nxstep = jumpvec[xstep];

    if(nxstep == 0 || xiter == zx)

      println("Xstop!!");
 
      xstop = true;

      continue;

      end;

    sc1 += s1[xstep, nxstep];

    sc2 += s2[xstep, nxstep];

    xstep = nxstep;

    end;

  return sc1, sc2, res;

  end;



function logtune2(targvecC, init, xscore)

  targvec = copy(targvecC);

  n = size(targvec, 1);

  retvec = zeros(Int, n);

  compvec = falses(n);

  xstop = false;

  xstep = init;

  xskip = 1;

  xiter = 1

  ziter = 1


  while (!xstop)

##    sleep(0.01);

    print("\b\b\b\b\b\b\b\b\b\b\b\b $xiter $ziter");

    xiter += 1;

    if(!compvec[targvec[xstep]])

      retvec[xstep] = targvec[xstep];

      compvec[xstep] = true;

      xstep = targvec[xstep];

      xskip = 1;

      ziter += 1;

      continue;

      end;

##    zstep = revtarg(targvec, xstep);

##    zstep = zstep[compvec[zstep].==false];

##    if(length(zstep) == 0) 

      targvec[xstep] = scoretargs(xscore, compvec, xstep);

      if(targvec[xstep] == 0)

        xstop = true;

        end;

##      end;

##    if(length(zstep) == 1) 

##      targvec[xstep] = zstep[1];

##      end;

##    if(length(zstep) > 1)

##      println("skipopt!?")

##      targvec[xstep] = zstep[xskip];

##      if(xskip > length(zstep))

##        xstop = true;

##        end;

##      xskip += 1; 

##      end;

    end;

  return retvec;

  end;  


