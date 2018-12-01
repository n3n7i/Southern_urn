
function santaScore1(data, Pathids, xprimes)

  xscore = 0;

  xstep = 0

  for iter = 2:length(Pathids)

    xstep +=1;

    stepc = pairwise(Euclidean(), reshape(data[Pathids[iter-1], :],2,1), reshape(data[Pathids[iter], :],2,1));

    stepc2 = 1.0;

    if( (xstep % 10 == 0) && (!any(x->x == Pathids[iter-1], xprimes)) ) stepc2 = 1.1; end;

    xscore += stepc[1] * stepc2;

    if(iter%10000 == 0) print(iter, " "); end;

    end;

  return xscore;

  end;



function santaScore2(data, Pathids, xprimes)

  xscore = 0;

  xstep = 0

  datax = @view data[Pathids, :];

  @time dists = sqrt.(sum((datax[1:end-1, :] .- datax[2:end, :]).^2, dims=2));

  @time for iter = 10:10:length(Pathids)-1

    if( !any(x->x == Pathids[iter-1], xprimes) ) dists[iter] *= 1.1; end;  

    end;

  return sum(dists);

  end;


function santaScore3(data, Pathids, xprimes)

  xscore = 0;

  xstep = 0;

  xp = BitArray(undef, maximum(xprimes));

  xp .= false;

  xp[xprimes] .= true;

  datax = @view data[Pathids, :];

  @time dists = sqrt.(sum((datax[1:end-1, :] .- datax[2:end, :]).^2, dims=2));

  @time for iter = 10:10:length(Pathids)-1

    if( !xp[Pathids[iter]-1] ) dists[iter] *= 1.1; end;  

    end;

  return sum(dists);

  end;


