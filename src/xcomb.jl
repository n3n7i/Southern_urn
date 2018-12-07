
function pointinsert(pathnodes, joinnodes, xdata)

  m = size(joinnodes,1);

  println(m);

  n = size(pathnodes,1);

  xdists = zeros(n, m);

  pathB = copy(pathnodes)

  for iter in 1:m

    xdists[:,iter] = xL6_dist(xdata[joinnodes[iter], :], xdata[pathnodes, :]);

    end;

##  return xdists;

  zdists = cartmap( findmin(sfx(xdists), dims=1)[2], 1); 

##  show(zdists);

  for iter in 1:m

    insert!(pathB, zdists[iter], joinnodes[iter]);

    end;

  return pathB;

  end;

