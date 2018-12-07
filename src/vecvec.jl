

mutable struct path;

  nodes::Array{Int};

  length::Int;

  score::Float64;

  cost::Float64;

  end;


function path(node::Int)

  path([node], 1, 0, 0);

  end;


function addnode(xpath::path, node::Int, score, cost)

  push!(xpath.nodes, node);

  xpath.length +=1;

  xpath.score += score;
   
  xpath.cost += cost;
 
  end;


function joinpath(head::path, tail::path, score, cost, rev=false)

  if(!rev) append!(head.nodes, tail.nodes); end;

  if(rev) append!(head.nodes, reverse(tail.nodes)); end;

  head.score += tail.score + score;

  head.cost += tail.cost + cost;

  end;



function pathends(xpath::path)

  return [xpath.nodes[1]; xpath.nodes[end]];

  end;



function pathinit(pathvec, distvec, skipdist = 100)

  patharr = path[];  

  push!(patharr, path(pathvec[1]));

  pathid = 1;

  n = length(pathvec);

  xpath = patharr[pathid];

  for iter in 1:n-1

    xdec = distvec[iter] >= skipdist;

    if(xdec) 

      push!(patharr, path(pathvec[iter+1]));

      pathid += 1;

      xpath = patharr[pathid];

      else

      addnode(xpath, pathvec[iter+1], distvec[iter], 0);

      end;

    end;

  return patharr;

  end;


function retrecon(xn, hvec, tvec, xid)

  i = 0;

  for iter in 1:xn

    i+= !hvec[iter];

    if(xid == i)

      return :h, iter;

      end;

    end;

  for iter in 1:xn

    i+= !tvec[iter];

    if(xid == i)

      return :t, iter;

      end;

    end;

  return :x, 0;

  end;
  

function retreconT(xn, tvec, xid)

  i = 0;

  for iter in 1:xn

    i+= !tvec[iter];

    if(xid == i)

      return :t, iter;

      end;

    end;

  return :x, 0;

  end;



function pathrecon(paths, data)

  n = length(paths);

  hcomp = falses(n);

  tcomp = falses(n);

  tcomp[1] = true;

  xjoins = zeros(Int, n,4);

  pathjump = 1;

  for iter in 1:n-1

    pnow = pathends(paths[pathjump])[1]; ## [ !hcomp[iter] ]!tcomp[iter]])

    if(length(pnow)==0) println("loop early!"); continue; end;

    targs = map(x->x.nodes[end], paths[.!tcomp]); ##map(x->x.nodes[1], paths[.!hcomp]); 

    print(length(targs));

    tdists = xL6_dist(data[pnow[1], :], data[targs, :]);

    x = findmin(sfx(tdists))[2][1];

    x2 = retreconT(n, tcomp, x);

    print(" >", x, " ", x2, "<  ")


    if(x2[1] !== :x && x2[2]!==pathjump)

      xjoins[iter, :] = [ pnow[1] targs[x] pathjump x2[2]];

      hcomp[pathjump] = true;

      tcomp[x2[2]] = true;

      pathjump = x2[2];

      end;

    println(x, " ", x2[1], " ", x2[2], " ", iter);

    end;

  return xjoins;

  end;


function gentour(xpaths, xjoins)

  tourinit = path(-1);

  n = length(xpaths);  

  for iter = n-1:-1:1

    joinpath(tourinit, xpaths[xjoins[iter, 4]],  0, 0);

    end;

  joinpath(tourinit, xpaths[1], 0, 0);
  
  return tourinit;

  end;  


function tourmissed(xtvec, xdata)

  n = size(xdata,1);

  xmiss = falses(n)

  for iter in xtvec

    xmiss[iter] = true;

    end;

  return collect(1:n)[.!xmiss];

  end;

