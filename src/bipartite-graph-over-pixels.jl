#using DataFrames

# Load in the data
#train = readtable("../input/train.csv")

#println("The columns in the train set are:\n")
#println(names(train))
#println(@sprintf("\nThere are %d rows in the training set", nrow(train)))

# It's yours to take from here!

## 
function loadDat(fname, opt = false)
  xred = readcsv(fname);
  xgreen = xred[1,:];
  xblue = xred[2:end, :];
  print("$fname $(size(xblue)) \n $xgreen \n");
  if(opt) 
    return (xblue, xgreen);
  else
    return xblue;
    end;
  end;

function logmaker(xvec, xwid)
n = size(xvec,1);
lm = zeros(n, length(xwid));
for(iter = 1:length(xwid))
  lm[:, iter] = xvec .== xwid[iter];
  end;
return lm;
end;


function wlogsim(uvec, cmat)
  return sum(min(repmat(uvec, size(cmat,1)),cmat), 2);
  end;

function wxlogsim(uvec, cmat)
  return sum(min(repmat(uvec, size(cmat,1)),cmat), 2) ./ sum(uvec);
  end;


function weightuc(loguc, clog, ulog)
cs = size(clog,1);
us = size(ulog, 1);
wuc = zeros(Float32, cs, us);
for(iter = 1:us)
  xvec = loguc[:, iter];
  wuc[xvec, iter] = wlogsim(ulog[iter,:], clog[xvec, :]);
  end;
return wuc;
end;

function weightedmean(inpmat, wvec, dim = 1)
return sum(inpmat .* wvec, dim) ./ (sum(wvec) .+ 1e-5);
end;

 
function wlogtunnelcw(loguc::BitArray, clog, wlog)

cwid = size(loguc,1);
nu = size(loguc,2);
cfeat = size(clog,2);
u_clog = zeros(Float32, nu, cfeat);

tic();
print("$(size(loguc)), $(size(clog)), > $(size(u_clog))\n");
print("$(cwid == size(clog,1)) !\n");
print(size(clog), size(wlog), size(loguc[:,1]));

for (iter = 1:nu)
  xvec = loguc[:,iter];
  u_clog[iter, :] = weightedmean(clog[xvec, :], wlog[xvec, iter], 1);
  end;

toc()
return u_clog;
end;

function wlearner(cu, clog, xiter = 1)
  wcu = zeros(Float32, size(cu,1), size(cu, 2));
  wcu[:] = 1;
  ug = wlogtunnelcw(cu, clog, wcu);
  for(iter = 1:xiter)
    wcu = weightuc(cu, clog, ug);
    ug = wlogtunnelcw(cu, clog, wcu);
    end;
  return (ug, wcu);
  end;



function xmaxind(datamat)
  kv = datamat[:, 1];
  ki = similar(kv);
  ki[:] = 1;
  for(iter = 2:size(datamat,2))
    kb = kv .< datamat[:, iter];
    ki[kb] = iter;
    kv[kb] = datamat[kb, iter];
    end;
  return ki;
  end;



function xS_getclust(data, cents, fxSd = wlogsim)
n = size(data, 1);
m = size(data, 2);

centers = size(cents,1);

dnow = ones(n, 1) * 1e25;
inow = zeros(n,1);

  for(countc = 1:centers)
    dist = -fxSd(cents[countc, :], data);

    ivec = dnow .> dist;
    xvec = inow .== countc;

    rvec = (ivec .+ xvec) .> 0;

    dnow[rvec] = dist[rvec];
    inow[rvec] = countc;

    print(".");
    end;

return (inow);
end;



function xS_getsim(data, cents, fxSd = wlogsim)
n = size(data, 1);
m = size(data, 2);
centers = size(cents,1);

sims = zeros(centers, n);
for(countc = 1:centers)
  sims[countc,:] = fxSd(cents[countc, :], data);
  end;
return (sims);
end;


function xS_normsim(sims)
  return (sims .- mean(sims, 2)) ./ std(sims, 2);
  end;


function predKit(data, cents)

  pksim = xS_getsim(data, cents, wxlogsim);
  resim = xS_normsim(pksim);
  return xmaxind(resim');
  end;



function trainErrX(cents, data, labels, fxsim = wlogsim)

  pred = xS_getclust(data, cents, fxsim) .- 1;
  
  return(sum(pred .== labels) ./ length(labels));
  
  end;



function tErrX(pred, labels)

  return(sum(pred .== labels) ./ length(labels));
  end;


function trainErrY(cents, data, loglabels, fxsim = wlogsim)

  pred = logmaker(xS_getclust(data, cents, fxsim) .- 1, [0:9]);
  #rec = sum((pred .+ loglabels) .== 2, 1) ./ sum(loglabels,1);
  #prec = sum((pred .+ loglabels) .== 2, 1) ./ sum(pred,1);
  ps = sum(pred,1);
  ts = sum(loglabels,1);
  tps = sum((pred .+ loglabels) .== 2, 1);
  return( [ps, ts, tps]);
  end;



function filtx(data, rf, rdim)

  n = size(data,1);
  rx, ry = size(rf);
  zvar = rdim .- [rx, ry];
  dx = reshape(data, n, rdim[1], rdim[2]);
  dret = zeros(n, zvar[1] + 1, zvar[2] + 1);
  for(itery = 1:ry)
    for(iterx = 1:rx)
      dret = dret .+ (dx[:, iterx:iterx+zvar[1], itery:itery+zvar[2]] .* rf[iterx, itery]);
      end;
    end;
  return reshape(dret, n, prod(zvar .+1));
  end;


using Gadfly


function multiplot(matr, geo = Geom.line, geo2 = geo)
  n = size(matr,1);
  m = size(matr,2);
  col = [1:m];
  col[:] = 1;
  k = layer(y = matr[1, :], x = [1:m], geo, color = [col], geo2);
  for(iter = 2:n)
    col[:] = iter;
    append!(k, layer(y = matr[iter, :], x = [1:m], geo, color = [col], geo2));
    end;
  return plot(k);
  end;

function pointplot(matr1, matr1b, matr2, geo = Geom.point)
  k = layer(y = matr1[:], x = matr1b'[:], geo, color = [matr2]);
  return plot(k);
  end;


function quickplot(k, dest = "def.png", path = "./")
  draw(PNG("$path$dest", 6inch, 4inch), k);
  end;

function scaff(x, z=1)
rvec = repmat([1:x[1]], 1, x[2])';
yvec = zeros(x[2],0);
xvec = zeros(x[1],0);
for(iter = 1:z)
  yvec = [yvec -rvec];
  xvec = [xvec rvec' .+ (iter-1) * x[1]];
  end;
return( yvec[:], xvec[:]);
end;


## startup ----------------------- ##

@time (data, header) = loadDat("../input/train.csv", true);

classlog = logmaker(data[:,1], [0:9]);

print(" $(sum(classlog)) ");

pixlog = data[:, 2:end] .> 128;

print(" $(size(pixlog)) $(size(classlog)) !\n");


## type A --------------

### multi feature  ---------------------------------  ##

shortening = randperm(size(pixlog,1))[1:2500];

quickplot(multiplot(sum(classlog[shortening,:], 1)),"sampclass.png");

f1 = [-1 0; 0 1];

f2 = [0 -1; 1 0];

rd = [28, 28];

fpixa = filtx(pixlog, f1, rd);

fpixb = filtx(pixlog, f2, rd);

xpixlog = [fpixa .> 0  fpixa .< 0  fpixb .> 0 fpixb .< 0];

(ugvalneg, weight1) = wlearner(xpixlog[shortening, :], classlog[shortening, :], 0);


simx = xS_getsim(xpixlog, ugvalneg');

quickplot(multiplot(simx[:, 1:100]), "simxplot.png");

quickplot(multiplot([mean(simx',1), maximum(simx',1), std(simx',1)]), "simxstat.png");


simz = xS_getsim(xpixlog, ugvalneg', wxlogsim);

quickplot(multiplot(simz[:, 1:100]), "simzplot.png");

quickplot(multiplot([mean(simz',1), maximum(simz',1), std(simz',1)]), "simzstat.png");


##
err = trainErrX(ugvalneg', xpixlog, data[:,1]);

print("\n\n train-a2 err measure X :: $(err) \n");


err2 = trainErrY(ugvalneg', xpixlog, classlog);

quickplot(multiplot(err2, Geom.line, Geom.point), "precis-rec_a2.png");


(ax, bx) = scaff([27, 27], 4);

quickplot(pointplot([ax, ax.-28] ,[bx, bx], ugvalneg[:, [3,4]]), "class34.png")




#errb = trainErrX(ugvalneg', xpixlog, data[:,1], wxlogsim);

#print("\n\n train-a2 err measure X(norm?) :: $(errb) \n");


#err2b = trainErrY(ugvalneg', xpixlog, classlog, wxlogsim);

#quickplot(multiplot(err2b, Geom.line, Geom.point), "precis-rec_simz.png");


## 

p2 = predKit(xpixlog, ugvalneg');

err3 = tErrX(p2 .-1, data[:,1]);

print("***\n Predkit error ", err3);










