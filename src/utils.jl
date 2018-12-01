## teardown

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
