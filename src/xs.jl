

## cluster imp.

function xS_cluster(data, centers, iter, fxSd = xL5_dist)
n = size(data, 1);
m = size(data, 2);
x = randperm(n);

cents, dnow = nlnmodeB(data, centers, mean(data, dims=1));

println(size(dnow), size(data), size(cents));

##dnow = ones(n, 1) * 1e25;

dnow = dnow .* 12.0; 

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

##println(size(dnow), "?");

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


function xL6_dist(cent, data)
##  print(size(cent), size(data));

  return sqrt.(sum((data .- cent').^2, dims=2));
  end;


