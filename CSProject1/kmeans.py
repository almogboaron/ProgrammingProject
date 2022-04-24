import math
import sys

def main(arr):
    
    try:
      if len(arr) == 4:
          K = int(arr[1])
          max_iter = 200
          filename = arr[2]
          fileout = arr[3]
      elif len(arr) == 5:
          K = int(arr[1])
          max_iter = int(arr[2])
          filename = arr[3]
          fileout = arr[4]
      else:
        raise Exception("Invalid Input!")
    
    except:
      print("Invalid Input!\n")
      sys.exit(1)

    assert K > 1, "Invalid Input!"
    
    centroids, datapoints, n = init(filename, K) ## n = number of datapoints
    assert K < n, "Invalid Input!"
    clusters = [[] for _ in range(K)]
    iter_number = 0
    epsilon = 0.001
    
    while iter_number < max_iter:
        for i in range(n):
            assign(datapoints[i], clusters, centroids, K)
        cnt = 0
        for k in range(K):
            delta, new = update(centroids[k], clusters[k])
            centroids[k] = new
            if norm_calc(delta) < epsilon:
                cnt += 1
        if cnt >= K:
            break
        iter_number += 1
    out_file = write_to_file(centroids, fileout)

def init(filename, K):
    ## reads the file and return array of centroids, array of data points and the number of lines ##
    centroids = []
    datapoints = []
    n = 0
    try:
        f = open(filename, "r")
        for line in f:
            tmplst = line.rstrip().split(",")
            tmpmap = map(float,tmplst)
            datapoints.append(list(tmpmap))
            if n < K:
                tmplst = line.rstrip().split(",")
                tmpmap = map(float,tmplst)
                centroids.append(list(tmpmap))
            n += 1
        f.close()
    except:
        print("An Error Has Occurred")
        sys.exit(1)
    return centroids, datapoints, n

def assign(datapoint, clusters, centroids, K):
    try:
        tmp = []
        for j in range(K):
            delta = []
            for i in range(len(datapoint)):
                delta.append(datapoint[i] - centroids[j][i])
            tmp.append(norm_calc(delta))
        idx_cluster = tmp.index(min(tmp))
        clusters[idx_cluster].append(datapoint)
    except:
        print("An Error Has Occurred")
        sys.exit(1)

def update(centroid, cluster):
    ## updates the centroid and calculate the change ##
    new_centroid = [0 for _ in range(len(centroid))]
    delta = []
    l = len(cluster)
    try:
        for x in cluster:
            for i in range(len(x)):
                new_centroid[i] += x[i]/l
        for i in range(len(centroid)):
            delta.append(centroid[i] - new_centroid[i])
    except:
        print("An Error Has Occurred")
        sys.exit(1)
    return delta, new_centroid

def norm_calc(delta):
    s = 0
    try:
        for x in delta:
            s += pow(x, 2)
    except:
        print("An Error Has Occurred")
        sys.exit(1)
    return math.sqrt(s)

def write_to_file(centroids,fileout):
    f = open(fileout,"w")
    try:
        for m in centroids:
            for i in range(len(m)):
                m[i] = round(m[i],4)
            f.write(str(m).strip("[]") + "\n")
        f.close()
    except:
        print("An Error Has Occurred")
        sys.exit(1)
    return f


if __name__ == "__main__":
    main(sys.argv)




