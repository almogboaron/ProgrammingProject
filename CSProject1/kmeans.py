import math

def main(K, filename, max_iter=200):
    centroids, datapoints, n = init(filename, K) ## n = number of datapoints
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
        if cnt == K:
            break
        iter_number += 1
    out_file = write_to_file(centroids)
    return out_file

def init(filename, K):
    ## reads the file and return array of centroids, array of data points and the number of lines ##
    centroids = []
    datapoints = []
    n = 0
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

    return centroids, datapoints, n

def assign(datapoint, clusters, centroids, K):
    tmp = []
    for j in range(K):
        delta = []
        for i in range(len(datapoint)):
            delta.append(datapoint[i] - centroids[j][i])
        tmp.append(norm_calc(delta))
    idx_cluster = tmp.index(min(tmp))
    clusters[idx_cluster].append(datapoint)

def update(centroid, cluster):
    ## updates the centroid and calculate the change ##
    new_centroid = [0 for _ in range(len(centroid))]
    delta = []
    l = len(cluster)
    for x in cluster:
        for i in range(len(x)):
            new_centroid[i] += x[i]/l
    for i in range(len(centroid)):
        delta.append(centroid[i] - new_centroid[i])
    return delta, new_centroid

def norm_calc(delta):
    s = 0
    for x in delta:
        s += pow(x, 2)
    return math.sqrt(s)

def write_to_file(centroids):
    f = open("output.txt","w")
    for m in centroids:
        f.write(str(m).strip("[]") + "\n")
    f.close()
    return f

if __name__ == "__main__":
    print(main(3,"input_1.txt",100))




