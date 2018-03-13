package cophymaru

import (
	"math/rand"
	"time"
)

//InitializeClusters will initialize branch length clusters for the cluster model
func InitializeClusters(chain *MCMC) {
	ntraits := len(chain.TREE.CONTRT)
	K := ntraits / 100
	if K == 0 {
		K = 5
	}
	for i := range chain.TREE.CONTRT {
		s1 := rand.NewSource(time.Now().UnixNano())
		r1 := rand.New(s1)
		clusterAssignment := r1.Intn(K)
		chain.CLUS = append(chain.CLUS, clusterAssignment)
		var uniqueK []int
		if _, ok := chain.CLUSTERSET[clusterAssignment]; ok {
			chain.CLUSTERSET[clusterAssignment] = append(chain.CLUSTERSET[clusterAssignment], i)
		} else {
			var curfill []int
			curfill = append(curfill, i)
			chain.CLUSTERSET[clusterAssignment] = curfill
			uniqueK = append(uniqueK, clusterAssignment)
		}
		chain.UNIQUEK = uniqueK
	}
	for _, n := range chain.NODES {
		for i := 0; i < K; i++ {
			n.ClustLEN[i] = n.LEN //append(n.ClustLEN[i], n.LEN)
		}
	}
	return
}
