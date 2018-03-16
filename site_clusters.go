package cophymaru

import (
	"fmt"
	"math/rand"
)

//InitializeClusters will initialize branch length clusters for the cluster model
func InitializeClusters(chain *MCMC) {
	ntraits := len(chain.TREE.CONTRT)
	K := ntraits / 100
	if K < 1 {
		K = 5
	}
	var uniqueK []int
	for k := 0; k < K; k++ {
		uniqueK = append(uniqueK, k)
		//break
	}
	chain.UNIQUEK = uniqueK
	for i := 0; i < int(chain.NSITES); i++ {
		clusterAssignment := uniqueK[rand.Intn(len(uniqueK))]
		if i == 0 {
			fmt.Println("site", i, clusterAssignment)
		}
		chain.CLUS = append(chain.CLUS, clusterAssignment)
		if _, ok := chain.CLUSTERSET[clusterAssignment]; ok {
			chain.CLUSTERSET[clusterAssignment] = append(chain.CLUSTERSET[clusterAssignment], i)
		} else {
			var curfill []int
			curfill = append(curfill, i)
			chain.CLUSTERSET[clusterAssignment] = curfill
		}
	}
	for _, n := range chain.NODES {
		n.ClustLEN = make(map[int]float64)
		for i := 0; i < K; i++ {
			n.ClustLEN[i] = n.LEN //append(n.ClustLEN[i], n.LEN)
		}
	}
	return
}

/*/SiteDistMatrix will calculate the distance matrix at each site
func SiteDistMatrix(nodes []*Node) {
	for i, ni := range nodes {
		if len(ni.CHLD) != 0 {
			continue
		}
		for j, nj := range nodes {
			if len(ni.CHLD) != 0 {
				continue
			} else if i == j{

			}
		}
	}
}
*/

//StartingSiteLen makes starting branch lengths for each site for clustering.
func StartingSiteLen(chain *MCMC) {
	for i := 0; i < int(chain.NSITES); i++ {
		chain.CLUS = append(chain.CLUS, i)
		if _, ok := chain.CLUSTERSET[i]; ok {
			chain.CLUSTERSET[i] = append(chain.CLUSTERSET[i], i)
		} else {
			var curfill []int
			curfill = append(curfill, i)
			chain.CLUSTERSET[i] = curfill
		}
	}
	chain.UNIQUEK = chain.CLUS
	for _, n := range chain.NODES {
		n.ClustLEN = make(map[int]float64)
		for i := 0; i < len(chain.CLUS); i++ {
			n.ClustLEN[i] = n.LEN //append(n.ClustLEN[i], n.LEN)
		}
	}
}

//AssignClustLens will assign the lengths associated with a particular cluster to the branch lengths
func AssignClustLens(chain *MCMC, cluster int) {
	for _, n := range chain.NODES {
		n.LEN = n.ClustLEN[cluster]
	}
}
