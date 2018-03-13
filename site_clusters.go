package cophymaru

import (
	"math/rand"
	"time"
)

//InitializeClusters will initialize branch length clusters for the cluster model
func InitializeClusters(tree *Node) (c []int) {
	ntraits := len(tree.CONTRT)
	K := ntraits / 100
	if K == 0 {
		K = 5
	}
	for range tree.CONTRT {
		s1 := rand.NewSource(time.Now().UnixNano())
		r1 := rand.New(s1)
		clusterAssignment := r1.Intn(K)
		c = append(c, clusterAssignment)
	}
	nodes := tree.PostorderArray()
	for _, n := range nodes {
		for i := 0; i < K; i++ {
			n.ClustLEN[i] = n.LEN //append(n.ClustLEN[i], n.LEN)
		}
	}
	return
}
