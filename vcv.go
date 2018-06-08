package cophycollapse

import (
	"gonum.org/v1/gonum/mat"
)

//PhyloVCV will return the variance-covariance matrix calculated from a tree with branch lengths
func PhyloVCV(nodes []*Node) *mat.Dense {
	var tips []*Node
	for _, n := range nodes {
		if len(n.CHLD) == 0 {
			tips = append(tips, n)
		}
	}
	dim := len(tips)
	vcv := mat.NewDense(dim, dim, nil)
	for i := range tips {
		for j := range tips {
			if i == j {
				variance := lengthToRoot(tips[i])
				vcv.Set(i, j, variance)
			} else {
				mrca := MRCA(tips[i], tips[j])
				variance := lengthToRoot(mrca)
				vcv.Set(i, j, variance)
			}
		}
	}
	return vcv
}

//MRCA will return the MRCA of two tips
func MRCA(n1, n2 *Node) *Node {
	traceback := make(map[*Node]int)
	cur := n1
	for {
		traceback[cur] = 0
		if cur.PAR == nil {
			break
		}
		cur = cur.PAR
	}
	cur = n2
	for {
		if _, ok := traceback[cur]; ok {
			break
		}
		cur = cur.PAR
	}
	return cur
}
