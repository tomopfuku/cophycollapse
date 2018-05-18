package cophymaru

import (
	"math"
)

//DistanceMatrix will store a continuous character matrix converted into a pairwise distance matrix for each site
type DistanceMatrix struct {
	MatSites map[int][]float64 // pairwise comparison is the key, value is a slice containing the comparison between taxa across all sites
	//MatSitesT map[int]mat.Matrix // same as above, but the site trait vectors are transposed
	//MatSitesProd map[int]*mat.Dense // stores A*A^T for each site A
	PWDistLabs []int // stores the order that traits are stored in for the untransposed trait vectors
	//Dim       int
	NSites    float64
	IntNSites int
}

/*
//DM will return the variance-covariance matrix calculated from a tree with branch lengths
func DM(nodes []*Node) *DistanceMatrix { //map[int][]float64 {
	var tips []*Node
	for _, n := range nodes {
		if len(n.CHLD) == 0 {
			tips = append(tips, n)
		}
	}
	//dim := len(tips)
	sites := make(map[int][]float64)
	nsites := len(tips[0].CONTRT)
	paircount := 0
	var pwlabs []int
	for i, n1 := range tips {
		for j, n2 := range tips {
			if i == j || i > j {
				continue
			} else {
				var dls []float64
				for t := 0; t < nsites; t++ {
					dist := math.Sqrt(math.Pow(n1.CONTRT[t]-n2.CONTRT[t], 2))
					dls = append(dls, dist)
					fmt.Println(n1.NAME, n2.NAME)
				}
				sites[paircount] = dls
				pwlabs = append(pwlabs, paircount)
				paircount++
			}
		}
	}
	dm := new(DistanceMatrix)
	dm.MatSites = sites
	//dm.PWDistLabs = pwlabs
	return dm
}
*/

//DM will return the variance-covariance matrix calculated from a tree with branch lengths
func DM(nodes []*Node) *DistanceMatrix { //map[int][]float64 {
	var tips []*Node
	for _, n := range nodes {
		if len(n.CHLD) == 0 {
			tips = append(tips, n)
		}
	}
	//dim := len(tips)
	sites := make(map[int][]float64)
	nsites := len(tips[0].CONTRT)
	paircount := 0
	var pwlabs []int
	for t := 0; t < nsites; t++ {
		var dls []float64
		for i, n1 := range tips {
			for j, n2 := range tips {
				if i == j || i > j {
					continue
				} else {
					dist := math.Sqrt(math.Pow(n1.CONTRT[t]-n2.CONTRT[t], 2))
					dls = append(dls, dist)
					//fmt.Println(n1.NAME, n2.NAME)
				}
			}
		}
		sites[t] = dls
		pwlabs = append(pwlabs, paircount)
		paircount++
	}
	dm := new(DistanceMatrix)
	dm.MatSites = sites
	dm.PWDistLabs = pwlabs
	dm.NSites = float64(nsites)
	dm.IntNSites = nsites
	return dm
}
