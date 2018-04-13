package cophymaru

import (
	"gonum.org/v1/gonum/mat"
)

//CharAlignment will store all characters for a collapsed Gibbs clustering analysis, along with their transposed vectors
type CharAlignment struct {
	MatSites     map[int]*mat.Dense // site index is the key, traits for each taxon are the values
	MatSitesT    map[int]mat.Matrix // same as above, but the site trait vectors are transposed
	MatSitesProd map[int]*mat.Dense // stores A*A^T for each site A
	SiteOrder    []string           // stores the order that traits are stored in for the untransposed trait vectors
	Dim          int
	NSites       int
}

//InitTraitMatrices will initialize the CharAlignment struct
func InitTraitMatrices(traits map[string][]float64) *CharAlignment {
	sites := make(map[int]*mat.Dense)
	tsites := make(map[int]mat.Matrix)
	prod := make(map[int]*mat.Dense)
	var taxOrder []string
	ntax := len(traits)
	tnum := 0
	for k := range traits {
		taxOrder = append(taxOrder, k)
		for i, tr := range traits[k] {
			if _, ok := sites[i]; ok {
				sites[i].Set(tnum, 0, tr)
			} else {
				sites[i] = mat.NewDense(ntax, 1, nil)
				sites[i].Set(tnum, 0, tr)
			}

		}
		tnum++
	}
	ALN := new(CharAlignment)
	ALN.MatSites = sites
	ALN.Dim = len(taxOrder)
	ALN.SiteOrder = taxOrder
	count := 0
	for k := range sites {
		cur := sites[k]
		tsites[k] = cur.T()
		p := mat.NewDense(ntax, ntax, nil)
		p.Product(cur, tsites[k])
		prod[k] = p
		count++
	}
	ALN.NSites = count
	ALN.MatSitesT = tsites
	ALN.MatSitesProd = prod
	return ALN
}
