package cophycollapse

import (
	"fmt"
	"io/ioutil"
	"math"
	"math/rand"
	"os"
	"strings"
	"time"

	"gonum.org/v1/gonum/mat"
)

//InternalNodeSlice will return a slice containing only internal nodes
func InternalNodeSlice(nodes []*Node) (inNodes []*Node) {
	for _, n := range nodes {
		if len(n.CHLD) == 2 {
			inNodes = append(inNodes, n)
		}
	}
	return
}

//InitParallelPRNLEN will set up empty slices for the prnlens
func InitParallelPRNLEN(nodes []*Node) {
	for _, n := range nodes {
		n.CONPRNLEN = make([]float64, len(nodes[0].CONTRT))
	}
}

//TreeLength will return the total length of a slice of nodes
func TreeLength(nodes []*Node) float64 {
	len := 0.
	for _, n := range nodes[1:] {
		len += n.LEN
	}
	return len
}

//MakeRandomStartingBranchLengths will initialize a tree with a set of random branch lengths
func MakeRandomStartingBranchLengths(tree *Node) {
	nodes := tree.PreorderArray()
	for _, n := range nodes {
		s1 := rand.NewSource(time.Now().UnixNano())
		r1 := rand.New(s1)
		u := r1.Float64()
		n.LEN = u
	}
}

//ReadLine is like the Python readline() and readlines()
func ReadLine(path string) (ln []string) {
	b, err := ioutil.ReadFile(path)
	if err != nil {
		fmt.Println(err)
		fmt.Println("There was an error when reading in the file:", path, ". Are you sure that it exists?")
		os.Exit(0)
	}
	ss := string(b)
	ln = strings.Split(ss, "\n")
	return
}

//Rexp will draw a random exponential number
func Rexp(lambda float64) (e float64) {
	u := rand.Float64()
	e = math.Log(1-u) / (-lambda)
	e = e / 2.
	return
}

//ReadFossils will read in a list of fossil tips one line at a time into a slice
//TODO: get this working
func ReadFossils(path string) (fos []string) {
	l := ReadLine(path)[0]
	fos = strings.Split(l, ",")
	return
}

//Max returns the maximum value in a map of ints used like a set
func Max(l map[int][]int) (biggest int) {
	biggest = -10000000
	for i := range l {
		if i > biggest {
			biggest = i
		}
	}
	return
}

//MaxClustLab returns the maximum value in a map of ints used like a set
func MaxClustLab(l map[int]float64) (biggest int) {
	biggest = -10000000
	for i := range l {
		if i > biggest {
			biggest = i
		}
	}
	return
}

func matPrint(X mat.Matrix) {
	fa := mat.Formatted(X, mat.Prefix(""), mat.Squeeze())
	fmt.Printf("%v\n", fa)
}

//SetIdentityMatrix will return an identiy matrix with dimensions ntax,ntax
func SetIdentityMatrix(dim int) *mat.Dense {
	matrix := mat.NewDense(dim, dim, nil)
	for i := 0; i < dim; i++ {
		for j := 0; j < dim; j++ {
			if i == j {
				matrix.Set(i, j, 1.0)
			} else {
				matrix.Set(i, j, 0.0)
			}
		}
	}
	return matrix
}

//ColumnMatrixToSlice will convert a Dx1 matrix to a slice
func ColumnMatrixToSlice(m *mat.Dense) []float64 {
	r, _ := m.Dims()
	var newSlice []float64
	for i := 0; i < r; i++ {
		cur := m.At(i, 0)
		newSlice = append(newSlice, cur)
	}
	return newSlice
}

//SymDenseConvert will convert a matrix of type *Dense to *SymDense
func SymDenseConvert(m *mat.Dense) *mat.SymDense {
	r, c := m.Dims()
	if r != c {
		fmt.Println("Matrix is not symmetric.")
		os.Exit(0)
	}
	new := mat.NewSymDense(r, nil)
	for i := 0; i < r; i++ {
		for j := 0; j < r; j++ {
			new.SetSym(i, j, m.At(i, j))
		}
	}
	return new
}

//LogGammaFn will return the log of the gamma function of x. X should be provided as a float64, but be a whole number
func LogGammaFn(x float64) float64 {
	sum := 0
	max := int(x)
	for i := 1; i < max; i++ {
		sum += i
	}
	return float64(sum)
}

func lengthToRoot(n *Node) float64 {
	if n.PAR == nil {
		return 0.
	}
	cur := n
	rootPath := 0.
	for {
		if cur.PAR == nil {
			//fmt.Println(cur.Newick(true))
			break
		}
		rootPath += cur.LEN
		//fmt.Println(cur.LEN, rootPath)
		cur = cur.PAR
	}
	return rootPath
}
