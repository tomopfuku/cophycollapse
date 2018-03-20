package cophymaru

import (
	"fmt"
	"io/ioutil"
	"math"
	"math/rand"
	"os"
	"strings"
	"time"
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
