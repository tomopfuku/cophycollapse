package cophycollapse

import (
	"fmt"
	"os"
)

//AssertNumChars will check to make sure that two character vectors are the same lengths (as they always should be)
func AssertNumChars(chars1, chars2 []float64) {
	if len(chars1) != len(chars2) {
		fmt.Println("ERROR. Something went wrong with the mapping of traits to the tips.")
		os.Exit(0)
	}
}

//AssertNumMis will make sure that the boolian slices are the same length
func AssertNumMis(chars1, chars2 []bool, node *Node) {
	if len(chars1) != len(chars2) {
		fmt.Println(len(chars1), len(chars2))
		fmt.Println(node.Newick(true))
		fmt.Println("ERROR. Something went wrong with the mapping of traits to the tips.")
		os.Exit(0)
	}
}
