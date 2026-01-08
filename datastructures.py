from dataclasses import dataclass

@dataclass
class Hairpin:
    """
    Describes a hairpin DNA structure

    stem1: A string representing the DNA bases from the 5' end of a hairpin up to the start of the loop (exclusive)
    loop: A string representing the DNA bases of the loop of a hairpin
    stem2: A string representing the DNA bases from the end of the loop (exclusive) to the 3' end of a hairpin
    """
    stem1: str
    loop: str
    stem2: str

@dataclass
class Primer:
    '''
    Describes a primer

    seq: A string representing the DNA bases of the primer
    start: A integer representing the start location of the primer on the original inputted sequence
    end: A integer representing the end location of the primer on the original inputted sequence
    meltTemp: A float representing the melting temperature of the primer
    gcContent: A float representing the gcContent of the primer
    gcInClamp: A integer representing the number of G's and C's in the GC clamp of the primer
    '''
    seq: str
    start: int
    end: int
    meltTemp: float
    gcContent: float
    gcInClamp: int

class Node:
    '''
    A generic trie node
    '''
    def __init__(self):
        self.position = None
        self.nextNodes = {}
    
    def getPosition(self):
        return self.position
    
    def getNextNode(self, key):
        return self.nextNodes.get(key)
    
    def getNextNodes(self):
        return self.nextNodes
    
    def setPosition(self, position):
        self.position = position
    
    def addNextNode(self, key, node):
        self.nextNodes[key] = node

    
class SuffixTrie:
    '''
    A generic suffix trie datastructure
    '''
    def __init__(self):
        self.root = Node()
        self.edges = 0

    def addWord(self, word, position):
        curNode = self.root
        # Loop through all the letters of the word
        for j in range(len(word)):
            # Check if the letter is already a key for a next node, if not, create a next node with the letter as the key.
            if (curNode.getNextNode(word[j]) == None):
                curNode.addNextNode(word[j], Node())
                self.edges += 1
            curNode = curNode.getNextNode(word[j])
        # Set the position of the last node indicating where the word can be found
        curNode.setPosition(position)
    
    def addSuffixes(self, word):
        #Add each suffix to the trie
        for i in range(len(word)):
            suffix = word[i:]
            self.addWord(suffix, i)
    
    def getEdges(self):
        return self.edges
    
    # Debugger method to print out paths of the trie
    def printPathsRoot(self):
        self.printPaths(self.root, "")
    
    # Helper method for printPathsRoot
    def printPaths(self, node, word):
        # Base case
        if len(node.getNextNodes()) == 0:
            print(word)
        for key in node.getNextNodes().keys():
            self.printPaths(node.getNextNode(key), word + key)


    
