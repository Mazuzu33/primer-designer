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

def calcMaxVocab(alphaSize, wordLen):
    '''
    Helper function for calculating linguistic complexity

    alphaSize: The amount of letters in the alphabet (int)
    wordLen: The length of the word (int)
    Returns: The maximum vocabulary of words of length 1 to m that can be formed by taking substrings of a word of wordLen (int)
    ''' 
    maxVocab = 0
    # Calculate the maximum vocabulary of words of length i
    for i in range(1, wordLen + 1):
        # Represents the max number of unique words of length i with regards to the alphabet. 
        val1 = alphaSize ** i
        # Represents the max number of unique words of length i that can be obtained by taking substrings of a word of wordLen. 
        # Invalid when wordLen > alphaSize
        val2 = wordLen - i + 1
        # Choose the minimum of the two values
        maxVocab += min(val1, val2)
    return maxVocab

def calcLinComp(seq):
    '''
    Calculates the linguistic complexity (LC) of a sequence where LC is defined as the ratio of the number of substrings of any
    length in the given sequence to the maximum possible number of substrings obtainable from a sequence of the given sequence's length. 
    
    seq: A DNA sequence (string)
    '''
    trie = SuffixTrie()
    trie.addSuffixes(seq)
    uniqueSubWords = trie.getEdges()
    maxVocab = calcMaxVocab(4, len(seq))
    return uniqueSubWords / maxVocab