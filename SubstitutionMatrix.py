class SubstitutionMatrix:
    """Score for substituting one letter for another.  Indexed by a pair of letters; e.g., blosum62['A','R'] => -1"""
    def __init__(self, alphabet, scores):
        """Create a new instance from a list of lists, scores, where the indices follow the order of letters in the alphabet."""
        self.matrix = {}
        for i in range(len(alphabet)):
            for j in range(len(alphabet)):
                self.matrix[alphabet[i],alphabet[j]] = scores[i][j]

    def __getitem__(self, pair):
        return self.matrix[pair]
