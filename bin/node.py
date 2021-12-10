class Node:
    def __init__(self, name,peak):
        """
        Args:
            name:
        """
        self.name = name
        self.children = []
        self.peak = peak
        self.data = [[0] for i in range(0, 288)]

    def addChild(self, name,peak):
        """
        Args:
            name:
        """
        node = Node(name, peak)
        self.children.append(node)
        return node

    def addData(self, dataList):
        """
        Args:
            dataList:
        """
        for i in range(0, 288):
            self.data[i] += dataList[i]
