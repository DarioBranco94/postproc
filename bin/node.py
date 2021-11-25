class Node:
    def __init__(self, name):
        """
        Args:
            name:
        """
        self.name = name
        self.children = []
        self.data = [[0] for i in range(0, 288)]

    def addChild(self, name):
        """
        Args:
            name:
        """
        node = Node(name)
        self.children.append(node)
        return node

    def addData(self, dataList):
        """
        Args:
            dataList:
        """
        for i in range(0, 288):
            self.data[i] += dataList[i]
