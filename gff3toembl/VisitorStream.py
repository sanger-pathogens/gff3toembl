import sys
from gt import CustomStream

class VisitorStream(CustomStream):

    def __init__(self, instream, visitor):
        CustomStream.__init__(self)
        self.instream = instream
        self.visitor = visitor

    def next(self):
        node = self.instream.next_tree()
        if node:
            node.accept(self.visitor)
        return node
