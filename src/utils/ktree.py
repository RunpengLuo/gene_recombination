from statistics import median

class KNode:
    def __init__(self, val="", ref_id="", rind=0) -> None:
        # self.prev = []
        # self.next = []
        self.val = val
        self.dist = [ref_id]
        self.pos = [rind]

    # def init_root(self):
    #     self.prev = None

class KTree:
    def __init__(self, k: int) -> None:
        self.k = k
        self.root = KNode()
        # self.root.init_root()

        # self.prefix_dict = {}
        self.key_dict = {}
        # self.suffix_dict = {}

        self.num_node = 0
        self.num_insert = 0
        self.num_merge = 0

        return
    def insert_by_val(self, val: str, ref_id: str, rind: int):
        self.num_insert += 1
        if val in self.key_dict:
            # collision, potential sharing sequence
            exist_node = self.key_dict[val]
            exist_node.dist.append(ref_id)
            exist_node.pos.append(rind)
            self.num_merge += 1
        else:
            self.num_node += 1

            node = KNode(val, ref_id, rind)
            prefix = val[:self.k - 1]
            suffix = val[1:]

            # connect all the existing suffix overlapping nodes
            # for next_node in self.prefix_dict.get(suffix, []):
            #     next_node.prev.append(node)

            # connect all the existing prefix overlapping nodes
            # for prev_node in self.suffix_dict.get(prefix, []):
            #     prev_node.next.append(node)

            # if prefix not in self.prefix_dict:
            #     self.prefix_dict[prefix] = []
            # self.prefix_dict[prefix].append(node)
            # if suffix not in self.suffix_dict:
            #     self.suffix_dict[suffix] = []
            # self.suffix_dict[suffix].append(node)

            self.key_dict[val] = node

        return
    
    def query_by_val(self):
        return

    def __str__(self) -> str:
        return """
            --kTree--\n
            k: {0}\n
            Num of node: {1}\n
            Num of insertion: {2}\n
            Num of merging: {3}\n""".format(self.k, self.num_node, self.num_insert, self.num_merge)