import json
import os
import sys

import numpy as np
import pandas as pd

from util.log_handler import logger

tree_name_separator = ">"

treefile_template = "{base}_tree.json"
metafile_template = "{base}_meta.json"
datafile_template = "{base}_data.json"


class TimingNode(dict):
    def __init__(self, name, ancestry=[]):
        super().__init__()
        self.__dict__ = self
        self.name = name
        self.ancestry = ancestry
        self.children = []

    def __eq__(self, other):
        return self.get_ancestry_name() == other.get_ancestry_name()

    def __hash__(self):
        return hash(self.get_ancestry_name())

    @staticmethod
    def from_dict(dict_):
        node = TimingNode(dict_["name"], dict_["ancestry"])
        node.children = list(map(TimingNode.from_dict, dict_["children"]))
        return node

    @staticmethod
    def name_from_ancestry_name(ancestry_name):
        return ancestry_name.split(tree_name_separator)[-1]

    def to_name_list(self):
        out = [c.get_name() for c in self.children]
        for c in self.children:
            out += c.to_name_list()
        return out

    def to_ancestry_name_list(self):
        out = [c.get_ancestry_name() for c in self.children]
        for c in self.children:
            out += c.to_ancestry_name_list()
        return out

    def to_list(self):
        out = [c for c in self.children]
        for c in self.children:
            out += c.to_list()
        return out

    def add_child(self, child):
        self.children.append(child)

    def remove_child(self, child):
        self.children.remove(child)

    def search_children(self, name):
        out = None
        if self.get_name() == name:
            out = self
        else:
            for node in self.children:
                out = node.search_children(name)
                if out:
                    break
        return out

    def set_children(self, children):
        self.children = children

    def get_ancestry_name(self):
        tmp = self.ancestry + [self.get_name()]
        return tree_name_separator.join(tmp)

    def get_name(self):
        return self.name

    def get_ancestry(self):
        return self.ancestry

    def intersection(self, other):
        left = set(self.to_list())
        right = set(other.to_list())

        intersection = left & right

        return list(intersection)

    def sub(self, other):
        left = set(self.to_list())
        right = set(other.to_list())

        sub = left - right

        all_children = set()
        for n in sub:
            all_children = all_children | set(n.to_list())

        independent_nodes = sub - all_children

        return list(independent_nodes)


class TimingTree:
    def __init__(self, root, data, meta_data):
        self.root = root
        self.data = data
        self.meta_data = meta_data

    @classmethod
    def from_logfile(cls, logfile, logfile_reader):
        logger.info("initializing TimingTree from file {}".format(logfile))
        timing_data, meta_data = logfile_reader(logfile)

        dfs = []
        roots = []
        for i, e in enumerate(timing_data):
            root, df = TimingTree.data_to_tree(e, meta_data, i)
            dfs.append(df)
            roots.append(root)

        return cls(roots, dfs, meta_data)

    @classmethod
    def from_json(cls, filename):
        with open(metafile_template.format(base=filename)) as f:
            meta_data = json.load(f)

        roots = []
        dfs = []
        for i in range(meta_data["n_tables"]):
            filename_i = "{}_{}".format(filename, i)
            with open(treefile_template.format(base=filename_i)) as f:
                json_dict = json.load(f)
            roots.append(TimingNode.from_dict(json_dict))

            dfs.append(
                pd.read_json(datafile_template.format(base=filename_i), orient="table")
            )

        return cls(roots, dfs, meta_data)

    @staticmethod
    def data_to_tree(timing_data, meta_data, i_table):
        n = meta_data["entries"][i_table]

        root = TimingNode("root")

        timer_names = timing_data["name"]
        date = [meta_data["finish_time"]]

        row_index = pd.MultiIndex.from_product([timer_names, date])
        column_index = [
            key for key in timing_data.keys() if key not in ["indent", "name"]
        ]

        # traverse the tree
        matrix = np.zeros((n, len(column_index)))
        for i in range(n):
            last = root
            ancestry = [last.get_name()]
            for _ in range(timing_data["indent"][i]):
                last = last.children[-1]
                ancestry.append(last.get_name())

            node = TimingNode(timing_data["name"][i], ancestry=ancestry)
            last.add_child(node)

            matrix[i, :] = [
                val[i]
                for key, val in timing_data.items()
                if key not in ["indent", "name"]
            ]

        return root, pd.DataFrame(matrix, index=row_index, columns=column_index)

    @staticmethod
    def input_exists(filename):
        treefile = treefile_template.format(base="{}_0".format(filename))
        metafile = metafile_template.format(base="{}".format(filename))
        datafile = datafile_template.format(base="{}_0".format(filename))

        tree_exists = os.path.exists(treefile)
        meta_exists = os.path.exists(metafile)
        data_exists = os.path.exists(datafile)

        if not tree_exists:
            logger.info("Treefile does not exist: {}".format(treefile))
        if not meta_exists:
            logger.info("Metafile does not exist: {}".format(metafile))
        if not data_exists:
            logger.info("Datafile does not exist: {}".format(datafile))

        return tree_exists and meta_exists and data_exists

    def json_dump(self, filename):
        for i in range(self.meta_data["n_tables"]):
            filename_i = "{}_{}".format(filename, i)
            with open(treefile_template.format(base=filename_i), "w") as f:
                json.dump(self.root[i], f, indent=2)

            self.data[i].to_json(
                datafile_template.format(base=filename_i), orient="table", indent=2
            )

        with open(metafile_template.format(base=filename), "w") as f:
            json.dump(self.meta_data, f, indent=2)

    def find(self, name, i_table):
        out = self.root[i_table].search_children(name)
        if not out:
            logger.info("did not find node {}".format(name))
        return out

    def find_ancestor(self, node, i_table, k=1):
        ancestor_name = node.ancestry[-k]
        return self.find(ancestor_name, i_table)

    def add_data(self, other):
        # add data
        if len(self.data) != len(other.data):
            logger.critical("cannot add reports with unequal number of tables")
            sys.exit(1)

        for i_table in range(len(self.data)):
            self.data[i_table], _ = self.data[i_table].align(other.data[i_table])
            self.data[i_table].loc[(slice(None), other.meta_data["finish_time"]), :] = (
                other.data[i_table]
            )

    def add_meta_data(self, other):
        # a single timing file will produce single values in the meta data
        # when adding multiple timing files together, these must be converted to lists
        convert = not isinstance(self.meta_data["finish_time"], list)
        # add meta data
        for key, val in self.meta_data.items():
            if key not in other.meta_data.keys():
                logger.critical(
                    "incompatible meta_data keys:\nfound: {}\nexpected: {}".format(
                        self.meta_data.keys(), other.meta_data.keys()
                    )
                )
                sys.exit(1)
            other_val = other.meta_data[key]
            # The only exception is n_tables which MUST be equal across timing files,
            # so it's not converted to a list
            if key == "n_tables":
                if val != other_val:
                    logger.critical(
                        "Trying to add a timing file with different number of tables\n"
                        + "expected {}, got {}".format(val, other_val)
                    )
                    sys.exit(1)
            else:
                if convert:
                    self.meta_data[key] = [val]
                self.meta_data[key].append(other_val)

    def add_tree(self, other):
        # manage tree
        new_nodes = other.root.sub(self.root)
        new_names = [n.get_name() for n in new_nodes]
        logger.info("growing the tree by new nodes: {}".format(" ".join(new_names)))
        if len(new_nodes) > 0:
            self.grow(new_nodes)

    def add(self, other):
        self.add_data(other)
        self.add_meta_data(other)
        self.root = other.root

    def grow(self, nodes, i_table):
        levels = [(len(nodes[i].get_ancestry()), i) for i in range(len(nodes))]
        levels.sort()
        _, permutation = zip(*levels)
        nodes = [nodes[i] for i in permutation]
        for n in nodes:
            # look for node in old tree (self.root)
            old_node = self.find(n.get_name(), i_table)
            if old_node:
                logger.info(
                    "trying to grow tree by node which is already present: "
                    + "{}. Remove.".format(old_node.get_ancestry_name())
                )
                old_parent = self.find_ancestor(old_node, i_table, k=1)
                old_parent.remove_child(old_node)
            logger.info(
                "adding new node {} to {}".format(
                    n.get_ancestry_name(), n.get_ancestry()[-1]
                )
            )
            new_parent = self.find_ancestor(n, i_table, k=1)
            new_parent.add_child(n)
        return
