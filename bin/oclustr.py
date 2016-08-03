#!/usr/bin/env python

import numpy as np
import networkx as nx
import argparse


def approx_intra_sim(v, g):
    card_v = g.degree(v) + 1.
    if card_v <= 1:
        raise RuntimeError('Cardinality of V* is less than two for node {0}'.format(v))
    return np.sum([dat['weight'] for u1, u2, dat in g.edges(v, data=True) if u2 is not v]) / (card_v - 1.)


def relative_density(v, g):
    v_star = g[v]
    deg_v = g.degree(v)
    density_v = len([u for u in v_star if g.degree(u) <= deg_v])
    return density_v / float(deg_v)


def relative_compactness(v, g):
    v_star = g[v]
    deg_v = g.degree(v)
    sim_v = approx_intra_sim(v, g)
    compact_v = len([u for u in v_star if u is not v and sim_v >= approx_intra_sim(u, g)])
    return compact_v / float(deg_v)


def relevance(v, g):
    return 0.5 * (relative_density(v, g) + relative_compactness(v, g))


class Cover:
    def __init__(self, g):
        self.g = g
        self.hubs = nx.DiGraph()

    def update(self, vlist):
        for v in vlist:
            self.add(v)

    def add(self, v):
        """
        For the given vertex, build its star graph and add this
        to the list of covering stars for graph.
        :param v: the hub vertex
        """

        # create new node or update existing one
        if v not in self.hubs:
            self.hubs.add_node(v, attr_dict={'analyzed': False, 'hub': True, 'owned_by': set()})
        self.hubs.node[v]['owned_by'].add(v)
        self.hubs.node[v]['hub'] = True

        v_adj = self.g[v].keys()
        for vi in v_adj:
            # for each satellite, create new node or update existing.
            # we might be covering the same node with more than one star
            if vi not in self.hubs:
                self.hubs.add_node(vi, attr_dict={'analyzed': False, 'hub': False, 'owned_by': set()})
            self.hubs.node[vi]['owned_by'].add(v)
            self.hubs.add_edge(v, vi)

    def shared(self, v):
        """
        Return the list of adjacent vertices (satellites) which are also shared with other hubs
        in cover.
        :param v: the hub vertex
        :return: the list of shared satellites
        """
        v_adj = self.hubs[v]
        # vertices with more than one owner
        return set([ui for ui in v_adj if self.hubs.degree(ui) > 1])

    def unshared(self, v):
        """
        Return the list of adjacent vertices (satellites) which are not covered by any other
        hub in cover.
        :param v: the hub vertex
        :return: the list of unshared satellites
        """
        v_adj = self.hubs[v]
        # vertices with only 1 owner
        return set([ui for ui in v_adj if self.hubs.degree(ui) <= 1])

    def is_useful(self, v):
        """
        A useful hub is one which has more unshared than shared satellites in the
        cover of G.
        :param v: the hub vertex to test.
        :return: True if this hub is a useful star.
        """
        ns = len(self.shared(v))
        nu = len(self.unshared(v))
        return ns <= nu

    # def highest_containing(self, u):
    #     """
    #     For a given vertex u, find which owning hub has the highest degree. If there is a
    #     tie between owning hubs, pick one at random.
    #     :param u: the vertex to test.
    #     :return: the hub vertex with high degree
    #     """
    #     owners = self.hubs.node[u]['owned_by']
    #
    #     # owners degrees, not including self
    #     owner_deg = dict((v, self.g.degree(v)) for v in owners if v != u)
    #
    #     if len(owner_deg) == 0:
    #         raise RuntimeWarning('{0} is an isolated vertex containing only itself'.format(u))
    #     elif len(owner_deg) == 1:
    #         # only one non-self owner
    #         return owner_deg.popitem()
    #     else:
    #         # sort owners by descending degree
    #         deg_desc = sorted(owner_deg, key=owner_deg.get, reverse=True)
    #         top = [owner_deg[deg_desc[0]], owner_deg[deg_desc[1]]]
    #         if top[0] == top[1]:
    #             # in event of tie, pick winner at random
    #             idx = np.random.choice(deg_desc[:2])
    #             return idx, owner_deg[idx]
    #         else:
    #             return deg_desc[0], top[0]

    # def demote(self, target):
    #     """
    #     Demote a node from being a hub. This means all subordinate nodes
    #     which recognised this node as a hub, are also updated.
    #     :param target: the target hub to remove
    #     """
    #     # for adjacent satellites vi of hub v
    #     adj_list = set(self.hubs[target])
    #     for adj in adj_list:
    #         adj_dat = self.hubs.node[adj]
    #         if target in adj_dat['owned_by']:
    #             # Remove target from list of owners of vadj
    #             adj_dat['owned_by'].remove(target)
    #
    #         # ?? this will cut off hub node from all its adjacent nodes
    #         self.hubs.remove_edge(target, adj)
    #     self.hubs.node[target]['owned_by'].remove(target)
    #     self.hubs.node[target]['hub'] = False
    #
    #     if self.hubs.degree(target) == 0:
    #         print '{0} was demoted and now has zero degree. Data: {1}'.format(target, self.hubs.node[target])

    # def migrate(self, vlist, from_v, to_v):
    #     """
    #     Update the edges and attributes of a list of nodes to
    #     reflect the change of ownership.
    #     :param vlist: the list of target vertices to move
    #     :param from_v: the origin hub vertex
    #     :param to_v: the destination hub vertex
    #     :return:
    #     """
    #     for vi in vlist:
    #         # avoid moving self
    #         print 'Moving {0} from {1} to {2}'.format(vi, from_v, to_v)
    #         if vi != from_v:
    #             try:
    #                 self.hubs.remove_edge(from_v, vi)
    #                 self.hubs.add_edge(to_v, vi)
    #                 self.hubs.node[vi]['owned_by'].remove(from_v)
    #                 self.hubs.node[vi]['owned_by'].add(to_v)
    #             except nx.NetworkXError as er:
    #                 print self.hubs[from_v]
    #                 raise er

    def exists(self, v):
        """
        The vertex v exists in cover graph.
        :param v: the target vertex.
        :return: True if the vertex exists in the cover graph.
        """
        return v in self.hubs

    def contains(self, vset):
        """
        The set of vertices exists in the cover graph.
        :param vset: the target set of vertices.
        :param no_isolates: include isolated vertices in output graph
        :return: True if the target set is contained by the cover graph.
        """
        return vset < set(self.hubs.nodes())

    def write(self, path, no_isolates=False):
        """
        Write cover graph to a file.
        :param path:
        :param no_isolates:
        :return:
        """
        # make a copy to modify
        gout = self.hubs.copy()

        # prune isolates if requested
        if no_isolates:
            iso_v = [v for v, dat in gout.nodes(data=True) if 'isolate' in dat and dat['isolate']]
            gout.remove_nodes_from(iso_v)

        for v in gout:
            # TODO networkx cannot write collections as attributes, so we're converting to a string for now.
            gout.node[v]['owned_by'] = ' '.join(gout.predecessors(v))
            gout.node[v].update(self.g.node[v])

        nx.write_graphml(gout, path)

    def write_mcl(self, path, no_isolates=False):
        # write MCL format cluster solution
        with open(path, 'w') as out:
            min_deg = 1 if no_isolates else 0
            # build a list of hubs and their degree, sort be descending degree for output
            hub_dict = dict((v, self.hubs.degree(v)) for v in self.hubs
                            if self.hubs.node[v]['hub'] and self.hubs.degree(v) >= min_deg)
            desc_deg = sorted(hub_dict, key=hub_dict.get, reverse=True)
            for vi in desc_deg:
                out.write('{0} {1}\n'.format(vi, ' '.join([str(vid) for vid in self.hubs[vi]])))
            print 'Final solution: {0} stars and {1} satellites'.format(len(hub_dict), sum(hub_dict.values()))

    def add_novel_adjacent(self, center, adj_list):
        """
        For a given center node, add all novel adjacent nodes supplied. Edges will
        be directed from the center hub to the node. Novelty is defined as nodes
        which currently do not exist in the cover graph.

        :param center: the vertex at the center (hub)
        :param adj_list: the list of adjacent nodes.
        """
        for vi in adj_list:
            if vi not in self.hubs:
                # add node with attributes
                self.hubs.add_node(vi, hub=False, analyzed=False)
            # directed edge from hub to node
            self.hubs.add_edge(center, vi, linked=False)

    def __str__(self):
        return str(self.hubs)

    def __repr__(self):
        return repr(self.hubs)

def relevance_by_node(g, desc=True):
    # TODO make this an ordered dict and dispense with one container
    rel_dict = dict((v, relevance(v, g)) for v in g.nodes() if len(g[v]) >= 1)
    desc_rel = sorted(dict((k, v) for k, v in rel_dict.iteritems() if v > 0.0), key=rel_dict.get, reverse=desc)
    return desc_rel, rel_dict

def ocluster(g, debug=False):
    print 'Graph contains {0} vertices, {1} edges'.format(g.order(), g.size())

    # prune self-loops before starting analysis
    for e in g.selfloop_edges():
        g.remove_edge(*e)
    print 'After removing self-loops, {0} edges remain'.format(g.size())

    desc_relevance, rel_dict = relevance_by_node(g, True)

    if debug:
        print 'Vertices sorted by decreasing relevance'
        for n, li in enumerate(desc_relevance, start=1):
            print '{0} {1} {2:.3f}'.format(n, li, rel_dict[li])

    C = Cover(g)

    # add isolated nodes first
    C.hubs.add_nodes_from(set(v for v in g.nodes_iter() if g.degree(v) == 0), hub=True, isolate=True)
    print 'There were {0} isolated vertices'.format(len(C.hubs))

    # traverse the set of vertices in descending order of relevance
    # looking for stars which cover new territory
    n_hub = 0
    for vi in desc_relevance:
        if vi not in C.hubs:
            # adding a new hub and adjacents
            print '{0} is a new hub'.format(vi)
            C.hubs.add_node(vi, hub=True, analyzed=False)
            C.add_novel_adjacent(vi, g[vi])
            n_hub += 1
        else:
            # check if adding this hub would cover new adjacent nodes
            # build set of adjacent nodes and compare against all
            vi_adj = set(g[vi])
            cover_set = set(C.hubs.nodes())
            if not vi_adj < cover_set:
                # set contained some extra (IE not a subset of cover_set)
                print '{0} covered another {1} adjacent'.format(vi, len(vi_adj - cover_set))
                # status is now a hub, add adjacents to it.
                C.hubs.node[vi]['hub'] = True
                C.add_novel_adjacent(vi, vi_adj)
                n_hub += 1

    print 'Finished step 1: {0} hubs and {1} vertices'.format(n_hub, len(C.hubs))

    # sort the ws-stars by descending degree
    deg_dict = dict((vi, g.degree(vi)) for vi in C.hubs if C.hubs.node[vi]['hub'])  # hubs only
    desc_degree = sorted(deg_dict, key=deg_dict.get, reverse=True)

    if debug:
        for n, v in enumerate(desc_degree, start=1):
            if v in rel_dict:
                print '{0} {1} {2:.3f} {3}'.format(n, v, rel_dict[v], deg_dict[v])
            else:
                print '{0} {1} NA {2}'.format(n, v, deg_dict[v])

    linked = {}
    demoted = []
    print 'Iterating through hubs by descending degree, eliminate useless ws_stars'
    # from largest hub to smallest, demote any hub which contains more shared
    # than unshared nodes (IE not useful).
    for vi in desc_degree:
        if deg_dict[vi] <= 0:
            continue

        print '\tFor hub {0} with degree {1}'.format(vi, deg_dict[vi])

        linked[vi] = set()
        v_adj = set(C.hubs[vi])
        for ui in v_adj:
            if C.hubs.node[ui]['hub'] and not C.hubs.node[ui]['analyzed']:
                if not C.is_useful(ui):
                    print '\t\tOverlapped hub {0} was not useful'.format(ui)
                    demoted.append(ui)
                    unshared = C.unshared(ui)
                    if ui in linked and len(linked[ui]) > 0:
                        print '\t\t\tForwarding {0} linked nodes from {1} to {2}'.format(len(linked[ui]), ui, vi)
                        unshared.update(linked[ui])
                        linked[ui].clear()
                    if len(unshared) > 0:
                        linked[vi].update(unshared)
                        print '\t\tThere were {0} unshared nodes to relink'.format(len(unshared))
                    u_adj = set(C.hubs[ui])
                    for wi in u_adj:
                        C.hubs.remove_edge(ui, wi)
                    C.hubs.node[ui]['hub'] = False
                else:
                    C.hubs.node[ui]['analyzed'] = True

    # Add re-link nodes that have become unconnected, attach them
    # to the grandparent hub of their demoted parent hub.
    print 'Relinking satellites'
    for li in linked:
        for vi in linked[li]:
            if not C.hubs.node[li]['hub']:
                raise RuntimeError('nodes linked to non-hub {0}'.format(li))
            C.hubs.add_edge(li, vi, linked=True)

    print 'Demoted {0} hub vertices, relinked {1} satellites'.format(len(demoted), len(linked))
    if debug:
        print 'Demoted hubs were: {0}'.format(demoted)
        print 'Relinked vertices were: {0}'.format(linked)

    return C

if __name__ == '__main__':

    parser = argparse.ArgumentParser('Cluster similarity graph using OClustR')
    parser.add_argument('-f', '--fmt', dest='format', choices=['graphml', 'mcl'], default='mcl', help='Output format')
    parser.add_argument('--debug', action='store_true', default=False, help='Enable debugging')
    parser.add_argument('--no-isolates', action='store_true', default=False, help='Ignore isolates in output')
    parser.add_argument('graph', nargs=1, help='Input graphml file')
    parser.add_argument('output', nargs=1, help='Output base file name')
    args = parser.parse_args()

    # read the graph
    g = nx.read_graphml(args.graph[0])

    # create cover of g
    cover = ocluster(g, args.debug)
    print cover

    if args.format == 'graphml':
        # write cover graph to file
        cover.write(args.output[0], args.no_isolates)
    else:
        # write MCL format cluster solution
        with open(args.output[0], 'w') as out:
            min_deg = 1 if args.no_isolates else 0
            # build a list of hubs and their degree, sort be descending degree for output
            hubs = cover.hubs
            hub_dict = dict((v, hubs.degree(v)) for v in hubs if hubs.node[v]['hub'] and hubs.degree(v) >= min_deg)
            desc_deg = sorted(hub_dict, key=hub_dict.get, reverse=True)
            for vi in desc_deg:
                out.write('{0} {1}\n'.format(vi, ' '.join(cover.hubs[vi])))
            print 'Final solution: {0} stars and {1} satellites'.format(len(hub_dict), sum(hub_dict.values()))
