/*
 * Graph utility functions
 */

class RetroGraph {
    constructor(nodes, edges, source, target, visOptions) {
        this.nodes = new vis.DataSet(visOptions);
        this.edges = new vis.DataSet(visOptions);
        this._source = source || 'from'
        this._target = target || 'to'
        this._succ = {};  // Precomputed successors for fast retrieval

        let graph = this
        this.edges.on('*', function(event, properties, senderId) {
            // Update the successors object when a new edge is added
            if (senderId === 'clear') {
                return  // _succ property will be manually reset
            }
            switch (event) {
                case 'add':
                    for (let id of properties.items) {
                        let edge = graph.edges.get(id)
                        graph._succ[edge[graph._source]] = graph._succ[edge[graph._source]] || {}
                        graph._succ[edge[graph._source]][edge[graph._target]] = edge.id
                    }
                    break;
                case 'remove':
                    for (let obj of properties.oldData) {
                        delete graph._succ[obj[graph._source]][obj[graph._target]]
                    }
                    break;
                default:
                    throw `Cannot handle ${event} event!`;
            }
        });

        if (nodes !== undefined) this.nodes.add(nodes);
        if (edges !== undefined) this.edges.add(edges);
    }
    clear() {
        this.nodes.clear('clear')
        this.edges.clear('clear')
        this._succ = {};
    }
    toJSON() {
        return {
            nodes: this.nodes.get(),
            edges: this.edges.get(),
        }
    }
    getPredecessors(node) {
        // Retrieve immediate predecessors of the specified node
        let predecessors = [];
        this.edges.forEach(edge => {
            if (edge != null && edge[this._target] === node) {
                predecessors.push(edge[this._source])
            }
        })
        return predecessors
    }
    getSuccessors(node) {
        // Retrieve immediate successors of the specified node
        let successors = this._succ[node]
        if (!!successors) {
            return Object.keys(successors)
        } else {
            return []
        }
    }
    getAllSuccessors(node) {
        // Retrieve all successors of the specified node
        let successors = this.getSuccessors(node);
        for (let succ of this.getSuccessors(node)) {
            successors.push(...this.getSuccessors(succ))
        }
        return successors
    }
    removeAllSuccessors(node) {
        // Remove all successors of the specified node
        this.nodes.remove(this.getAllSuccessors(node));
        this.trimDanglingEdges()
    }
    trimDanglingEdges() {
        // Remove any edges which are only connected on one end
        let nodeIds = this.nodes.getIds();
        let dangling = this.edges.getIds({
            filter: edge => !nodeIds.includes(edge[this._source]) || !nodeIds.includes(edge[this._target])
        });
        this.edges.remove(dangling)
    }
    getChemPaths(node, chemPath, maxDepth) {
        let successors = this.getSuccessors(node)
        let paths = [];
        if (successors.length === 0 || maxDepth !== undefined && chemPath.length >= maxDepth) {
            paths.push({nodes: [node], edges: [], depth: chemPath.length})
        } else {
            for (let rxn of successors) {
                for (let subPath of this.getRxnPaths(rxn, [...chemPath, node], maxDepth)) {
                    subPath.nodes.push(node);
                    subPath.edges.push(this._succ[node][rxn])
                    paths.push(subPath);
                }
            }
        }
        return paths
    }
    getRxnPaths(node, chemPath, maxDepth) {
        let successors = this.getSuccessors(node)
        let paths = [];
        if (chemPath.some(item => successors.includes(item))) {
            return paths
        }
        let subPaths = successors.map(c => this.getChemPaths(c, chemPath, maxDepth))
        for (let pathCombo of product(subPaths)) {
            let subPath = mergePaths(pathCombo)
            subPath.nodes.push(node)
            successors.forEach(c => subPath.edges.push(this._succ[node][c]))
            paths.push(subPath)
        }
        return paths
    }
    getPaths(root, maxDepth, maxTrees, validatePaths = false) {
        let paths = [];
        for (let path of this.getChemPaths(root, [], maxDepth)) {
            if (maxTrees !== undefined && paths.length >= maxTrees) {
                break
            }
            if (!validatePaths || validatePaths && this.isValid(path)) {
                paths.push(path)
            }
        }
        return paths
    }
    isValid(path) {
        for (let node of path.nodes) {
            if (!this._succ[node] && !this.nodes.get(node).terminal) {
                return false
            }
        }
        return true
    }
}

function product(sets) {
    // Cartesian product generator, similar to itertools.product in Python
    // Performance tested to be faster than most alternatives
    let max = sets.length - 1;
    let lens = sets.map(set => set.length);
    let results = [];
    let combo = [];
    function build(n) {
        let set = sets[n];
        let len = lens[n];
        if (n === max) {
            for (let i = 0; i < len; ++i) {
                combo[n] = set[i];
                results.push([...combo]);
            }
        } else {
            for (let i = 0; i < len; ++i) {
                combo[n] = set[i];
                build(n + 1);
            }
        }
        combo.pop()
    }
    build(0)
    return results
}

function mergePaths(paths) {
    // Combine multiple path objects
    let result = {'nodes': [], 'edges': [], 'depth': 0}
    for (let path of paths) {
        result.nodes.push(...path.nodes)
        result.edges.push(...path.edges)
        result.depth = Math.max(result.depth, path.depth)
    }
    return result
}
