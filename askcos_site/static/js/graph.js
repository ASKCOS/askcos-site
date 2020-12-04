/*
 * Graph utility functions
 */

class RetroGraph {
    constructor(nodes, edges, options) {
        this.nodes = new vis.DataSet(options);
        this.edges = new vis.DataSet(options);
        this.edges.on('*', this.addEdgeCallback);

        this.succ = {};  // Precomputed successors for fast retrieval
        this.nodes.add(nodes)
        this.edges.add(edges)
    }
    addEdgeCallback(event, properties, senderId) {
        // Update the successors object when a new edge is added
        switch (event) {
            case 'add':
                for (let id of properties.items) {
                    let edge = this.edges.get(id)
                    this.succ[edge.from] = this.succ[edge.from] || {}
                    this.succ[edge.from][edge.to] = edge.id
                }
                break;
            case 'remove':
                for (let id of properties.items) {
                    let edge = this.edges.get(id)
                    delete this.succ[edge.from][edge.to]
                }
                break;
            default:
                throw `Cannot handle ${event} event!`;
        }
    }
    getPredecessors(node) {
        // Retrieve immediate predecessors of the specified node
        let predecessors = [];
        this.edges.forEach(edge => {
            if (edge != null && edge.to === node) {
                predecessors.push(edge.from)
            }
        })
        return predecessors
    }
    getSuccessors(node) {
        // Retrieve immediate successors of the specified node
        return Object.keys(this.succ[node])
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
        let dangling = this.edges.get({
            filter: edge => !nodeIds.includes(edge.from) || !nodeIds.includes(edge.to)
        });
        this.edges.remove(dangling)
    }
    *getChemPaths(node, chemPath, maxDepth) {
        let successors = this.getSuccessors(node)
        if (successors.length === 0 || maxDepth !== undefined && chemPath.length >= maxDepth) {
            yield [node]
        } else {
            for (let rxn of successors) {
                for (let subPath of this.getRxnPaths(rxn, chemPath.concat([node]))) {
                    subPath.push(rxn);
                    yield subPath
                }
            }
        }
    }
    *getRxnPaths(node, chemPath, maxDepth) {
        let successors = this.getSuccessors(node)
        if (chemPath.some(item => successors.includes(item))) {
            return
        }
        for (let pathCombo of product(this.getChemPaths(c))) {
            let subPath = concat(pathCombo)
            subPath.push(node)
            yield subPath
        }
    }
    *getPaths(maxDepth, maxTrees, validatePaths, isTerminal) {
        let numPaths = 0;
        for (let path of this.getChemPaths(NIL_UUID, [], maxDepth)) {
            if (maxTrees !== undefined && numPaths >= maxTrees) {
                break
            }
            if (validatePaths && isTerminal instanceof Function && isTerminal(path)) {
                numPaths += 1
                yield path
            }
        }
    }
}
