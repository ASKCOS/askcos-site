/*
 * Graph utility functions
 */

class RetroGraph {
    constructor(nodes, edges, options) {
        this.nodes = new vis.DataSet(nodes, options);
        this.edges = new vis.DataSet(edges, options);
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
        let successors = [];
        this.edges.forEach(edge => {
            if (edge != null && edge.from === node) {
                successors.push(edge.to)
            }
        })
        return successors
    }
    getAllSuccessors(node) {
        // Retrieve all successors of the specified node
        let successors = [];
        this.edges.forEach(edge => {
            if (edge != null && edge.from === node) {
                successors.push(edge.to)
                successors.push(...this.getAllSuccessors(edge.to))
            }
        })
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
