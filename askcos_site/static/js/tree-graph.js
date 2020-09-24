var container = document.getElementsByClassName('container')[0];
container.classList.remove('container')
container.classList.add('container-fluid')
container.style.width = null;

function hideNetwork(n) {
    var networkDiv = document.querySelectorAll('.tree-graph')[n];
    var hideDiv = document.querySelectorAll('.hideNetwork')[n];
    var showDiv = document.querySelectorAll('.showNetwork')[n];
    networkDiv.style.display = 'none';
    hideDiv.style.display = 'none';
    showDiv.style.display = '';
}

function hideAllNetworks() {
    for (var n = 0; n < document.querySelectorAll('.tree-graph').length; n++) {
        hideNetwork(n)
    }
}

function showNetwork(n) {
    var networkDiv = document.querySelectorAll('.tree-graph')[n];
    var hideDiv = document.querySelectorAll('.hideNetwork')[n];
    var showDiv = document.querySelectorAll('.showNetwork')[n];
    networkDiv.style.display = '';
    hideDiv.style.display = '';
    showDiv.style.display = 'none';
}

function showAllNetworks() {
    for (var n = 0; n < document.querySelectorAll('.tree-graph').length; n++) {
        showNetwork(n)
    }
}

function colorOf(child) {
    if (child['ppg']) {
        if (child['as_reactant'] || child['as_product']) {
            return "#1B5E20" // green
        } else {
            return '#FFC400' // yellow
        }
    } else {
        if (child['as_reactant'] || child['as_product']) {
            return '#E65100' // orange
        } else {
            return '#B71C1C' // red
        }
    }
}

function loadNodeLinkGraph(data, showDetail = true) {
    /* Load tree in node link format into visjs and add visualization related attributes. */
    /* Add image data to original object */
    data.nodes.forEach(node => {
        if (node.type === 'chemical') {
            node.image = `/api/v2/draw/?smiles=${encodeURIComponent(node.smiles)}`;
            node.shape = 'image';
        }
    })
    let nodes
    if (showDetail) {
        /* For detail view, make a shallow copy of the nodes and add extra visual attributes */
        nodes = []
        data.nodes.forEach(node => {
            let newNode = Object.assign({}, node)
            if (newNode.type === 'chemical') {
                const buyableString = newNode.ppg !== 0 ? `$${newNode.ppg}/g` : 'not buyable';
                newNode.title = `${newNode.smiles}<br>${newNode.as_reactant} precedents as reactant<br>${newNode.as_product} precedents as product<br>${buyableString}`;
                newNode.borderWidth = 2;
                newNode.color = {
                    border: colorOf(newNode)
                }
            } else {
                newNode.label = `${newNode.num_examples} examples
FF score: ${Number(newNode.plausibility).toFixed(3)}
Template score: ${Number(newNode.template_score).toFixed(3)}`;
                newNode['font'] = {align: 'center'}
            }
            nodes.push(newNode)
        })
    } else {
        /* Otherwise, just display the original data */
        nodes = data.nodes;
    }
    return {
        nodes: new vis.DataSet(nodes),
        edges: new vis.DataSet(data.edges),
    }
}

function treeStats(tree) {
    let numReactions = 0
    let avgScore = 0
    let avgPlausibility = 0
    let minScore = 1.0
    let minPlausibility = 1.0
    for (node of tree.nodes) {
        if (node.type === 'reaction') {
            if (node.smiles.includes(app.settings.smiles)) {
                /* This is the first reaction step */
                tree.firstStepScore = node.template_score
            }
            numReactions += 1
            avgScore += node.template_score
            avgPlausibility += node.plausibility
            minScore = Math.min(minScore, node.template_score)
            minPlausibility = Math.min(minPlausibility, node.plausibility)
        }
    }
    avgScore /= numReactions
    avgPlausibility /= numReactions

    tree.numReactions = numReactions
    tree.avgScore = avgScore
    tree.avgPlausibility = avgPlausibility
    tree.minPlausibility = minPlausibility

    if ('score' in tree.graph) {
        tree.score = tree.graph.score
    }
    if ('cluster_id' in tree.graph) {
        tree.cluster_id = tree.graph.cluster_id
    }
    if ('depth' in tree.graph) {
        tree.depth = tree.graph.cluster_id
    }
}

function sortObjectArray(arr, prop, ascending) {
    arr.sort(function (a, b) {
        if (ascending) {
            return a[prop] - b[prop]
        } else {
            return b[prop] - a[prop]
        }
    })
}

function initializeNetwork(data, container, showDetail = true) {
    const options = showDetail ?  defaultVisOptions : condensedVisOptions
    return new vis.Network(container, data, options);
}

function buildTreeList() {
    /* Callback used by list view panel for drawing trees after panel creation */
    app.trees.forEach((tree, index) => {
        let elem = document.getElementById(`treeList-${index}`);
        let networkData = loadNodeLinkGraph(tree, false)
        let network = initializeNetwork(networkData, elem, false);
    })
}

const defaultVisOptions = {
    nodes: {
        color: {
            background: '#FFFFFF',
            border: '#000000',
        },
        shapeProperties: {
            useBorderWithImage: true,
            useImageSize: true,
        }
    },
    edges: {
        length: 1,
    },
    interaction: {
        dragNodes: false,
        dragView: true,
        hover: true,
        multiselect: false,
        selectConnectedEdges: false,
        tooltipDelay: 0,
        zoomView: true,
    },
    layout: {
        hierarchical: {
            direction: 'LR',
            levelSeparation: 250,
            nodeSpacing: 175,
            sortMethod: 'directed',
        }
    },
    physics: false,
};

const condensedVisOptions = {
    nodes: {
        color: {
            background: '#FFFFFF',
            border: '#000000',
        },
        shapeProperties: {
            useImageSize: true
        }
    },
    edges: {
        length: 1,
    },
    interaction: {
        dragNodes: false,
        dragView: true,
        hover: false,
        multiselect: false,
        selectable: false,
        selectConnectedEdges: false,
        tooltipDelay: 0,
        zoomView: true,
    },
    layout: {
        hierarchical: {
            direction: 'LR',
            levelSeparation: 200,
            nodeSpacing: 175,
            sortMethod: 'directed',
        }
    },
    clickToUse: true,
    physics: false,
};

var csrftoken = getCookie('csrftoken');

var app = new Vue({
    el: '#app',
    data: {
        resultId: "",
        numChemicals: 0,
        numReactions: 0,
        alltrees: [],
        settings: {},
        tbVersion: null,
        showInfoPanel: true,
        showListView: false,
        selected: null,
        currentTreeId: 0,
        networkData: {},
        cluster: false,
        currentClusterId: 0,
        sortOrderAscending: false,
        treeSortOption: 'numReactions',
        infoPanelOptions: {
            id: 'infoPanel',
            headerTitle: 'Info',
            headerControls: {size: 'sm'},
            position: {my: 'left-top', at: 'left-top', of: '#graph'},
            panelSize: {width: 500, height: 500},
        },
        detailPanelOptions: {
            id: 'detailPanel',
            headerTitle: 'Node Details',
            headerControls: {size: 'sm'},
            position: {my: 'right-top', at: 'right-top', of: '#graph'},
            panelSize: {width: 500, height: 500},
        },
        listPanelOptions: {
            id: 'listPanel',
            headerTitle: 'List View',
            headerControls: {size: 'sm'},
            position: {my: 'center-top', at: 'center-top', of: '#graph'},
            panelSize: {width: () => window.innerWidth, height: 600},
            callback: buildTreeList,
        },
    },
    mounted: function () {
        this.resultId = this.$el.getAttribute('data-id');
        this.getResult(this.resultId);
    },
    methods: {
        getResult: function (id) {
            fetch(`/api/v2/results/${id}/tree/`)
                .then(resp => resp.json())
                .then(json => {
                    var result = json['result'];
                    var stats = result['result']['status'];
                    var trees = result['result']['paths'];
                    this.numChemicals = stats[0];
                    this.numReactions = stats[1];
                    this.alltrees = trees;
                    this.settings = result['settings'];
                    if (!!this.settings.buyables_source
                        && (this.settings.buyables_source.includes(null) || this.settings.buyables_source.includes(''))) {
                        const to_remove = [null, '']
                        this.settings.buyables_source = this.settings.buyables_source.filter(item => !to_remove.includes(item))
                        this.settings.buyables_source.push('(no source)')
                    }
                    // If version is not present in the result, then it is version 1
                    this.tbVersion = result['settings']['version'] || 1;
                    this.networkContainer = document.getElementById('graph')
                    if (this.trees.length) {
                        this.allTreeStats()
                        this.setDefaultSortOrder()
                        this.sortTrees()
                        this.buildTree();
                    }
                    setTimeout(() => {
                        document.querySelector('#splash').classList.replace("d-flex", "d-none")
                    }, 1000)
                })
        },
        buildTree: function () {
            this.networkData = loadNodeLinkGraph(this.trees[this.currentTreeId], true)
            this.network = initializeNetwork(this.networkData, this.networkContainer, true);
            this.network.on('selectNode', function (params) {
                app.showNode(params.nodes[0])
            });
            this.network.on('deselectNode', this.clearSelection);
        },
        sortTrees: function () {
            sortObjectArray(this.trees, this.treeSortOption, this.sortOrderAscending)
            this.currentTreeId = 0
            this.buildTree()
            if (this.showListView) {
                buildTreeList()
            }
        },
        setDefaultSortOrder: function() {
            this.sortOrderAscending = ['numReactions', 'depth'].includes(this.treeSortOption)
        },
        nextTree: function () {
            if (this.currentTreeId < this.trees.length - 1) {
                this.clearSelection()
                this.currentTreeId += 1
                this.buildTree()
            }
        },
        prevTree: function () {
            if (this.currentTreeId > 0) {
                this.clearSelection()
                this.currentTreeId -= 1
                this.buildTree()
            }
        },
        firstTree: function () {
            this.clearSelection()
            this.currentTreeId = 0
            this.buildTree()
        },
        lastTree: function () {
            this.clearSelection()
            this.currentTreeId = this.trees.length - 1
            this.buildTree()
        },
        nextCluster: function () {
            if (this.currentClusterId < this.maxClusterId) {
                this.clearSelection()
                this.currentClusterId += 1
                this.currentTreeId = 0
                this.buildTree()
                if (this.showListView) {
                    buildTreeList()
                }
            }
        },
        prevCluster: function () {
            if (this.currentClusterId > this.minClusterId) {
                this.clearSelection()
                this.currentClusterId -= 1
                this.currentTreeId = 0
                this.buildTree()
                if (this.showListView) {
                    buildTreeList()
                }
            }
        },
        firstCluster: function () {
            this.clearSelection()
            this.currentClusterId = this.minClusterId
            this.currentTreeId = 0
            this.buildTree()
            if (this.showListView) {
                buildTreeList()
            }
        },
        lastCluster: function () {
            this.clearSelection()
            this.currentClusterId = this.maxClusterId
            this.currentTreeId = 0
            this.buildTree()
            if (this.showListView) {
                buildTreeList()
            }
        },
        allTreeStats: function () {
            this.trees.forEach(treeStats)
        },
        banItem: function () {
            var nodeId = this.network.getSelectedNodes();
            if (!nodeId.length) {
                return
            }
            var desc = prompt("Please enter a reason (for your records only)", "no reason");
            var node = this.networkData.nodes.get(nodeId[0]);
            if (node['type'] === 'chemical') {
                var url = '/api/v2/banlist/chemicals/'
                var speciesName = 'chemical'
            } else {
                var url = '/api/v2/banlist/reactions/'
                var speciesName = 'reaction'
            }
            const body = {
                smiles: node.smiles,
                description: desc,
            };
            fetch(url, {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                    'X-CSRFToken': csrftoken
                },
                body: JSON.stringify(body)
            })
                .then(resp => resp.json())
                .then(json => {
                    const datetime = dayjs(json.created).format('MMMM D, YYYY h:mm A');
                    alert(`Banned ${speciesName} ${node.smiles} at ${datetime}`)
                });
        },
        showNode: function (nodeId) {
            var node = this.networkData.nodes.get(nodeId)
            this.selected = node
            if (node.type == 'chemical' && !!!node.source) {
                fetch(`/api/v2/buyables/?q=${encodeURIComponent(node.smiles)}&source=${this.settings.buyables_source}`)
                    .then(resp => resp.json())
                    .then(json => {
                        if (json.result.length) {
                            this.networkData.nodes.update({id: node.id, source: json.result[0].source})
                            this.$set(this.selected, 'source', json.result[0].source)
                        }
                    })
            }
        },
        clearSelection: function () {
            this.selected = null
        },
    },
    computed: {
        trees: function () {
            if (this.cluster) {
                return this.alltrees.filter(tree => tree.cluster_id === this.currentClusterId)
            } else {
                return this.alltrees
            }
        },
        maxClusterId: function () {
            if (!!this.alltrees.length && 'cluster_id' in this.alltrees[0]) {
                return Math.max(...this.alltrees.map(tree => tree.cluster_id))
            } else {
                return 0
            }
        },
        minClusterId: function () {
            if (!!this.alltrees.length && 'cluster_id' in this.alltrees[0]) {
                return Math.min(...this.alltrees.map(tree => tree.cluster_id))
            } else {
                return 0
            }
        },
    },
    watch: {
        cluster: function () {
            this.currentClusterId = this.minClusterId;
            this.currentTreeId = 0;
            this.buildTree()
            if (this.showListView) {
                buildTreeList()
            }
        },
    },
    delimiters: ['%%', '%%'],
});
