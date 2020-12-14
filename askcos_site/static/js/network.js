var container = document.getElementsByClassName('container')[0];
container.classList.remove('container')
container.classList.add('container-fluid')
container.style.width=null;

const NIL_UUID = '00000000-0000-0000-0000-000000000000'

function updateObj(dest, src) {
    // take properties of src and overwrite matching properties of dest
    // ignores properties in src if they do not exist in dest
    // modifies dest object in place
    for (let p in src) {
        if (src.hasOwnProperty(p) && dest.hasOwnProperty(p)) {
            if ( typeof dest[p] === 'object' && !Array.isArray(dest[p]) ) {
                updateObj(dest[p], src[p])
            } else {
                dest[p] = src[p];
            }
        }
    }
}

function num2str(n, len, exp = false) {
    if (len === undefined) {
        return n === undefined || isNaN(n) ? 'N/A' : n.toString();
    } else if (exp) {
        return n === undefined || isNaN(n) ? 'N/A' : n.toExponential(len);
    } else {
        return n === undefined || isNaN(n) ? 'N/A' : n.toFixed(len);
    }
}

function edgeScaling(min, max, total, value) {
    if (value >= 0.25) {
        return 1.0
    } else {
        return 16 * value * value
    }
}

/* DnD */
function disable_dragstart_handler(e) {
    e.preventDefault();
}

function clusteredit_dragover_handler(event) {
    event.preventDefault(); // important
    event.stopPropagation();
    event.dataTransfer.dropEffect = 'move'; // important
}

function clusteredit_dragenter_handler(event) {
    event.preventDefault(); // important
    event.stopPropagation();
    event.target.classList.add('dragover');
}

function clusteredit_dragleave_handler(event) {
    event.target.classList.remove('dragover');
}

const tbSettingsDefault = {
    quick: "normal",
    version: 1,
    maxDepth: 5,
    maxBranching: 20,
    expansionTime: 60,
    maxChemicals: null,
    maxReactions: null,
    maxIterations: null,
    buyableLogic: 'and',
    maxPPGLogic: 'none',
    maxPPG: 100,
    maxScscoreLogic: 'none',
    maxScscore: 0,
    chemicalPropertyLogic: 'none',
    chemicalPropertyC: 0,
    chemicalPropertyN: 0,
    chemicalPropertyO: 0,
    chemicalPropertyH: 0,
    chemicalPopularityLogic: 'none',
    chemicalPopularityReactants: 0,
    chemicalPopularityProducts: 0,
    buyablesSource: [],
    buyablesSourceAll: true,
    returnFirst: false,
    maxTrees: 500,
    templateSet: "reaxys",
    templateSetVersion: 1,
    precursorScoring: "RelevanceHeuristic",
    numTemplates: 1000,
    maxCumProb: 0.999,
    minPlausibility: 0.1,
    allowSelec: true,
    attributeFilter: []
};

const visjsOptionsDefault = {
    edges: {
        length: 1
    },
    nodes: {
        mass: 1,
        size: 25,
        font: {
            size: 14,
        },
        color: {
            border: '#000000',
            background: '#FFFFFF'
        },
        shapeProperties: {
            useBorderWithImage: true
        }
    },
    layout: {
        hierarchical: {
            enabled:false,
            levelSeparation: 150,
            nodeSpacing: 100,
            treeSpacing: 200,
            blockShifting: true,
            edgeMinimization: true,
            parentCentralization: true,
            direction: 'UD',
            sortMethod: 'directed',
            shakeTowards: 'roots',
        }
    },
    interaction:{
        dragNodes:true,
        dragView: true,
        hideEdgesOnDrag: false,
        hideNodesOnDrag: false,
        hover: false,
        hoverConnectedEdges: true,
        keyboard: {
        enabled: false,
        speed: {x: 10, y: 10, zoom: 0.02},
        bindToWindow: true
        },
        multiselect: true,
        navigationButtons: false,
        selectable: true,
        selectConnectedEdges: true,
        tooltipDelay: 300,
        zoomView: true
    },
    physics:{
        enabled: true,
        barnesHut: {
            gravitationalConstant: -2000,
            centralGravity: 0.3,
            springLength: 95,
            springConstant: 0.04,
            damping: 0.09,
            avoidOverlap: 0
        },
        maxVelocity: 50,
        minVelocity: 0.1,
        solver: 'barnesHut',
        stabilization: {
            enabled: true,
            iterations: 1000,
            updateInterval: 100,
            onlyDynamicEdges: false,
            fit: true
        },
        timestep: 0.5,
        adaptiveTimestep: true
    }
};

function getVisjsUserOptions(obj) {
    // extract user adjustable options from the full visjs options object
    return {
        nodes: {
            mass: obj.nodes.mass,
            size: obj.nodes.size,
            font: {
                size: obj.nodes.font.size,
            },
        },
        layout: {
            hierarchical: {
                enabled: obj.layout.hierarchical.enabled,
                levelSeparation: obj.layout.hierarchical.levelSeparation,
                direction: obj.layout.hierarchical.direction,
            }
        },
        physics: {
            barnesHut: {
                springConstant: obj.physics.barnesHut.springConstant,
            }
        }
    }
}

const ippSettingsDefault = {
    allowCluster: true,
    filterReactingAtoms: false,
    allowResolve: false,
    isHighlightAtom: true,
    reactionLimit: 5,
    sortingCategory: "retroScore",
    sortOrderAscending: false,
    selectivityModel: "qm_GNN",
    clusterOptions: {
        feature: 'original',
        fingerprint:'morgan',
        fpRadius: 1, fpBits: 512,
        cluster_method: 'kmeans',
        isAlternatingColor: false,
    },
};

var app = new Vue({
    el: '#app',
    data: {
        isAuth: isAuth,
        window: {
            width: 0,
            height: 0,
        },
        target: '',
        dataGraph: new RetroGraph(),
        dispGraph: new RetroGraph(),
        templateSets: {},
        templateAttributes: {},
        buyablesSources: [],
        templateNumExamples: {},
        allowCluster: ippSettingsDefault.allowCluster,
        filterReactingAtoms: ippSettingsDefault.filterReactingAtoms,
        refreshFilter: 1,
        allowResolve: ippSettingsDefault.allowResolve,
        showSettingsModal: false,
        showLoadModal: false,
        showDownloadModal: false,
        showClusterPopoutModal: false,
        showClusterEditModal: false,
        showAddNewPrecursorModal: false,
        downloadName: "network.json",
        tb: {
            settings: JSON.parse(JSON.stringify(tbSettingsDefault)),
            modes: {
                quickest: {
                    maxDepth: 4,
                    expansionTime: 30,
                    returnFirst: true,
                    maxBranching: 20,
                    numTemplates: 100,
                    maxCumProb: 0.995,
                    minPlausibility: 0.001
                },
                shallow: {
                    maxDepth: 4,
                    expansionTime: 30,
                    returnFirst: false,
                    maxBranching: 20,
                    numTemplates: 100,
                    maxCumProb: 0.995,
                    minPlausibility: 0.75
                },
                normal: {
                    maxDepth: 5,
                    expansionTime: 60,
                    returnFirst: false,
                    maxBranching: 20,
                    numTemplates: 1000,
                    maxCumProb: 0.999,
                    minPlausibility: 0.1
                },
                deep: {
                    maxDepth: 6,
                    expansionTime: 120,
                    returnFirst: false,
                    maxBranching: 25,
                    numTemplates: 1000,
                    maxCumProb: 0.9909,
                    minPlausibility: 0.01
                }
            },
            redirectToGraph: false
        },
        clusterPopoutModalData: {
            optionsDisplay : {
                showScore: false,
                showSCScore: false,
                showNumExample: true,
                showTemplateScore: false,
                showPlausibility: true,
                showClusterId: false,
            },
        },
        clusterEditModalData: {
            optionsDisplay : {
                showScore: false,
                showNumExample: false,
                showTemplateScore: false,
                showPlausibility: false,
                showClusterId: false,
            },
        },
        addNewPrecursorModal: {},
        clusterOptions: JSON.parse(JSON.stringify(ippSettingsDefault.clusterOptions)),
        selected: null,
        isHighlightAtom: ippSettingsDefault.isHighlightAtom,
        reactionLimit: ippSettingsDefault.reactionLimit,
        selectivityModel: ippSettingsDefault.selectivityModel,
        sortingCategory: ippSettingsDefault.sortingCategory,
        sortOrderAscending: ippSettingsDefault.sortOrderAscending,
        networkOptions: JSON.parse(JSON.stringify(visjsOptionsDefault)),
        recompute: 0,  // Dummy property to trigger computed properties depending on dataGraph or dispGraph
    },
    beforeMount: function() {
        this.enableResolve = this.$el.querySelector('[ref="enableResolve"]').checked;
        this.allowResolve = this.$el.querySelector('[ref="allowResolve"]').checked;
    },
    created: function() {
        window.addEventListener('resize', this.handleResize);
        this.handleResize();

        this.loadNetworkOptions()
        this.loadTarget()
        this.loadTbSettings()
        this.loadIppSettings()
        var urlParams = new URLSearchParams(window.location.search);
        let urlTarget = urlParams.get('target')
        if (urlTarget) {
            this.target = urlTarget
        }
        let run = urlParams.get('run')
        if (run && JSON.parse(run)) {
            this.changeTarget()
        }
        let loadTreeBuilder = urlParams.get('tb')
        let numTrees = urlParams.get('view')
        if (loadTreeBuilder) {
            this.loadFromTreeBuilder(loadTreeBuilder, numTrees)
        }
        fetch('/api/v2/template/sets/')
            .then(resp => resp.json())
            .then(json => {
                this.templateAttributes = json.attributes
                for (templateSet of json.template_sets) {
                    this.templateSets[templateSet] = { versions: [] }
                    fetch('/api/v2/retro/models/?template_set='+templateSet)
                        .then(resp => resp.json())
                        .then(json => {
                            if (json.versions) {
                                this.templateSets[json.request.template_set] = json.versions.map(x => Number(x))    
                            }
                        })
                }
            })
        fetch('/api/v2/buyables/sources/')
            .then(resp => resp.json())
            .then(json => {this.buyablesSources = json.sources})
    },
    mounted: function() {
        setTimeout(() => {
            document.querySelector('#splash').classList.replace("d-flex", "d-none")
        }, 1000)
    },
    destroyed: function() {
        window.removeEventListener('resize', this.handleResize);
    },
    methods: {
        initializeNetwork(data) {
            var container = document.getElementById('network');
            this.network = new vis.Network(container, data, this.networkOptions);
            this.network.on("beforeDrawing",  function(ctx) {
                ctx.save();
                ctx.setTransform(1, 0, 0, 1, 0, 0);
                ctx.fillStyle = '#ffffff';
                ctx.fillRect(0, 0, ctx.canvas.width, ctx.canvas.height)
                ctx.restore();
            })
        },
        centerGraph() {
            if (!!this.network) {
                this.network.fit()
            }
        },
        saveNetworkOptions() {
            if (!storageAvailable('localStorage')) return
            localStorage.setItem('visjsOptions', encodeURIComponent(JSON.stringify(getVisjsUserOptions(this.networkOptions))))
        },
        saveTarget() {
            if (!storageAvailable('localStorage')) return
            localStorage.setItem('target', this.target)
        },
        saveTbSettings() {
            if (!storageAvailable('localStorage')) return
            localStorage.setItem('tbSettings', encodeURIComponent(JSON.stringify(this.tb.settings)))
        },
        saveIppSettings() {
            if (!storageAvailable('localStorage')) return
            const obj = {
                allowCluster: this.allowCluster,
                filterReactingAtoms: ippSettingsDefault.filterReactingAtoms,
                allowResolve: this.allowResolve,
                isHighlightAtom: this.isHighlightAtom,
                reactionLimit: this.reactionLimit,
                selectivityModel: this.selectivityModel,
                sortingCategory: this.sortingCategory,
                sortOrderAscending: this.sortOrderAscending,
                clusterOptions: this.clusterOptions,
            }
            localStorage.setItem('ippSettings', encodeURIComponent(JSON.stringify(obj)))
        },
        loadNetworkOptions() {
            if (!storageAvailable('localStorage')) return
            const settings = localStorage.getItem('visjsOptions')
            if (!settings) return
            const obj = JSON.parse(decodeURIComponent(settings))
            const userOptions = getVisjsUserOptions(this.networkOptions)
            updateObj(userOptions, obj)
            updateObj(this.networkOptions, userOptions)
        },
        loadTarget() {
            if (!storageAvailable('localStorage')) return
            const target = localStorage.getItem('target')
            if (!target) return
            this.target = target
        },
        loadTbSettings() {
            if (!storageAvailable('localStorage')) return
            const settings = localStorage.getItem('tbSettings')
            if (!settings) return
            const obj = JSON.parse(decodeURIComponent(settings))
            updateObj(this.tb.settings, obj)
        },
        loadIppSettings() {
            if (!storageAvailable('localStorage')) return
            const settings = localStorage.getItem('ippSettings')
            if (!settings) return
            const obj = JSON.parse(decodeURIComponent(settings))
            updateObj(this, obj)
        },
        handleResize: function() {
            this.window.width = window.innerWidth;
            this.window.height = window.innerHeight;
        },
        tbSettings(mode) {
            if ((mode != "quickest") && (mode != "shallow") && (mode != "normal") && (mode != "deep")) {
                return
            }
            this.tb.settings.quick = mode
            this.tb.settings.maxDepth = this.tb.modes[mode].maxDepth
            this.tb.settings.expansionTime = this.tb.modes[mode].expansionTime
            this.tb.settings.returnFirst = this.tb.modes[mode].returnFirst
            this.tb.settings.maxBranching = this.tb.modes[mode].maxBranching
            this.tb.settings.numTemplates = this.tb.modes[mode].numTemplates
            this.tb.settings.maxCumProb = this.tb.modes[mode].maxCumProb
            this.tb.settings.minPlausibility = this.tb.modes[mode].minPlausibility
        },
        isTbQuickSettingsMode(mode) {
            if (this.tb.settings.maxDepth != this.tb.modes[mode].maxDepth) return false
            if (this.tb.settings.expansionTime != this.tb.modes[mode].expansionTime) return false
            if (this.tb.settings.returnFirst != this.tb.modes[mode].returnFirst) return false
            if (this.tb.settings.maxBranching != this.tb.modes[mode].maxBranching) return false
            if (this.tb.settings.numTemplates != this.tb.modes[mode].numTemplates) return false
            if (this.tb.settings.maxCumProb != this.tb.modes[mode].maxCumProb )return false
            if (this.tb.settings.minPlausibility != this.tb.modes[mode].minPlausibility) return false
            return true
        },
        resetSettings() {
            this.$set(this.tb, 'settings', JSON.parse(JSON.stringify(tbSettingsDefault)));
            this.$set(this, 'networkOptions', JSON.parse(JSON.stringify(visjsOptionsDefault)));
            for (let key in JSON.parse(JSON.stringify(ippSettingsDefault))) {
                this.$set(this, key, ippSettingsDefault[key])
            }
            this.saveTbSettings();
            this.saveNetworkOptions();
            this.saveIppSettings();
        },
        sendTreeBuilderJob() {
            if (!isAuth) {
                alert('Error: must be logged in to start tree builder')
                return
            }
            if (this.tb.settings.name === '') {
                this.tb.settings.name = this.target
            }
            this.validatesmiles(this.target, !this.allowResolve)
            .then(isvalidsmiles => {
                if (isvalidsmiles) {
                    return this.target
                } else {
                    return this.resolveChemName(this.target)
                }
            })
            .then(smiles => {
                this.target = smiles
                this.mctsTreeBuilderAPICall()
            })
        },
        mctsTreeBuilderAPICall: function() {
            var url = '/api/v2/tree-builder/'
            this.saveNetworkOptions()
            this.saveTbSettings()
            this.saveIppSettings()
            this.saveTarget()
            var description = this.tb.settings.name ? this.tb.settings.name : this.target
            var body = {
                description: description,
                smiles: this.target,
                version: this.tb.settings.version,
                template_set: this.tb.settings.templateSet,
                template_prioritizer_version: this.tb.settings.templateSetVersion,
                max_depth: this.tb.settings.maxDepth,
                max_branching: this.tb.settings.maxBranching,
                expansion_time: this.tb.settings.expansionTime,
                max_chemicals: this.tb.settings.maxChemicals,
                max_reactions: this.tb.settings.maxReactions,
                max_iterations: this.tb.settings.maxIterations,
                buyable_logic: this.tb.settings.buyableLogic,
                max_ppg_logic: this.tb.settings.maxPPGLogic,
                max_ppg: this.tb.settings.maxPPG,
                max_scscore_logic: this.tb.settings.maxScscoreLogic,
                max_scscore: this.tb.settings.maxScscore,
                num_templates: this.tb.settings.numTemplates,
                max_cum_prob: this.tb.settings.maxCumProb,
                filter_threshold: this.tb.settings.minPlausibility,
                return_first: this.tb.settings.returnFirst,
                max_trees: this.tb.settings.maxTrees,
                store_results: true,
                chemical_property_logic: this.tb.settings.chemicalPropertyLogic,
                max_chemprop_c: this.tb.settings.chemicalPropertyC,
                max_chemprop_n: this.tb.settings.chemicalPropertyN,
                max_chemprop_o: this.tb.settings.chemicalPropertyO,
                max_chemprop_h: this.tb.settings.chemicalPropertyH,
                chemical_popularity_logic: this.tb.settings.chemicalPopularityLogic,
                min_chempop_reactants: this.tb.settings.chemicalPopularityReactants,
                min_chempop_products: this.tb.settings.chemicalPopularityProducts,
                json_format: 'nodelink',
            }
            if (!this.tb.settings.buyablesSourceAll) {
                body.buyables_source = this.tb.settings.buyablesSource
            }
            fetch(url, {
                method: 'POST', 
                headers: {
                    'Content-Type': 'application/json', 
                    'X-CSRFToken': getCookie('csrftoken')
                }, 
                body: JSON.stringify(body)
            })
                .then(resp => resp.json())
                .then(json => {
                    if (json.error) {
                        alert('Error: could not start tree builder. Try again later.')
                        return
                    }
                    else {
                        this.tb.taskID = json.task_id
                        this.tb.poll = setTimeout(() => this.pollForTbResult(), 1000)
                        notificationOptions = {
                            requireInteraction: true,
                            body: "The job will run in the background. You will see a new notification when the job completes."
                        }
                        app = this
                        this.makeNotification("Tree builder job submitted!", notificationOptions, (event) => {
                            this.close()
                        })
                    }
                })
        },
        makeNotification(title, options, callback) {
            if (!("Notification" in window)) {
                alert("This browser does not support desktop notifications! Notifications about tree builder submission and completion will not show.")
            }

            // Let's check whether notification permissions have already been granted
            else if (Notification.permission === "granted") {
                // If it's okay let's create a notification
                var notification = new Notification(title, options);
                app = this
                notification.onclick = callback
            }

            // Otherwise, we need to ask the user for permission
            else if (Notification.permission !== "denied") {
                Notification.requestPermission().then(function (permission) {
                // If the user accepts, let's create a notification
                if (permission === "granted") {
                    var notification = new Notification(title, options);
                    app = this
                    notification.onclick = callback
                }
                });
            }
        },
        pollForTbResult() {
            fetch('/api/v2/celery/task/'+this.tb.taskID)
                .then(resp => resp.json())
                .then(json => {
                    notificationOptions = {
                        requireInteraction: true,
                    }
                    if (json.complete) {
                        notifyMessage = "Tree builder job complete! Click to view results in a new tab."
                        app = this
                        notifyCallback = function(event) {
                            event.preventDefault(); // prevent the browser from focusing the Notification's tab
                            if (app.tb.redirectToGraph) {
                                window.open('/retro/network/?view=25&tb='+app.tb.taskID, '_blank')
                                this.close()
                            }
                            else {
                                window.open('/view-tree-graph/?id='+app.tb.taskID, '_blank')
                                this.close()
                            }
                        }
                        notificationOptions.body = "Click here to open a new tab with results."
                        this.makeNotification("Tree builder results", notificationOptions, notifyCallback)
                        this.tb.poll = null
                    }
                    else if (json.failed) {
                        notificationOptions.body = "Job failed. Try submitting a new job."
                        this.makeNotification("Tree builder results", notificationOptions, (event) => {this.close()})
                    }
                    else {
                        setTimeout(() => this.pollForTbResult(), 1000)
                    }
                })
                .catch(error => {
                    if (error instanceof TypeError) {
                        console.log('Unable to fetch tree builder results due to connection error. Will keep trying.')
                        setTimeout(() => this.pollForTbResult(), 2000)
                    } else {
                        console.error('There was a problem fetching results:', error);
                    }
                });
        },
        addAttributeFilter: function() {
            this.tb.settings.attributeFilter.push({
                name: this.templateAttributes[this.tb.settings.templateSet][0],
                logic: '>',
                value: 0.5
            })
        },
        requestRetro: function(smiles, callback) {
            showLoader()
            const url = '/api/v2/retro/';
            const body = {
                target: smiles,
                template_set: this.tb.settings.templateSet,
                template_prioritizer_version: this.tb.settings.templateSetVersion,
                precursor_prioritizer: this.tb.settings.precursorScoring,
                num_templates: this.tb.settings.numTemplates,
                max_cum_prob: this.tb.settings.maxCumProb,
                filter_threshold: this.tb.settings.minPlausibility,
                cluster_method: this.clusterOptions.cluster_method,
                cluster_feature: this.clusterOptions.feature,
                cluster_fp_type: this.clusterOptions.fingerprint,
                cluster_fp_length: this.clusterOptions.fpBits,
                cluster_fp_radius: this.clusterOptions.fpRadius,
                selec_check: this.tb.settings.allowSelec,
                attribute_filter: this.tb.settings.attributeFilter,
            };
            fetch(url,{
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                    'X-CSRFToken': getCookie('csrftoken')
                },
                body: JSON.stringify(body)
            })
                .then(resp => {
                    if (!resp.ok) {
                        throw Error(resp.statusText)
                    }
                    return resp
                })
                .then(resp => resp.json())
                .then(json => {
                    setTimeout(() => this.pollCeleryResult(json.task_id, callback), 1000)
                })
                .catch(error => {
                    hideLoader();
                    alert('There was an error predicting precursors for this target: '+error)
                })
        },
        pollCeleryResult: function(taskId, callback) {
            fetch(`/api/v2/celery/task/${taskId}/`)
            .then(resp => resp.json())
            .then(json => {
                if (json.complete) {
                    callback(json.output);
                    hideLoader();
                }
                else if (json.failed) {
                    hideLoader();
                    throw Error('Celery task failed.');
                }
                else {
                    setTimeout(() => {this.pollCeleryResult(taskId, callback)}, 1000)
                }
            })
            .catch(error => {
                if (error instanceof TypeError && error.message === 'Failed to fetch') {
                    console.log('Unable to fetch celery results due to connection error. Will keep trying.')
                    setTimeout(() => {this.pollCeleryResult(taskId, callback)}, 2000)
                } else {
                    console.error('There was a problem fetching results:', error);
                }
            });
        },
        resolveChemName: function(name) {
            if (this.enableResolve && this.allowResolve) {
                let url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'+encodeURIComponent(name)+'/property/IsomericSMILES/txt'
                return fetch(url)
                    .then(resp => {
                        if (resp.status == 404) {
                            throw Error(resp.statusText);
                        } else {
                            return resp.text()
                        }
                    })
                    .catch(err => {
                        throw Error('Cannot resolve "'+name+'" to smiles: '+err.message);
                    })
            } else {
                throw Error('Resolving chemical name using external server is not allowed.');
            }
        },
        validatesmiles: function(smiles, iswarning) {
            return fetch(
                '/api/v2/rdkit/smiles/validate/',
                {
                    method: 'POST',
                    headers: {
                        'Content-Type': 'application/json',
                        'X-CSRFToken': getCookie('csrftoken'),
                    },
                    body: JSON.stringify({
                        smiles: smiles
                    })
                }
            )
                .then(resp => {
                    if (!resp.ok) {
                        throw 'Unable to connect to server: error code '+resp.status;
                    } else {
                        return resp.json()
                    }
                })
                .then(res_json => {
                    if (!res_json['correct_syntax']) {
                        if(iswarning) alert('Input SMILES string: invalid syntax')
                        return false
                    } else if (!res_json['valid_chem_name']) {
                        if(iswarning) alert('Input SMILES string: invalid chemical')
                        return false
                    } else {
                        return true
                    }
                })
        },
        lookupPrice: function(smiles) {
            let url = `/api/v2/buyables/?q=${encodeURIComponent(smiles)}&canonicalize=True`
            if (!this.tb.settings.buyablesSourceAll) {
                if (this.tb.settings.buyablesSource.length) {
                    this.tb.settings.buyablesSource.forEach(source => {url += '&source=' + source})
                } else{
                    url += '&source=[]'
                }
            }
            return fetch(url)
                .then(resp => resp.json())
                .then(json => {
                    let result = {
                        'search': json.search,
                        'ppg': 'not buyable',
                        'source': '',
                    }
                    if (json.result.length > 0) {
                        result.ppg = json.result[0].ppg;
                        result.source = json.result[0].source;
                        for (entry of json.result) {
                            if (entry.ppg < result.ppg) {
                                result.ppg = entry.ppg;
                                result.source = entry.source;
                            }
                        }
                    }
                    return result
                })
        },
        changeTarget: function() {
            showLoader();
            this.saveTbSettings()
            this.saveNetworkOptions()
            this.saveIppSettings()
            this.validatesmiles(this.target, !this.allowResolve)
            .then(isvalidsmiles => {
                if (isvalidsmiles) {
                    return this.target
                } else {
                    return this.resolveChemName(this.target)
                }
            })
            .then(smiles => this.canonicalize(smiles, 'target'))
            .then(() => {
                this.saveTarget()
                if (this.target !== undefined) {
                    const smi = this.target;
                    const app = this;
                    function callback(precursors) {
                        app.dataGraph.clear()
                        app.dispGraph.clear()
                        app.initTargetDataNode()
                        app.initTargetDispNode()
                        app.initializeNetwork(app.dispGraph)
                        app.network.on('selectNode', app.showInfo)
                        app.network.on('deselectNode', app.clearSelection)
                        let addedReactions = app.addRetroResultToDataGraph(precursors, smi)
                        let reactionsToAdd = addedReactions.filter(reactionSmiles => {
                            return !app.allowCluster || app.dataGraph.nodes.get(reactionSmiles).clusterRep
                        }).slice(0, app.reactionLimit)
                        app.addRetroResultToDispGraph(reactionsToAdd, NIL_UUID)
                        app.network.selectNodes([NIL_UUID])
                        app.showInfo({'nodes': [NIL_UUID]})
                    }
                    this.requestRetro(smi, callback);
                } else {
                    hideLoader();
                }
            })
            .catch(error => {
                hideLoader();
                var error_msg = 'unknown error'
                if ('message' in error) {
                    error_msg = error.name+':'+error.message
                } else if (typeof(error) == 'string') {
                    error_msg = error
                }
                alert('There was an error fetching precursors for this target with the supplied settings: '+error_msg)
            })
        },
        initTargetDataNode() {
            this.dataGraph.nodes.add({
                id: this.target,
                type: 'chemical',
            })
            this.lookupPrice(this.target)
                .then(result => {
                    let node = this.dataGraph.nodes.get(this.target);
                    node.ppg = result.ppg;
                    node.source = result.source;
                })
        },
        initTargetDispNode() {
            this.dispGraph.nodes.add({
                id: NIL_UUID,
                smiles: this.target,
                borderWidth: 3,
                color: {
                    border: '#000088'
                },
                shape: 'image',
                image: this.getMolDrawEndPoint(this.target),
                type: 'chemical',
            })
        },
        addRetroResultToDataGraph(data, parentSmiles) {
            // Add results as reaction and chemical nodes under the specified parent chemical
            // Arguments should be list of outcome objects and the SMILES of the parent node
            let clusterTracker = new Set(this.clusteredResultsIndex[parentSmiles]);
            let addedReactions = [];
            for (let reaction of data) {
                let reactionSmiles = reaction['smiles'] + '>>' + parentSmiles
                addedReactions.push(reactionSmiles)
                let node = {
                    id: reactionSmiles,
                    rank: reaction['rank'],
                    ffScore: reaction['plausibility'],
                    retroScore: reaction['score'],
                    templateScore: reaction['template_score'],
                    templateRank: reaction['template_rank'],
                    templateIds: reaction['templates'],
                    clusterId: reaction['group_id'],
                    clusterRep: !clusterTracker.has(reaction['group_id']),  // Tag first result in each cluster as the representative
                    precursors: reaction['smiles_split'],
                    precursorSmiles: reaction['smiles'],
                    numExamples: reaction['num_examples'],
                    necessaryReagent: reaction['necessary_reagent'],
                    mappedSmiles: reaction['mapped_smiles'],
                    reactingAtoms: reaction['reacting_atoms'],
                    numRings: reaction['num_rings'],
                    rmsMolwt: reaction['rms_molwt'],
                    scscore: reaction['scscore'],
                    type: 'reaction',
                    inVis: {},  // Tracks associated nodes in display graph (key: parent ID, value: node ID)
                }

                clusterTracker.add(reaction['group_id'])
                if (node.templateIds) {
                    node.templateIds.forEach(template => this.apiTemplateCount(template))
                }

                if ('outcomes' in reaction) {
                    node['outcomes'] = reaction['outcomes'].split('.')
                    node['selectivity'] = new Array(node.outcomes.length)
                    node['mappedPrecursors'] = reaction['mapped_precursors']
                    node['mappedOutcomes'] = reaction['mapped_outcomes']
                } else if ('selec_error' in reaction) {
                    node['selecError'] = reaction['selec_error']
                }

                this.dataGraph.nodes.add(node)
                this.dataGraph.edges.add({
                    id: uuidv4(),
                    from: parentSmiles,
                    to: reactionSmiles,
                })

                for (let precursorSmiles of node.precursors) {
                    if (!this.dataGraph.nodes.get(precursorSmiles)) {
                        this.dataGraph.nodes.add({
                            id: precursorSmiles,
                            type: 'chemical',
                        })
                        this.lookupPrice(precursorSmiles)
                            .then(result => {
                                let node = this.dataGraph.nodes.get(precursorSmiles);
                                node.ppg = result.ppg;
                                node.source = result.source;
                            })
                    }
                    this.dataGraph.edges.add({
                        id: uuidv4(),
                        from: reactionSmiles,
                        to: precursorSmiles,
                    })
                }
            }
            this.recompute += 1  // Force recompute of properties
            return addedReactions
        },
        addRetroResultToDispGraph(data, parentId) {
            // Add reaction and chemical nodes with display properties under the specified parent chemical
            // Arguments should be list of reaction smiles to add and the ID of the parent node
            for (let reactionSmiles of data) {
                let reactionObj = this.dataGraph.nodes.get(reactionSmiles)
                let reactionId = uuidv4()
                reactionObj['inVis'][parentId] = reactionId

                let node = {
                    id: reactionId,
                    smiles: reactionSmiles,
                    label: '#' + reactionObj['rank'],
                    type: 'reaction',
                }

                if ('outcomes' in reactionObj) {
                    node['borderWidth'] = 2
                    node['color'] = { border: '#ff4444' }
                    node['title'] = "Selectivity warning! Select this node to see more details"
                } else if ('selecError' in reactionObj) {
                    node['borderWidth'] = 2
                    node['color'] = { border: '#ffbb00' }
                }

                this.dispGraph.nodes.add(node)
                this.dispGraph.edges.add({
                    id: uuidv4(),
                    from: parentId,
                    to: reactionId,
                    scaling: {
                        min: 1,
                        max: 5,
                        customScalingFunction: edgeScaling,
                    },
                    color: {
                        color: '#000000',
                        inherit: false,
                    },
                    value: reactionObj['templateScore'],
                })

                for (let precursorSmiles of reactionObj.precursors) {
                    let precursorObj = this.dataGraph.nodes.get(precursorSmiles)
                    let precursorId = uuidv4()
                    this.dispGraph.nodes.add({
                        id: precursorId,
                        smiles: precursorSmiles,
                        borderWidth: 2,
                        color: {
                            border: (precursorObj.ppg === 'not buyable') ? '#880000' : '#008800'
                        },
                        shape: 'image',
                        image: this.getMolDrawEndPoint(precursorSmiles),
                        type: 'chemical',
                    })
                    this.dispGraph.edges.add({
                        id: uuidv4(),
                        from: reactionId,
                        to: precursorId,
                        scaling: {
                            min: 1,
                            max: 5,
                            customScalingFunction: edgeScaling,
                        },
                        color: {
                            color: '#000000',
                            inherit: false,
                        },
                        value: reactionObj['templateScore'],
                    })
                }
            }
            this.recompute += 1  // Force recompute of properties
        },
        addTreeBuilderResultToDataGraph(data) {
            this.dataGraph.nodes.add(data.nodes.map(node => {
                if (node.type === 'chemical') {
                    return {
                        id: node['id'],
                        asReactant: node['as_reactant'],
                        asProduct: node['as_product'],
                        ppg: (node['purchase_price'] === 0) ? 'not buyable' : node['purchase_price'],
                        terminal: node['terminal'],
                        type: node['type'],
                    }
                } else {
                    node['tforms'].forEach(template => this.apiTemplateCount(template))
                    return {
                        id: node['id'],
                        rank: node['rank'],
                        ffScore: node['plausibility'],
                        retroScore: node['template_score'],
                        templateScore: node['template_score'],
                        templateIds: node['tforms'],
                        precursors: node['precursor_smiles'].split('.'),
                        precursorSmiles: node['precursor_smiles'],
                        numExamples: node['num_examples'],
                        necessaryReagent: node['necessary_reagent'],
                        numRings: node['num_rings'],
                        rmsMolwt: node['rms_molwt'],
                        scscore: node['scscore'],
                        type: node['type'],
                        inVis: {},
                    }
                }
            }))
            this.dataGraph.edges.add(data.links.map(edge => {
                return {
                    id: edge['id'],
                    from: edge['source'],
                    to: edge['target'],
                }
            }))
            this.recompute += 1  // Force recompute of properties
        },
        addTreeBuilderResultToDispGraph(data) {
            this.dispGraph.nodes.add(data.nodes.map(node => {
                let dataObj = this.dataGraph.nodes.get(node['smiles'])
                if (node.type === 'chemical') {
                    return {
                        id: node['id'],
                        smiles: node['smiles'],
                        borderWidth: (node['smiles'] === this.target) ? 3 : 2,
                        color: {
                            border: (node['smiles'] === this.target) ? '#000088' : (dataObj['ppg'] === 'not buyable') ? '#880000' : '#008800'
                        },
                        shape: 'image',
                        image: this.getMolDrawEndPoint(node['smiles']),
                        type: 'chemical',
                    }
                } else {
                    return {
                        id: node['id'],
                        smiles: node['smiles'],
                        label: '#' + node['rank'],
                        type: 'reaction',
                    }
                }
            }))
            this.dispGraph.edges.add(data.edges.map(edge => {
                let from = this.dispGraph.nodes.get(edge['from'])
                let to = this.dispGraph.nodes.get(edge['to'])
                let reactionObj
                if (from['type'] === 'reaction') {
                    reactionObj = this.dataGraph.nodes.get(from['smiles'])
                } else {
                    reactionObj = this.dataGraph.nodes.get(to['smiles'])
                    reactionObj.inVis[from['id']] = to['id']
                }
                return {
                    id: edge['id'],
                    from: edge['from'],
                    to: edge['to'],
                    scaling: {
                        min: 1,
                        max: 5,
                        customScalingFunction: edgeScaling,
                    },
                    color: {
                        color: '#000000',
                        inherit: false,
                    },
                    value: reactionObj['templateScore'],
                }
            }))
            this.recompute += 1  // Force recompute of properties
        },
        updateNetworkOptions() {
            if (typeof(this.network) != 'undefined') {
                this.network.setOptions(JSON.parse(JSON.stringify(this.networkOptions)))
            }
            
        },
        toggleHierarchical: function() {
            this.networkOptions.layout.hierarchical.enabled = !this.networkOptions.layout.hierarchical.enabled
            this.updateNetworkOptions()
        },
        expandNode: function() {
            if (this.isModalOpen() || typeof(this.network) == "undefined") {
                return
            }
            showLoader();
            var selected = this.network.getSelectedNodes();
            if (selected.length !== 1) {
              hideLoader();
              if (selected.length === 0) {
                  alert('Please select a terminal chemical node to expand')
              }
              else {
                  alert('Please only select 1 node at a time to expand')
              }
              return
            }
            var nodeId = selected[0];
            if (typeof(nodeId) == 'string' && nodeId.startsWith('cluster')) {
                alert('Cannot expand collpased node! To toggle collpased state, click collapse toggle button again with collapsed cluster selected.')
                hideLoader();
                return
            }
            var node = this.dispGraph.nodes.get(nodeId)
            if (node.type !== 'chemical') {
                alert('Cannot expand reaction; try expanding with a chemical node selected');
                hideLoader();
                return
            }
            var childrenOfSelected = this.dispGraph.getSuccessors(nodeId)
            if (childrenOfSelected.length !== 0) {
                alert("You've already expanded this node. If you would like to re-expand, please use the 'Remove children nodes' button to clear results in the visualization for this chemical. Please note that this will replace the previously predicted results for this chemical (for example, if you've changed any settings)")
                hideLoader();
                return
            }
            const smi = node.smiles;
            const app = this;
            function callback(precursors) {
                if (precursors.length === 0) {
                    alert('No precursors found!')
                } else {
                    let addedReactions = app.addRetroResultToDataGraph(precursors, smi)
                    let reactionsToAdd = addedReactions.filter(reactionSmiles => {
                        return !app.allowCluster || app.dataGraph.nodes.get(reactionSmiles).clusterRep
                    }).slice(0, app.reactionLimit)
                    app.addRetroResultToDispGraph(reactionsToAdd, nodeId)
                    app.showInfo({'nodes': [nodeId]})
                    app.network.fit()
                }
            }
            this.requestRetro(smi, callback);
        },
        apiTemplateCount: function(templateId) {
            if (typeof(this.templateNumExamples[templateId]) == 'undefined') {
                fetch('/api/template/?id='+templateId)
                .then(resp => resp.json())
                .then(json => {
                    var id = json["request"]["id"][0];
                    var count = json["template"]["count"];
                    this.templateNumExamples[id] = count;
                })
            }
        },
        deleteChoice: function() {
            // for all selected nodes, delete reaction nodes and delete children of chemical nodes
            let res = confirm('This will delete all selected reaction nodes and all children of the selected chemical nodes. Continue?')
            if (!res) {
                return
            }
            let selected = this.network.getSelectedNodes();
            for (let nodeId of selected) {
                let node = this.dispGraph.nodes.get(nodeId)
                if (node === null) {
                    // the node does not exist, it may have already been deleted
                    continue
                }
                if (node.type === 'chemical') {
                    this.deleteChildren(node)
                } else {
                    this.deleteNode(node)
                }
            }
            this.selected = null;
            this.network.unselectAll();
        },
        deleteNode: function(node) {
            // Delete the specified node and its children from the graph
            if (node.type !== 'reaction') {
                console.warn('Attempted to delete chemical node. Use deleteChildren instead.')
                return
            }
            // Delete children of children
            const childrenIds = this.dispGraph.getSuccessors(node.id)
            const children = this.dispGraph.nodes.get(childrenIds)
            children.forEach(this.deleteChildren)
            // Delete children
            this.dispGraph.nodes.remove(childrenIds)
            // Update inVis property of corresponding node in dataGraph
            const dataObj = this.dataGraph.nodes.get(node.smiles)
            const parentNodeId = this.dispGraph.getPredecessors(node.id)[0]
            delete dataObj.inVis[parentNodeId]
            // Delete selected node
            this.dispGraph.nodes.remove(node.id)
            this.dispGraph.trimDanglingEdges()
        },
        deleteChildren: function(node) {
            // Delete children of the specified node from the graph
            if (node.type !== 'chemical') {
                console.warn('Attempted to delete children of reaction node. Use deleteNode instead.')
                return
            }
            const childrenIds = this.dispGraph.getSuccessors(node.id)
            const children = this.dispGraph.nodes.get(childrenIds)
            children.forEach(this.deleteNode)
            this.dispGraph.trimDanglingEdges()
        },
        toggleResolver: function() {
            if (this.allowResolve) {
                this.allowResolve = false
            }
            else {
                this.allowResolve = true
            }
        },
        download: function() {
            let data = {
                dataGraph: this.dataGraph,
                dispGraph: this.dispGraph,
                version: 1.0,
            }
            let dataStr = "data:text/json;charset=utf-8," + encodeURIComponent(JSON.stringify(data));
            let dlAnchorElem = document.getElementById('downloadAnchorElem');
            dlAnchorElem.setAttribute("href", dataStr);
            dlAnchorElem.setAttribute("download", this.downloadName);
            dlAnchorElem.click();
        },
        hasUndefinedClusterId: function() {
            for (let rxn of this.dataGraph.nodes.get({filter: node => node.type === 'reaction'})) {
                if (rxn.clusterId === undefined) {
                    return true
                }
            }
            return false
        },
        load: function() {
            if (this.clear()) {
                let file = document.getElementById("loadNetwork").files[0];
                let reader = new FileReader();
                reader.readAsText(file)
                reader.onload = function(e) {
                    let data = JSON.parse(e.target.result)
                    if (data.version === 1.0) {
                        app.importDataV1(data)
                    } else {
                        app.importDataV0(data)
                    }
                }
            }
        },
        importDataV0: function(data) {
            // Parse old data download format, before version numbers were introduced
            // Top level properties should be nodes, edges, and results
            // This converts the data to the current format for subsequent downloads
            this.target = data.nodes[0].smiles
            this.initTargetDataNode()
            for (let [chem, precursors] of Object.entries(data.results)) {
                this.addRetroResultToDataGraph(precursors, chem)
            }
            this.dispGraph.nodes.add(data.nodes.map(node => {
                if (node.type === 'chemical') {
                    let dataObj = this.dataGraph.nodes.get(node['smiles'])
                    // Transfer properties which were not in the results object
                    dataObj.ppg = node['ppg']
                    dataObj.source = node['source']
                    return {
                        id: node['id'],
                        smiles: node['smiles'],
                        borderWidth: node['borderWidth'],
                        color: node['color'],
                        shape: 'image',
                        image: this.getMolDrawEndPoint(node['smiles']),
                        type: 'chemical',
                    }
                } else {
                    let dataObj = this.dataGraph.nodes.get(node['reactionSmiles'])
                    // Transfer properties which were not in the results object
                    dataObj.selectivity = node['selectivity']
                    let newNode = {
                        id: node['id'],
                        smiles: node['reactionSmiles'],
                        label: node['label'],
                        type: 'reaction',
                    }
                    for (let key of ['borderWidth', 'color', 'title']) {
                        if (key in node) {
                            newNode[key] = node[key]
                        }
                    }
                    return newNode
                }
            }))
            this.dispGraph.edges.add(data.edges.map(edge => {
                let from = this.dispGraph.nodes.get(edge['from'])
                let to = this.dispGraph.nodes.get(edge['to'])
                let reactionObj
                if (from['type'] === 'reaction') {
                    reactionObj = this.dataGraph.nodes.get(from['smiles'])
                } else {
                    reactionObj = this.dataGraph.nodes.get(to['smiles'])
                    reactionObj.inVis[from['id']] = to['id']
                }
                return {
                    id: edge['id'],
                    from: edge['from'],
                    to: edge['to'],
                    scaling: {
                        min: 1,
                        max: 5,
                        customScalingFunction: edgeScaling,
                    },
                    color: edge['color'],
                    value: edge['value'],
                }
            }))
            if (this.hasUndefinedClusterId()) {
                if (confirm('The uploaded json file does not have reaction cluster information for some precursors. Select "Ok" to re-cluster all of them. Select "Cancel" to disable clustering.')) {
                    this.updateAllClusters()
                } else {
                    this.allowCluster = false;
                }
            }
            this.initializeNetwork(this.dispGraph)
            this.network.on('selectNode', this.showInfo)
            this.network.on('deselectNode', this.clearSelection)
        },
        importDataV1: function(data) {
            // Parse data format version 1.0
            // Top level properties should be dataGraph, dispGraph, and version
            this.dataGraph.nodes.add(data.dataGraph.nodes)
            this.dataGraph.edges.add(data.dataGraph.edges)
            this.dispGraph.nodes.add(data.dispGraph.nodes)
            this.dispGraph.edges.add(data.dispGraph.edges)
            this.target = this.dispGraph.nodes.get(NIL_UUID).smiles
            if (this.hasUndefinedClusterId()) {
                if (confirm('The uploaded json file does not have reaction cluster information for some precursors. Select "Ok" to re-cluster all of them. Select "Cancel" to disable clustering.')) {
                    this.updateAllClusters()
                } else {
                    this.allowCluster = false;
                }
            }
            this.initializeNetwork(this.dispGraph)
            this.network.on('selectNode', this.showInfo)
            this.network.on('deselectNode', this.clearSelection)
        },
        updateAllClusters: function () {
            // Recompute cluster IDs for all results
            for (let target of this.dataGraph.nodes.getIds({filter: node => node.type === 'chemical'})) {
                if (this.dataGraph.getSuccessors(target).length > 0) {
                    this.requestClusterId(target)
                }
            }
        },
        clear: function(skipConfirm = false) {
            // Returns true or false depending on whether results were cleared
            if (skipConfirm || confirm('This will clear all of your current results. Continue anyway?')) {
                this.target = '';
                this.selected = null;
                this.dataGraph.clear();
                this.dispGraph.clear();
                return true
            } else {
                return false
            }
        },
        clearSelection: function() {
            this.selected = null;
        },
        copySelectedSmiles: function() {
            let copyTooltip = document.querySelector('#copy-tooltip')
            copyToClipboard(this.selected.smiles)
            copyTooltip.innerHTML = 'Copied!'
            setTimeout(() => {copyTooltip.innerHTML = "Click to copy!"}, 2000)
        },
        collapseNode: function() {
            let selected = this.network.getSelectedNodes();
            selected.forEach(node => {
                if (this.network.clustering.isCluster(node)) {
                    this.network.openCluster(node)
                } else {
                    let forCluster = this.dispGraph.getAllSuccessors(node);
                    let options = {
                        joinCondition: (nodeOptions) => {
                            return forCluster.includes(nodeOptions.id) || nodeOptions.id === node
                        }
                    };
                    this.network.clustering.cluster(options);
                }
            })
        },
        addFromResults: function(selected, reaction) {
            if (selected.id in reaction.inVis) {
                return
            }
            this.addRetroResultToDispGraph([reaction.id], selected.id)
            document.querySelector('.addRes[data-rank="'+Number(reaction.rank)+'"]').style.display='none';
            document.querySelector('.remRes[data-rank="'+Number(reaction.rank)+'"]').style.display='';
        },
        remFromResults: function(selected, reaction) {
            let node = this.dispGraph.nodes.get(reaction.inVis[selected.id])
            this.deleteNode(node)
            document.querySelector('.addRes[data-rank="'+Number(reaction.rank)+'"]').style.display='';
            document.querySelector('.remRes[data-rank="'+Number(reaction.rank)+'"]').style.display='none';
        },
        resetSortingCategory: function() {
            this.sortingCategory = 'retroScore'
            this.selectSortingOrder()
        },
        selectSortingOrder: function() {
            this.sortOrderAscending = ["rmsMolwt", "numRings", "scscore", "templateRank"].includes(this.sortingCategory) || (this.sortingCategory === 'retroScore' && this.tb.settings.precursorScoring === 'SCScore');
        },
        applyFilterReactingAtoms: function() {
            // Not in use right now, but could come in handy later.
            return 
        },
        checkFilter: function(result)
        {
            // If Ketcher is not active, there is no way for user to interact with filter;
            // the filter should pass
            if (!$('#ketcher-iframe-min')[0]) {
                return true;
            }

            var reactingAtoms = result.reactingAtoms.map((el) => el-1);
            var reactingAtomsArray = Array.from(reactingAtoms.values());

            var selection = $('#ketcher-iframe-min')[0].contentWindow.ketcher.editor.selection();
            if (selection && selection.atoms) {
                var selectionAtomsArray = Array.from(selection.atoms.values());
            }
            else {
                // No atoms are selected in Ketcher
                var selectionAtomsArray = []
            }

            var filterResult = reactingAtomsArray.some(element => selectionAtomsArray.includes(element));
            return filterResult;
        },
        showInfo: function(obj) {
            let nodeId = obj.nodes[obj.nodes.length-1];
            let dispObj = this.dispGraph.nodes.get(nodeId);
            if (!dispObj) return
            let dataObj = this.dataGraph.nodes.get(dispObj.smiles);
            if (!dataObj) return

            this.selected = {
                'id': dispObj.id,
                'smiles': dispObj.smiles,
                'type': dispObj.type,
                'data': dataObj,
                'disp': dispObj,
            };

            if (dispObj.type === 'chemical') {
                setSmilesDrawingKetcherMin(dataObj.id);
                if (!dataObj.source) {
                    this.lookupPrice(dataObj.id)
                        .then(result => {
                            dataObj.ppg = result.ppg;
                            dataObj.source = result.source;
                        })
                }
            }
        },
        openModal: function(modalName) {
            /*
            this.clearSelection();
            if (network) {
                network.unselectAll();
            }
            */
            if (modalName == "settings") {
                this.showSettingsModal = true
            }
            else if (modalName == "download") {
                this.showDownloadModal = true
            }
            else if (modalName == "load") {
                this.showLoadModal = true
            }
        },
        openClusterPopoutModal: function(selected, res) {
            if(selected === undefined) {
                alert('No target molecule selected. Please select a molecule in the tree.')
                return
            }
            this.$set(this.clusterPopoutModalData, 'selected', selected);
            this.$set(this.clusterPopoutModalData, 'selectedSmiles', selected.smiles);
            this.$set(this.clusterPopoutModalData, 'res', res);
            this.$set(this.clusterPopoutModalData, 'clusterId', res.clusterId);
            this.showClusterPopoutModal = true;
        },
        closeClusterPopoutModal: function() {
            this.showClusterPopoutModal = false;
            this.clusterPopoutModalData['selected'] = undefined;
            this.clusterPopoutModalData['selectedSmiles'] = undefined;
            this.clusterPopoutModalData['res'] = undefined;
            this.clusterPopoutModalData['clusterId'] = undefined;
        },
        clusterPopoutModalIncGroupID: function() {
            let allIds = this.clusteredResultsIndex[this.clusterPopoutModalData['selectedSmiles']];
            let idx = allIds.indexOf(this.clusterPopoutModalData['clusterId']);
            if (idx !== allIds.length-1) {
                this.$set(this.clusterPopoutModalData, 'clusterId', allIds[idx+1]);
            }
        },
        clusterPopoutModalDecGroupID: function() {
            let allIds = this.clusteredResultsIndex[this.clusterPopoutModalData['selectedSmiles']];
            let idx = allIds.indexOf(this.clusterPopoutModalData['clusterId']);
            if (idx !== 0) {
                this.$set(this.clusterPopoutModalData, 'clusterId', allIds[idx-1]);
            }
        },
        openClusterEditModal: function(selected, clusterId) {
            if(selected === undefined) {
                alert('No target molecule selected. Please select a molecule in the tree.')
                return
            }
            if(clusterId === undefined) {
                clusterId = 0
            }
            this.$set(this.clusterEditModalData, 'selected', selected);
            this.$set(this.clusterEditModalData, 'selectedSmiles', selected.smiles);
            this.$set(this.clusterEditModalData, 'clusterId', clusterId);
            this.showClusterEditModal = true;
        },
        closeClusterEditModal: function() {
            this.showClusterEditModal = false;
            this.clusterEditModalData['selected'] = undefined;
            this.clusterEditModalData['selectedSmiles'] = undefined;
            this.clusterEditModalData['clusterId'] = undefined;
        },
        clusteredit_dragstart_handler: function(precursor, event) {
            event.target.style.opacity = '0.4';
            event.dataTransfer.setData('text/plain', precursor.id);
            var img = new Image();
            img.src = this.getMolDrawEndPoint(precursor.precursorSmiles);
            // set opacity does not work..
            event.dataTransfer.setDragImage(img, 10, 10);
            event.dataTransfer.effectAllowed = 'all';
            // disable all buttons on dragging
            var buttons = document.querySelectorAll("button");
            buttons.forEach(function(e){e.style.pointerEvents = "none";});
        },
        clusteredit_dragend_handler: function(event) {
            event.target.style.opacity = '1';
            // enable all buttons
            var buttons = document.querySelectorAll("button");
            buttons.forEach(function(e){e.style.pointerEvents = "all";});
        },
        clusteredit_drop_handler: function(target, event) {
            event.preventDefault(); // important
            let smi = event.dataTransfer.getData('text/plain')  // precursor.id
            let obj = this.dataGraph.nodes.get(smi)
            let oldId = obj.clusterId
            obj.clusterId = target.clusterId
            this.clusteredit_dragend_handler(event);
            clusteredit_dragleave_handler(event);
            this.updateClusterReps(this.clusterEditModalData['selectedSmiles'], [oldId, obj.clusterId])
        },
        clusteredit_drop_handler_newcluster: function(event) {
            event.preventDefault(); // important
            let smi = event.dataTransfer.getData('text/plain');  // precursor.id
            let obj = this.dataGraph.nodes.get(smi)
            let allIds = this.clusteredResultsIndex[this.clusterEditModalData['selectedSmiles']]
            let oldId = obj.clusterId
            obj.clusterId = allIds[allIds.length - 1] + 1;
            this.clusteredit_dragend_handler(event);
            clusteredit_dragleave_handler(event);
            this.updateClusterReps(this.clusterEditModalData['selectedSmiles'], [oldId, obj.clusterId])
        },
        updateClusterReps(target, clusterIds) {
            // Update the cluster representatives for the specified clusterIds
            for (let clusterId of clusterIds) {
                let options = {filter: item => item.clusterId === clusterId}
                let precursorSmiles = this.dataGraph.getSuccessors(target)
                let precursors = this.dataGraph.nodes.get(precursorSmiles, options)
                precursors.forEach((item, index) => item.clusterRep = index === 0)  // Set first item as cluster rep
            }
            this.recompute += 1
        },
        clusterEditModalIncGroupID: function() {
            let allIds = this.clusteredResultsIndex[this.clusterEditModalData['selectedSmiles']];
            let idx = allIds.indexOf(this.clusterEditModalData['clusterId']);
            if (idx !== allIds.length-1) {
                this.$set(this.clusterEditModalData, 'clusterId', allIds[idx+1]);
            }
        },
        clusterEditModalDecGroupID: function() {
            let allIds = this.clusteredResultsIndex[this.clusterEditModalData['selectedSmiles']];
            let idx = allIds.indexOf(this.clusterEditModalData['clusterId']);
            if (idx !== 0) {
                this.$set(this.clusterEditModalData, 'clusterId', allIds[idx-1]);
            }
        },
        clusterEditModalAddPrecursor: function(selectedSmiles, smiles, clusterId) {
            // clusterId == undefined is to add a new cluster
            let successors = this.dataGraph.getSuccessors(selectedSmiles)
            let allIds = this.clusteredResultsIndex[selectedSmiles]
            if (clusterId === undefined) {
                if (successors.length === 0 || allIds.length === 0) {
                    clusterId = 0
                } else {
                    clusterId = allIds[allIds.length - 1] + 1
                }
            }
            let rank = Math.max(...this.dataGraph.nodes.get(successors).map(s => s.rank)) + 1;
            let res = {
                smiles: smiles,
                smiles_split: smiles.split('.'),
                rank: rank,
                group_id: clusterId,
            }
            this.addRetroResultToDataGraph([res], selectedSmiles)
        },
        openAddNewPrecursorModal: function(selectedSmiles, clusterId) {
            // if clusterId == undefined, add to a new group
            this.showAddNewPrecursorModal = true;
            this.$set(this.addNewPrecursorModal, 'selectedSmiles', selectedSmiles === undefined ? this.selected.smiles : selectedSmiles);
            this.$set(this.addNewPrecursorModal, 'clusterId', clusterId);
            this.$set(this.addNewPrecursorModal, 'newPrecursorSmiles', '');
            this.$set(this.addNewPrecursorModal, 'noDupCheck', false);
        },
        closeAddNewPrecursorModal: function() {
            this.showAddNewPrecursorModal = false;
            this.addNewPrecursorModal['selectedSmiles'] = '';
            this.addNewPrecursorModal['clusterId'] = '';
            this.addNewPrecursorModal['newPrecursorSmiles'] = '';
            this.addNewPrecursorModal['noDupCheck'] = false;
        },
        checkDuplicatePrecursor: function(selectedSmiles, p) {
            let existing = this.dataGraph.getSuccessors(selectedSmiles)
            return (existing.includes(p)) ? this.dataGraph.nodes.get(p) : undefined
        },
        addNewPrecursorModalSubmit: function() {
            let gid = this.addNewPrecursorModal['clusterId']
            let smi = this.addNewPrecursorModal['newPrecursorSmiles']
            this.validatesmiles(smi, !this.allowResolve)
                .then(isValid => {
                    return isValid ? smi : this.resolveChemName(smi)
                })
                .then(smi => this.canonicalize(smi, (res) => this.addNewPrecursorModal['newPrecursorSmiles'] = res))
                .then(() => {
                    if (this.addNewPrecursorModal['newPrecursorSmiles'] !== undefined) {
                        if (!this.addNewPrecursorModal['noDupCheck']) {
                            let s = this.checkDuplicatePrecursor(this.addNewPrecursorModal['selectedSmiles'], this.addNewPrecursorModal['newPrecursorSmiles'])
                            if (s !== undefined) {
                                alert('There may be a duplicated precursor: rank: '+s.rank+' cluster: '+s.clusterId+'. If you still want to proceed, please select "No duplicate check" option.')
                                return
                            }
                        }
                        this.clusterEditModalAddPrecursor(this.addNewPrecursorModal['selectedSmiles'], this.addNewPrecursorModal['newPrecursorSmiles'], gid)
                        this.closeAddNewPrecursorModal();
                    }
                })
                .catch(error => {
                    var error_msg = 'unknown error'
                    if ('message' in error) {
                        error_msg = error.name+':'+error.message
                    } else if (typeof(error) == 'string') {
                        error_msg = error
                    }
                    alert('There was an error fetching precursors for this target with the supplied settings: '+error_msg)
                })
        },
        getMolDrawEndPoint: function(precursor, isHighlight, isTransparent) {
            //  precursor can be
            //      1) a smiles string,
            //      2) a object with properties "reactingAtoms" and "mappedSmiles"
            //      3) a object with property "smiles"
            //      4) an object with property "precursorSmiles"
            //  isTransparent is false by default
            //  isHighlight is set to this.isHighlight by default, but can be overidden
            if (isHighlight === undefined) {
                isHighlight = this.isHighlightAtom;
            }
            if (isTransparent === undefined) {
                isTransparent = false;
            }
            let smiles;
            let mappedSmiles;
            let reactingAtoms;
            let url;
            if (typeof(precursor) == "string") {
                smiles = precursor;
                isHighlight = false;
            } else if (typeof(precursor) == "object") {
                if (precursor.mappedSmiles !== undefined && precursor.reactingAtoms !== undefined) {
                    mappedSmiles = precursor.mappedSmiles;
                    reactingAtoms = precursor.reactingAtoms;
                }
                if (precursor.smiles !== undefined) {
                    smiles = precursor.smiles;
                } else if (precursor.precursorSmiles !== undefined) {
                    smiles = precursor.precursorSmiles
                }
            }
            if (isHighlight && mappedSmiles !== undefined && reactingAtoms !== undefined) {
                url = `/api/v2/draw/?smiles=${encodeURIComponent(mappedSmiles)}&highlight=true`
                for (let ra of reactingAtoms) {
                    url += `&reacting_atoms=${ra}`
                }
            } else {
                if (smiles === undefined) {
                    console.log('Error: cannot plot precursor=' + precursor)
                    return ''
                }
                url = `/api/v2/draw/?smiles=${encodeURIComponent(smiles)}`
            }
            if (isTransparent) {
                url += '&transparent=true';
            }
            return url;
        },
        isModalOpen: function() {
            var res = false;
            res = res || this.showSettingsModal;
            res = res || this.showDownloadModal;
            res = res || this.showLoadModal;
            res = res || this.showClusterPopoutModal;
            res = res || this.showClusterEditModal;
            res = res || this.showAddNewPrecursorModal;
            return res;
        },
        startTour: function() {
            if (confirm('Starting the tutorial will clear all of your current results. Continue anyway?')) {
                this.clear(true);
                tour.restart();
            }
        },
        requestClusterId: function(smiles) {
            showLoader();
            let rxns = this.dataGraph.nodes.get(this.dataGraph.getSuccessors(smiles))
            let outcomes = rxns.map(rxn => rxn.precursorSmiles)
            let scores = rxns.map(rxn => rxn.score || 0)
            let url = '/api/v2/cluster/'
            let body = {
                original: smiles,
                outcomes: outcomes,
                feature: this.clusterOptions.feature,
                fp_name: this.clusterOptions.fingerprint,
                fpradius: this.clusterOptions.fpRadius,
                fpnbits: this.clusterOptions.fpBits,
                cluster_method: this.clusterOptions.cluster_method,
                scores: scores,
            }

            fetch(url, {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                    'X-CSRFToken': getCookie('csrftoken')
                },
                body: JSON.stringify(body),
            })
                .then(resp => {
                    if (!resp.ok) {
                        throw resp.statusText;
                    }
                    return resp
                })
                .then(resp => resp.json())
                .then(json => {
                    let clusterTracker = new Set()
                    rxns.forEach((rxn, i) => {
                        let cid = json.output[i]
                        rxn.clusterId = cid
                        rxn.clusterRep = !clusterTracker.has(cid)
                        clusterTracker.add(cid)
                    })
                    hideLoader()
                })
                .catch(error => {
                    hideLoader()
                    alert('There was an error fetching cluster results for this target with the supplied settings: '+error)
                })
        },
        predictSelectivity: function(){
            showLoader();
            let data = this.selected.data
            let url = '/api/v2/general-selectivity/';
            let body = {
                reactants: data.mappedPrecursors,
                product: data.mappedOutcomes,
                mapped: true,
                all_outcomes: true,
                verbose: false,
                mode: this.selectivityModel,
            }
            fetch(url, {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                    'X-CSRFToken': getCookie('csrftoken')
                },
                body: JSON.stringify(body)
            })
                .then(resp => {
                    if (!resp.ok) {
                        throw Error(resp.statusText)
                    }
                    return resp
                })
                .then(resp => resp.json())
                .then(json => {
                    function callback(result) {
                        if (Array.isArray(result)) {
                            data.selectivity = result;
                        } else {
                            alert('Could not predict selectivity for this reaction.')
                        }
                    }
                    setTimeout(() => this.pollCeleryResult(json.task_id, callback), 1000)
                })
                .catch(error => {
                    hideLoader();
                    alert('There was an error predicting selectivity for this reaction: '+error)
                })
        },
        loadFromTreeBuilder: function(objectId, numTrees) {
            this.allowCluster = false
            showLoader()
            let url = `/api/v2/results/${objectId}/ipp/`
            if (numTrees !== 'all') {
                url += `?num=${numTrees}`
            }
            fetch(url)
                .then(resp => {
                    if (!resp.ok) {
                        throw Error(resp.statusText)
                    }
                    return resp.json()
                })
                .then(json => {
                    if (json.error) {
                        alert(json.error)
                    }
                    let result = json['result'];
                    this.target = result['settings']['smiles'];
                    let graph = result['result']['graph']
                    let tree = result['result']['tree']
                    this.addTreeBuilderResultToDataGraph(graph)
                    if (!!tree) {
                        this.addTreeBuilderResultToDispGraph(tree)
                    } else {
                        this.initTargetDispNode()
                    }
                    this.networkOptions.layout.hierarchical.enabled = true
                    this.initializeNetwork(this.dispGraph);
                    this.network.on('selectNode', this.showInfo);
                    this.network.on('deselectNode', this.clearSelection);
                    this.network.on('afterDrawing', hideLoader);
                })
                .catch(error => {
                    hideLoader()
                    console.error('There was a problem loading tree builder results:', error);
                });
        },
        canonicalize(smiles, input) {
            return fetch(
                '/api/v2/rdkit/smiles/canonicalize/',
                {
                    method: 'POST',
                    headers: {
                        'Content-Type': 'application/json',
                        'X-CSRFToken': getCookie('csrftoken'),
                    },
                    body: JSON.stringify({
                        smiles: smiles
                    })
                }
            )
            .then(resp => resp.json())
            .then(json => {
                if (json.smiles) {
                    if (typeof input === 'string') {
                        this[input] = json.smiles
                    } else if (input instanceof Function) {
                        input(json.smiles)
                    }
                }
            })
        },
        updateSmilesFromKetcher() {
            let smiles = ketcher.getSmiles();
            this.canonicalize(smiles, drawBoxId)
        },
        resetTemplateSetVersion(event) {
            this.tb.settings.attributeFilter = []
            this.tb.settings.templateSetVersion = this.templateSets[event.target.value][0]
        }
    },
    computed: {
        clusteredResults: function() {
            // Mapping of target smiles to list of cluster representatives
            // {'targetSmiles0':[{cluster0_rep}, {cluster1_rep},...], ...}
            let recompute = this.recompute
            let res = {}
            let options = {
                filter: (item) => item.clusterRep,
                order: (a, b) => a.clusterId - b.clusterId,
            }
            for (let target of this.dataGraph.nodes.getIds({filter: node => node.type === 'chemical'})) {
                res[target] = this.dataGraph.nodes.get(this.dataGraph.getSuccessors(target), options)
            }
            return res
        },
        clusteredResultsIndex: function() {
            // Mapping of target smiles to list of cluster IDs in ascending order
            // {'targetSmiles0':[unique cluster Ids in ascending order], ...}
            let recompute = this.recompute
            let res = {}
            for (let target of this.dataGraph.nodes.getIds({filter: node => node.type === 'chemical'})) {
                let precursors = this.dataGraph.nodes.get(this.dataGraph.getSuccessors(target))
                res[target] = [...new Set(precursors.map(rxn => rxn.clusterId))]
                res[target].sort((a, b) => a - b)
            }
            return res
        },
        currentClusterViewPrecursors: function() {
            let recompute = this.recompute
            let smi = this.clusterPopoutModalData.selectedSmiles
            let cid = this.clusterPopoutModalData.clusterId
            let precursorSmiles = this.dataGraph.getSuccessors(smi)
            let options = {
                filter: (item) => item.clusterId === cid
            }
            return this.dataGraph.nodes.get(precursorSmiles, options)
        },
        currentClusterEditPrecursors: function () {
            let recompute = this.recompute
            let smi = this.clusterEditModalData.selectedSmiles
            let cid = this.clusterEditModalData.clusterId
            let precursorSmiles = this.dataGraph.getSuccessors(smi)
            let options = {
                filter: (item) => item.clusterId === cid
            }
            return this.dataGraph.nodes.get(precursorSmiles, options)
        },
        resultsAvailable: function() {
            // Boolean of whether the selected chemical has been expanded
            // Returns false if a reaction node is selected
            let recompute = this.recompute
            if (this.selected.type !== 'chemical') {
                return false
            } else {
                return this.dataGraph.getSuccessors(this.selected.smiles).length !== 0
            }
        },
        currentPrecursors: function() {
            // Array of precursors corresponding to the selected chemical
            let recompute = this.recompute
            let refreshFilter = this.refreshFilter
            if (this.selected.type !== 'chemical') {
                return []
            }
            let precursorSmiles = this.dataGraph.getSuccessors(this.selected.smiles)
            let cmp = (this.sortOrderAscending) ? (a, b) => {return a - b} : (a, b) => {return b - a};
            let options = {
                filter: (item) => {
                    if (this.allowCluster && !item.clusterRep) {
                        return false
                    }
                    if (this.filterReactingAtoms && !this.checkFilter(item)) {
                        return false
                    }
                    return true
                },
                order: (a, b) => {
                    let a_ = a[this.sortingCategory] === undefined ? 0 : a[this.sortingCategory];
                    let b_ = b[this.sortingCategory] === undefined ? 0 : b[this.sortingCategory];
                    if (a_ === b_) {
                        return a.rank - b.rank
                    }
                    return cmp(a_, b_);
                }
            }
            return this.dataGraph.nodes.get(precursorSmiles, options)
        }
    },
    delimiters: ['%%', '%%'],
});

var tour = new Tour({
    framework: 'bootstrap4',
    storage: false,
    steps: [
        {
            title: "A guided tour through retrosynthesis",
            content: "Welcome to this guided tour through retrosynthesis planning using our interactive path planning tool. This will demonstrate the purpose of the tool and explain the user interface using a real example. Thanks to <a href='https://github.com/IGreatlyDislikeJavascript/bootstrap-tourist' target='_blank'>bootstrap-tourist</a> for the great guided tour JavaScript package making it very easy to provide this tour to you!",
            orphan: true,
            backdropContainer: '#body'
        },
        {
            element: "#target",
            title: "Start with a target compound",
            content: "You start the retrosynthetic planning by entering a target compounds SMILES formatted string here. If the name resolver is enabled (server icon to the left is green) you can also enter a chemical name. The name will be resolved using a third-party server (PubChem). If you wish to turn the name resolving feature off, click the server icon and it will turn red. For this tutorial we're going to explore an example reaction for <a href='https://en.wikipedia.org/wiki/Fluconazole' target='_blank'>Fluconazole</a>. Press 'Next' to continue!",
            placement: "bottom",
            onNext: function() {
                app.target = 'OC(Cn1cncn1)(Cn2cncn2)c3ccc(F)cc3F'
            }
        },
        {
            element: "#target",
            title: "Fluconazole",
            content: "Here's the SMILES string for Fluconazole. If you're unfamiliar with the SMILES format, click on the edit icon to open the drawing tool or try using software like ChemDraw to draw a structure and copy it's SMILES string (right click -> molecule -> copy as -> SMILES). Click 'Next to continue!",
            placement: "bottom",
            onNext: function() {
                if (app.data.nodes.length == null | app.data.nodes.length == 0) {
                    app.changeTarget();
                }
            }
        },
        {
            element: "#network",
            title: "One-step retrosynthesis results",
            content: "When the results are ready, they will be shown in the main window on the left. The target molecule you entered will be shown in the center inside a <span class='blue-text'>blue</span> box (it is currently selected). You can click and drag on empty space in the window to navigate through the entire network when zoomed in. Scrolling inside the window will zoom in and out. You can rearrange nodes by clicking and dragging on them. Take a second to enjoy the inverted gravity model, courtesy of <a href='http://visjs.org' target='_blank'>vis.js</a>. Click 'Next' to continue.",
            placement: 'right',
            backdropContainer: '#network'
        },
        {
            element: "#network",
            title: "Predicted reactions",
            content: "The children node(s) of your target molecule represent predicted <b>reactions</b> that may result in your target molecule. The number inside this node is the rank of the precursor, scored by the precursor prioritization method currently selected (more on this later). On the left you can see that the highest ranked prediction is highlighted.",
            onShown: function () {
                app.data.nodes.forEach(function(n) {
                    if (n.reactionSmiles === 'Fc1ccc(C2(Cn3cncn3)CO2)c(F)c1.c1nc[nH]n1>>OC(Cn1cncn1)(Cn2cncn2)c3ccc(F)cc3F') {
                        app.network.selectNodes([n.id])
                        app.selected = app.data.nodes.get(n.id);
                    }
                })
            },
            placement: 'right',
        },
        {
            element: "#network",
            title: "Reactants",
            content: "The children node(s) of <b>reactions</b> represent <b>chemicals</b>, and are the predicted reactants for this reaction. Chemicals in a <span class='red-text'>red</span> box were not found in the buyables database. <b>Chemicals</b> in a <span class='green-text'>green</span> box were found in the database and are buyable.",
            placement: 'right',
            onNext: function() {
                app.data.nodes.forEach(function(n) {
                    if (n.smiles === 'Fc1ccc(C2(Cn3cncn3)CO2)c(F)c1') {
                        app.network.selectNodes([n.id])
                        app.selected = app.data.nodes.get(n.id);
                    }
                })
            }
        },
        {
            element: "#network",
            placement: 'right',
            title: "Reactants",
            content: "In this example, we'll see if we can predict a reaction to make the non-buyable reactant in the rank 1 reaction prediction from buyable starting materials (it's been selected for you)."
        },
        {
            element: '#expand-btn',
            title: "Expanding chemical nodes",
            content: "To make a new prediction for the non-buyable chemical, which been highlighted in blue below, you would normally click the <b>Expand Node</b> button. As you are in the tutorial, please click 'Next' to see what would happen when you click this button.",
            placement: "bottom",
            reflex: true,
            onNext: function() {
                app.expandNode();
            }
        },
        {
            element: '#network',
            title: "Expanding chemical nodes",
            content: "A new prediction was made for this non-buyable chemical, and when everything is ready, the results will be added to the network visualization. It might look hectic at first, but the appropriate node positions should resolve quickly (thanks again to the <a href='http://visjs.org' target='_blank'>vis.js</a> inverted gravity!). If not, click and drag a node to give it a jiggle.",
            placement: "right"
        },
        {
            element: '#details',
            title: "Result details",
            content: "You may have noticed there's been a lot going on on the right side of the screen in addition to the changes in the network visualization. On this side, details of the currently selected node are shown. In this case, a <b>chemical</b> node is selected. At the top you can see its SMILES string, its cost in $/g (if buyable) and a 2d rendering of its structure.",
            placement: "left"
        },
        {
            element: '#details',
            title: "Precursors",
            content: "Additionally, if you've already made a retrosynthetic prediction for the currently selected <b>chemical</b>, you'll see list of the precursor results. Each entry shows the reactants for the reaction to make the currently selected chemical. Additional information such as a relative score and the number of examples there were for the templates that support the suggested reaction are also shown. If you haven't performed a retrosynthetic prediction for the selected chemical, the same <b>Expand Node</b> button you used before will be shown.",
            placement: "left"
        },
        {
            element: '#details',
            title: "Adding and removing reactions",
            content: "You may also notice there are many more precursor results shown on the right side here than were added into the network visualization (it's a scrolling list) - this is to keep things tidy in the visualization. By default, only the top 5 results (scored by retro 'score') are added to the visualization (this can be changed in the settings menu). The plus (+) and minus (-) buttons, when shown below the reaction, can be used to add or remove that reaction to or from the visualization. Go ahead and give it a try if you'd like. Click 'Next' to continue.",
            placement: "left",
            onNext: function() {
                app.data.nodes.forEach(function(n) {
                    if (n.reactionSmiles === 'Fc1ccc(C2(Cn3cncn3)CO2)c(F)c1.c1nc[nH]n1>>OC(Cn1cncn1)(Cn2cncn2)c3ccc(F)cc3F') {
                        app.network.selectNodes([n.id])
                        app.selected = app.data.nodes.get(n.id);
                    }
                })
            }
        },
        {
            element: '#details',
            title: "Viewing reaction details",
            content: "If you have a reaction node selected, Rank number 1 in this example, the right side of your screen will show you details for that reaction. At the top you can see the reaction SMILES, a 2d rendering, and similar reaction scores that you have seen before. You will also see a list of links to templates that support the reaction. Clicking one will open a new tab with more details about each template. There is also a link to 'Evaluate reaction in new tab', which will let you predict reaction conditions and evaluate the reaction in the forward direction.",
            placement: "left",
            onNext: function() {
                app.data.nodes.forEach(function(n) {
                    if (n.smiles === 'Fc1ccc(C2(Cn3cncn3)CO2)c(F)c1') {
                        app.network.selectNodes([n.id])
                        app.selected = app.data.nodes.get(n.id);
                    }
                })
            }
        },
        /*  Start of cluster section */
        {
            element: '#details',
            title: "Result clustering, Group similar",
            content: "You may quickly notice as you scroll down through the results, some are not shown. For example, the second result in the list for the currently selected chemical is ranked #8. This is because the 'Group similar' checkbox is checked. With the option enabled, results perceived to be the same by an unsupervised machine learning clustering algorithm are collapsed into the same group. In this way, only 1 representative example for the top 5 groups are added to the visualization by default, making it easier to browse the meaningfully different transformations. Click on 'Next' to uncheck 'Group similar' to reveal the hidden results.",
            placement: "left" ,
            onNext: function() {
                app.allowCluster = false ; 
            } 
        },
        {
            element: '#details',
            title: "Result clustering sorting options",
            content: "Now that you have confirmed that all the results are there, please also notice the drop-down box that appears. This drop-down box allows you to re-order the results depending on the Score, the number of examples used in that prediction, the template score, plausibility, root mean square of the molecular weight and finally the number of rings. Please select a different scoring order and see how it changes the results. Click 'Next' to re-enable clustering and to continue.",
            placement: "left",
            onNext: function() {
                app.allowCluster = true ; 
            }
        },
       {
            element: '#details',
            title: "Viewing clusters",
            content: "You may want to view the clusters to see which predictions have been grouped together. Each cluster set is displayed in a box that shows the reactants, Rank, Score etc. You will also notice the Red or Green box and another button with 4 squares in it, this is the view cluster button. Click 'Next' to view the clusters.",
            placement: "left",
            onNext: function() {
		        app.openClusterPopoutModal(app.selected, app.results[app.selected.smiles][0]);
            }
        },
        {
            title: "Cluster UI",
            content: "Here you can see all of the different reactions that were grouped together in the same cluster. The green (+) or red (-) buttons can be used to choose a different variant of the reaction type and add it or remove it from the graph visualization. Click 'Next' to close the cluster UI popup and continue with the tutorial.",
            orphan: true,
            onNext: function() {
                app.closeClusterPopoutModal()
            }
        },
	/* End of cluster section */
        {
            element: "#network",
            title: "Understanding the network",
            content: "We can see that the prediction gave us a few reactions that use buyable starting materials (<b>reaction</b> nodes with children <b>chemical</b> nodes highlighted in green) to make the new target we were interested in. Now we have a full path to our original target, Fluconazole, starting with buyable compounds.",
            placement: "right"
        },
        {
            title: "What do these buttons do?",
            element: "#expand-btn",
            content: "Lets quickly discuss what the buttons in this row do.",
            placement: "left"
        },
        {
            title: "Expand button",
            element: "#expand-btn",
            content: "As you have seen before, the Expand button will perform a prediction on the selected node and display the results in the network visualisation.",
            placement: "bottom"
        },
        {
            element: "#delete-btn",
            title: "Delete button",
            content: "In addition to expanding nodes, you can easily delete selected nodes, or children of a selected node using this button.",
            placement: "bottom"
        },
        {
            element: "#collapse-btn",
            title: "Collapse children button",
            content: "The 'Collapse children' button will group the currently selected node and its children into one cluster (this may be useful to keep things organized). Clicking this button with a cluster selected will expand it to show all of the nodes again.",
            placement: "bottom"
        },
        {
            element: "#clear-reactions-btn",
            title: "Clear reactions button",
            content: "This button will reset your target and clear all the reactions from the visualization.",
            placement: "right"
        },
        {
            element: "#settings-btn",
            title: "Settings",
            content: "There are various advanced settings for one-step/tree builder prediction parameters, as well as graph visualization options. Use this button to open the settings UI, where you can read more about each by hovering your mouse over the the tooltip icon (i).",
            placement: "right"
        },
        {
            element: "#hierarchical-button",
            title: "Hierarchical/Graph button",
            content: "Clicking on this button changes how the results are displayed below. The default mode is graphical, G, where the target is displayed in the center and the child nodes fan out in all directions. Clicking on this button will change the display to hierarchical mode where the target appears at the top of the tree and the child node(s) project downwards. Click on the button to try it out.",
            placement: "right"
        },
        {
            element: "#center-graph-button",
            title: "Center graph button",
            content: "This button will fit the network visualization inside the canvas. This is useful if you have zoomed in on a specific region but would like to reset the view.",
            placement: "right"
        },
        {
            element: "#download-btn",
            title: "Saving results",
            content: "You can save the network structure (JSON of nodes and edges) and download it to your computer. You may also share these JSON files with your colleagues so that they can get excited about the molecules you are working on.",
            placement: "right"
        },
        {
            element: "#load-btn",
            title: "Restoring results",
            content: "You can restore a previously saved network here.",
            placement: "right"
        },
        {
            element: "#tb-submit",
            title: "Tree builder button",
            content: "You can start a tree builder job, using the target SMILES string, by clicking on this button. This is an asynchronous job, so you can continue examining the predictions in the main window below. Once the job has completed, a browser popup will appear informing you that the job has finished. Clicking on that popup will bring you to the tree builder visualization page.",
            placement: "left"
        },
        {
            element: "#tb-submit-settings",
            title: "Tree builder button",
            content: 'You can give your tree builder job a name using this dropdown menu, as well as choose between a few different preset "quick" settings. The tree builder job that gets submitted will show up in your saved results accessible from the "My Results" link all the way at the top of the page, and a name can help identify this job later. If not provided, it will default to the SMILES string of your target. If you want more control over the tree builder parameters, you can go into the advanced settings and look for the "MCTS Tree Builder Settings" section.',
            placement: "right"
        },
        {
            title: "End of tour",
            content: "That's the end of the guided tour. Go ahead and change the target, and build your own reaction networks, or continue to expand this one further.",
            orphan: true
        }
    ]
});

function closeAll() {
    app.showSettingsModal = false;
    app.showLoadModal = false;
    app.showDownloadModal = false;
    app.showClusterPopoutModal = false;
    app.showClusterEditModal = false;
    app.showAddNewPrecursorModal = false;
}

/* key binding */
var keys = keycharm();
keys.bind("esc", closeAll, 'keyup');
