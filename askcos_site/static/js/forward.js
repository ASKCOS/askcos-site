function showLoader() {
    var loader = document.getElementById("pageLoader");
    loader.style.display = "block";
}

function hideLoader() {
    var loader = document.getElementById("pageLoader");
    loader.style.display = "none";
}

var app = new Vue({
    el: '#app',
    data: {
        reactants: '',
        product: '',
        reagents: '',
        solvent: '',
        evaluating: false,
        showSettings: false,
        numForwardResults: 100,
        numContextResults: 10,
        forwardResults: [],
        contextResults: [],
        impurityResults: [],
        reactionScore: null,
        mode: 'forward',
        forwardModel: 'wln',
        inspectionModel: 'fastFilter',
        atomMappingModel: 'wln',
        impurityTopk: 3,
        inspectionThreshold: 0.75,
        impurityProgress: {
            percent: 0,
            message: ''
        }
    },
    mounted: function() {
        var urlParams = new URLSearchParams(window.location.search)
        let mode = urlParams.get('mode')
        let reactants = urlParams.get('reactants')
        if (!!reactants) {
            this.reactants = reactants
        }
        let product = urlParams.get('product')
        if (!!product) {
            this.product = product
        }
        let reagents = urlParams.get('reagents')
        if (!!reagents) {
            this.reagents = reagents
        }
        let solvent = urlParams.get('solvent')
        if (!!solvent) {
            this.solvent = solvent
        }
        if (!!mode) {
            this.changeMode(mode)
            this.predict()
        }
    },
    methods: {
        clear() {
            this.reactants = ''
            this.product = ''
            this.reagents = ''
            this.solvent = ''
            this.forwardResults = []
            this.contextResults = []
            this.impurityResults = []
            this.impurityProgress = {
                percent: 0,
                message: ''
            }
        },
        changeMode(mode) {
            this.mode = mode
        },
        constructForwardQuery(reagents, solvent) {
            var query = `reactants=${encodeURIComponent(this.reactants)}`
            if (!!reagents) {
                query += `&reagents=${encodeURIComponent(reagents)}`
            }
            if (!!solvent) {
                query += `&solvent=${encodeURIComponent(solvent)}`
            }
            query += `&num_results=${this.numForwardResults}`
            return query
        },
        constructuContextQuery() {
            var reactants = encodeURIComponent(this.reactants)
            var product = encodeURIComponent(this.product)
            return `reactants=${reactants}&products=${product}&return_scores=True&num_results=${this.numContextResults}`
        },
        constructFastFilterQuery() {
            var reactants = encodeURIComponent(this.reactants)
            var product = encodeURIComponent(this.product)
            return `reactants=${reactants}&products=${product}`
        },
        constructImpurityQuery() {
            var query = `reactants=${encodeURIComponent(this.reactants)}`
            query += `&products=${encodeURIComponent(this.product)}`
            if (!!this.reagents) {
                query += `&reagents=${encodeURIComponent(this.reagents)}`
            }
            if (!!this.solvent) {
                query += `&solvent=${encodeURIComponent(this.solvent)}`
            }
            return query
        },
        predict() {
            switch(this.mode) {
                case 'forward':
                    this.forwardPredict()
                    break;
                case 'context':
                    this.contextPredict()
                    break;
                case 'evaluate':
                    this.evaluatePredict()
                    break;
                case 'impurity':
                    this.impurityPredict()
                    break;
                default:
                    alert('unsupported mode')
            }
        },
        forwardPredict() {
            showLoader()
            this.forwardResults = []
            var query = this.constructForwardQuery(this.reagents, this.solvent)
            fetch('/api/forward/?'+query)
            .then(resp => resp.json())
            .then(json => {
                this.forwardResults = json['outcomes']
                hideLoader()
            })
        },
        goToContext(smiles) {
            this.product = smiles
            this.mode = 'context'
            this.contextPredict()
        },
        goToForward(index) {
            var context = this.contextResults[index]
            var reagents = ''
            if (context['reagent']) {
                reagents += context['reagent']
            }
            if (context['catalyst']) {
                reagents += '.'+context['catalyst']
            }
            this.reagents = reagents
            if (context['solvent']) {
                this.solvent = context['solvent']
            }
            this.mode = 'forward'
            this.forwardPredict()
        },
        goToImpurity(index) {
            var context = this.contextResults[index]
            var reagents = ''
            if (context['reagent']) {
                reagents += context['reagent']
            }
            if (context['catalyst']) {
                reagents += '.'+context['catalyst']
            }
            this.reagents = reagents
            if (context['solvent']) {
                this.solvent = context['solvent']
            }
            this.mode = 'impurity'
            this.impurityPredict()
        },
        contextPredict() {
            showLoader()
            this.contextResults = []
            var query = this.constructuContextQuery()
            fetch('/api/context/?'+query)
            .then(resp => resp.json())
            .then(json => {
                this.contextResults = json['contexts']
                hideLoader()
            })
        },
        evaluateIndex(index) {
            var reagents = this.contextResults[index]['reagent']
            if (!!this.contextResults[index]['catalyst']) {
                if (!!reagents) {
                    reagents += '.'
                }
                reagents += this.contextResults[index]['catalyst']
            }
            var solvent = this.contextResults[index]['solvent']
            var query = this.constructForwardQuery(reagents, solvent)
            return fetch('/api/forward/?'+query)
            .then(resp => resp.json())
            .then(json => {
                for (var n in json['outcomes']) {
                    var outcome = json['outcomes'][n]
                    if (outcome['smiles'] == this.product) {
                        this.$set(this.contextResults[index], 'evaluation', Number(n)+1)
                        break
                    }
                }
                if (!this.contextResults[index]['evaluation']) {
                    this.$set(this.contextResults[index], 'evaluation', 0)
                }
            })
        },
        evaluate() {
            this.evaluating = true
            var query = this.constructFastFilterQuery()
            fetch('/api/fast-filter/?'+query)
            .then(resp => resp.json())
            .then(json => {
                this.reactionScore = json['score']
            })
            var promises = []
            for (index in this.contextResults) {
                promises.push(this.evaluateIndex(index))
            }
            Promise.all(promises).then(() => {this.evaluating = false})
        },
        updateImpurityProgress(taskId) {
            hideLoader()
            fetch(`/api/celery/task/?task_id=${taskId}`)
            .then(resp => resp.json())
            .then(json => {
                if (json['complete']) {
                    this.impurityProgress.percent = 1.0
                    this.impurityProgress.message = 'Prediction complete!'
                    this.impurityResults = json.results.predict_expand
                }
                else if (json['failed']) {
                    this.impurityProgress.percent = 0.0
                    this.impurityProgress.message = 'impurity prediction failed!'
                }
                else {
                    this.impurityProgress.percent = json['percent']
                    this.impurityProgress.message = json['message']
                    setTimeout(() => {this.updateImpurityProgress(taskId)}, 1000)
                }
            })
        },
        impurityPredict() {
            showLoader()
            this.impurityResults = []
            var query = this.constructImpurityQuery()
            fetch('/api/impurity/?'+query)
            .then(resp => resp.json())
            .then(json => {
                this.updateImpurityProgress(json['task_id'])
            })
        },
        startTour() {
            this.clear()
            tour.init()
            tour.restart()
        }
    },
    delimiters: ['%%', '%%'],
});

var tour = new Tour({
    storage: false,
    steps: [
        {
            title: "A guided tour through synthesis prediction",
            content: "Welcome to this guided tour",
            orphan: true,
            backdropContainer: '#body'
        }
    ]
});