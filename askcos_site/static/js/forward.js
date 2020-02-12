function showLoader() {
    var loader = document.getElementsByClassName("loader")[0];
    loader.style.display = "block";
}

function hideLoader() {
    var loader = document.getElementsByClassName("loader")[0];
    loader.style.display = "none";
}

var app = new Vue({
    el: '#app',
    data: {
        reactants: '',
        product: '',
        reagents: '',
        solvent: '',
        numForwardResults: 100,
        numContextResults: 10,
        forwardResults: [],
        contextResults: [],
        evaluationResults: [],
        mode: 'forward',
        forwardModel: 'wln',
        inspectionModel: 'fastFilter',
        atomMappingModel: 'wln',
        impurityTopk: 3,
        inspectionThreshold: 0.75
    },
    methods: {
        constructForwardQuery() {
            var reactants = encodeURIComponent(this.reactants)
            return `reactants=${reactants}&num_results=${this.numForwardResults}`
        },
        constructEvaluationQuery(reagents, solvent) {
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
            var query = this.constructForwardQuery()
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
        contextPredict() {
            showLoader()
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
            var query = this.constructEvaluationQuery(reagents, solvent)
            fetch('/api/forward/?'+query)
            .then(resp => resp.json())
            .then(json => {
                this.evaluationResults[index]['found'] = false
                for (outcome of json['outcomes']) {
                    if (outcome['smiles'] == this.product) {
                        this.evaluationResults[index]['found'] = true
                        break
                    }
                }
            })
        },
        evaluate() {
            var query = this.constructFastFilterQuery()
            fetch('/api/fast-filter/?'+query)
            .then(resp => resp.json())
            .then(json => {
                this.reactionScore = json['score']
            })
            for (index in this.contextResults) {
                this.evaluationResults.push({})
                this.evaluateIndex(index)
            }
        },
        showEvaluation(index) {
            return (!!this.evaluationResults) && (!!this.evaluationResults[index]) && (!!this.evaluationResults[index]['found'])
        }
    },
    delimiters: ['%%', '%%'],
});
