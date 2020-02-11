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
            var reagents = encodeURIComponent(this.reagents)
            var solvent = encodeURIComponent(this.solvent)
            return `reactants=${reactants}&reagents=${reagents}&solvent=${solvent}&num_results=${this.numForwardResults}`
        },
        constructuContextQuery() {
            var reactants = encodeURIComponent(this.reactants)
            var product = encodeURIComponent(this.product)
            return `reactants=${reactants}&products=${product}&return_scores=True&num_results=${this.numContextResults}`
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
        }
    },
    delimiters: ['%%', '%%'],
});
