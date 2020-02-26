function showLoader() {
    var loader = document.getElementById("pageLoader");
    loader.style.display = "block";
}

function hideLoader() {
    var loader = document.getElementById("pageLoader");
    loader.style.display = "none";
}

function getCookie(cname) {
    var name = cname + "=";
    var cookie_str = document.cookie;
    if (cookie_str && cookie_str != '') {
        var cookie_splitted = cookie_str.split(';');
        for(var i = 0; i <cookie_splitted.length; i++) {
            var c = cookie_splitted[i].trim();
            if (c.indexOf(name) == 0) {
                return decodeURIComponent(c.substring(name.length, c.length));
            }
        }
    }
  return undefined;
}

Vue.component('modal', {
    template: '#modal-template'
})

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
        mode: 'context',
        forwardModel: 'wln',
        inspectionModel: 'fastFilter',
        atomMappingModel: 'wln',
        impurityTopk: 3,
        inspectionThreshold: 0.75,
        impurityCheckMapping: true,
        impurityProgress: {
            percent: 0,
            message: ''
        }
    },
    mounted: function() {
        var urlParams = new URLSearchParams(window.location.search)
        let mode = urlParams.get('mode')

        let rxnsmiles = urlParams.get('rxnsmiles')
        if (!!rxnsmiles) {
            var split = rxnsmiles.split('>>')
            this.reactants = split[0]
            this.product = split[split.length-1]
        }
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
        clearContext() {
            this.contextResults = []
            this.reactionScore = null
        },
        clearEvaluation() {
            this.reactionScore = null
            for (var res of this.contextResults) {
                res.evaluation = undefined
            }
        },
        clearForward() {
            this.forwardResults = []
        },
        clearImpurity() {
            this.impurityResults = []
            this.impurityProgress = {
                percent: 0,
                message: ''
            }
        },
        clearSmiles() {
            this.reactants = ''
            this.product = ''
            this.reagents = ''
            this.solvent = ''
        },
        clear() {
            this.clearSmiles()
            this.clearContext()
            this.clearForward()
            this.clearImpurity()
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
            query += `&top_k=${this.impurityTopk}&threshold=${this.inspectionThreshold}&check_mapping=${this.impurityCheckMapping}`
            return query
        },
        predict() {
            switch(this.mode) {
                case 'forward':
                    this.clearForward()
                    this.forwardPredict()
                    break;
                case 'context':
                    this.clearContext()
                    this.contextPredict()
                    break;
                case 'impurity':
                    this.clearImpurity()
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
        goToImpurity(smiles) {
            this.product = smiles
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
            if (this.evaluating) {
                return
            }
            this.clearEvaluation()
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
        updateSmilesFromJSME() {
            var smiles = jsmeApplet.smiles();
            fetch(
                '/api/rdkit/canonicalize/',
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
                this[drawBoxId] = json.smiles
            })
        },
        startTour() {
            var res = confirm('Starting the tutorial will clear all current results, continue?')
            if (res) {
                this.clear()
                this.mode = 'context'
                tour.init()
                tour.restart()
            }
        }
    },
    delimiters: ['%%', '%%'],
});

var tour = new Tour({
    storage: false,
    steps: [
        {
            title: "A guided tour through synthesis prediction",
            content: `
Welcome to this guided tour through synthesis prediction in ASKCOS, which will demonstrate how to use the different parts of the user interface (UI). 
Thanks to <a href='http://bootstraptour.com/' target='_blank'>bootstrap-tour</a> for the JavaScript package making it very easy to provide this tour to you!
`,
            orphan: true,
            backdropContainer: '#body'
        },
        {
            element: "#contextArrowStep",
            placement: "bottom",
            title: "Entrypoint to a synthesis prediction",
            content: `
The entrypoint is predicting possible reaction conditions given known reactants and products, for example, after making a retrosynthetic prediction for a given target.
`,

        },
        {
            title: "Providing molecular structures",
            element: "#reactants",
            placement: "left",
            content: `
Ultimately, the software will need SMILES strings for each compound. 
These can be entered directly (i.e. - copy and pasted from ChemDraw) or drawn using the JSME molecular editor we have provided in the UI.
When you start entering a SMILES string in each input field, the structure will be rendered, on-the-fly, as you type. 
Don't be alarmed if you are in the middle of writing a ring structure and the image looks broken. 
Give it a try by writing a simple molecule like "CCOCC" in the reactants input field to the right.
`
        },
        {
            title: "Drawing structures",
            element: "#reactants-edit-icon",
            placement: "bottom",
            content: `
Alternatively, structures can be drawn using a simple molecular editor by clicking this edit (pencil) button. 
Each input field has it's own button you should click to draw a structure for that field.
If a SMILES string is already present in the input field, the drawing interface will be prepopulated with that structure for you to edit.
Give it a try now if you'd like.
`,
            onNext: () => {
                app.reactants = "Brc1ccccc1.OB(O)c1ccccc1"
                app.product = "c1ccc(-c2ccccc2)cc1"
            }
        },
        {
            title: "An example prediction",
            element: "#submit-button",
            placement: "top",
            content: `
For this tutorial, let's take a look at an example suzuki coupling reaction. 
Reactants and products have been prepopulated for you, and you can run the prediction by clicking submit (or click continue and we'll pretend you clicked submit)
`,
            relfex: true,
            onNext: () => {
                if (app.contextResults.length == 0) {
                    app.predict()
                }
            }
        },
        {
            title: "Condition results",
            element: "#context-results",
            placement: "top",
            content: `
After the prediction on the server has finished, the top 10 results will be displayed below.
`,
            onNext: () => {
                app.evaluate()
            }
        },
        {
            title: "Evaluating reactions and conditions",
            element: "#context-results",
            placement: "top",
            content: `
You can evaluate the reaction and the reaction conditions using this evaluation button here (we just clicked it for you). 
For evaluation, each set of reaction conditions is sent through a forward reaction prediction. 
If the product is found in the top 10 forward prediction results, a green checkmark will appear next to the recommendation (with the rank assigned to the product).
Additionally, the reaction evluator (which does not currently consider reaction conditions) will give a reaction score for the transformation from reactants to products.
`
        },
        {
            title: "Moving to the forward prediction",
            element: "#predict-conditions-0",
            placement: "top",
            content: `
When you've picked a set of conditions you'd like to explore further, you can click this button to make a forward prediction and see all of the possible results. 
Either click this button or click 'Next >>' to continue.
`,
            reflex: true,
            onNext: () => {
                app.goToForward(0)
            }
        },
        {
            title: "Viewing results",
            orphan: true,
            backdropContainer: '#body',
            content: `
All of the results from the prediction after including the reaction conditions will be displayed here. 
This is always a good check to see if the product you entered appears at the top of this list, and with what probability. 
`
        },
        {
            title: "Predicting impurities",
            element: "#predict-impurities-0",
            placement: "top",
            content: `
When you've picked a major product from the list you can move on to predict impurities that may come from a variety of imurity modes (i.e. - over reaction, dimerization). 
Click this button or click "Next >>" to continue and make the impurity prediction.
`,
            reflex: true,
            onNext: () => {
                app.goToImpurity(app.product)
            }
        },
        {
            title: "Impurity prediction",
            orphan: true,
            backdropContainer: '#body',
            content: `
The impurity prediction will start on the server, and the progress will be displayed here. 
Once the task finishes, the results wil be shown here. That ends this tour of the interface!
`
        },
    ]
});