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
        contextModel: 'neuralnetwork',
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
        }
        if (!!this.reactants) {
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
            window.history.pushState({mode: mode}, mode, '?mode='+mode)
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
            showLoader()
            this.canonicalizeAll()
            .then(() => {
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
            })
        },
        forwardPredict() {
            showLoader()
            this.forwardResults = []
            if (this.reactants.length < 4) {
                alert('Please enter a reactant with at least 4 atoms.')
                hideLoader()
                return
            }
            var query = this.constructForwardQuery(this.reagents, this.solvent)
            fetch('/api/forward/?'+query)
            .then(resp => resp.json())
            .then(json => {
                this.forwardResults = json['outcomes']
                hideLoader()
            })
        },
        goToForward(index) {
            window.history.pushState({mode: 'forward'}, 'forward', '?mode=forward')
            this.canonicalizeAll()
            .then(() => {
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
            })
        },
        goToImpurity(smiles) {
            window.history.pushState({mode: 'impurity'}, 'impurity', '?mode=impurity')
            this.canonicalizeAll()
            .then(() => {
                this.product = smiles
                this.mode = 'impurity'
                this.impurityPredict()
            })
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
        canonicalize(smiles, input) {
            return fetch(
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
                this[input] = json.smiles
            })
        },
        canonicalizeAll() {
            var promises = []
            for (var smi of ['reactants', 'product', 'reagents', 'solvent']) {
                if (!!this[smi]) {
                    promises.push(
                        this.canonicalize(this[smi], smi)
                    )
                }
            }
            return Promise.all(promises)
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
            this.canonicalize(smiles, drawBoxId)
        },
        downloadForwardResults() {
            if (!!!this.forwardResults) {
                alert('There are no forward predictor results to download!')
            }
            var downloadData = "Rank,SMILES,Probability,Score,MolWt\n"
            this.forwardResults.forEach((res) => {
                downloadData += `${res.rank},${res.smiles},${res.prob},${res.score},${res.mol_wt}\n`
            })
            var dataStr = "data:text/json;charset=utf-8," + downloadData
            var dlAnchorElem = document.getElementById('downloadForwardAnchorElem')
            dlAnchorElem.setAttribute("href",     dataStr     )
            dlAnchorElem.setAttribute("download", "askcos_forward_export.csv")
            dlAnchorElem.click()
        },
        downloadImpurityResults() {
            if (!!!this.impurityResults) {
                alert('There are no impurity predictor results to download!')
            }
            var downloadData = "No.,SMILES,Mechanism,InspectorScore,SimilarityScore\n"
            this.impurityResults.forEach((res) => {
                downloadData += `${res.no},${res.prd_smiles},${res.modes_name},${res.avg_insp_score},${res.similarity_to_major}\n`
            })
            var dataStr = "data:text/json;charset=utf-8," + downloadData
            var dlAnchorElem = document.getElementById('downloadImpurityAnchorElem')
            dlAnchorElem.setAttribute("href",     dataStr     )
            dlAnchorElem.setAttribute("download", "askcos_impurity_export.csv")
            dlAnchorElem.click()
        },
        startTour() {
            var res = confirm('Starting the tutorial will clear all current results, continue?')
            if (res) {
                this.clear()
                this.mode = 'context'
                tour.restart()
            }
        }
    },
    delimiters: ['%%', '%%'],
});

var tour = new Tour({
    framework: 'bootstrap4',
    storage: false,
    steps: [
        {
            title: "A guided tour through synthesis prediction",
            content: `
Welcome to this guided tour through synthesis prediction in ASKCOS, which will demonstrate how to use the different parts of the user interface (UI). 
Thanks to <a href='https://github.com/IGreatlyDislikeJavascript/bootstrap-tourist' target='_blank'>bootstrap-tourist</a> for the JavaScript package making it very easy to provide this tour to you!
`,
            orphan: true,
            backdropContainer: '#body'
        },
        {
            element: "#contextArrowStep",
            placement: "bottom",
            title: "Entrypoint to a synthesis prediction",
            content: `
The entrypoint is predicting possible reaction conditions given known reactants and products. For example, after making a retrosynthetic prediction for a given target (look for a link saying "evaluate reaction in a new tab").
`,

        },
        {
            title: "Providing molecular structures",
            element: "#reactants",
            placement: "left",
            content: `
The predictor needs SMILES strings as inputs in the appropriate boxes. 
SMILES strings can be entered directly (i.e. - copy and paste from another source) or deduced from a drawing made in the provided molecular editor (more in the next tour popup).
When you start entering a SMILES string in each input field, the structure will be rendered, dynamically, as you type. 
Don't be alarmed if you are in the middle of writing a ring structure and the image looks broken. 
Give it a try by writing a simple molecule like "CCOCC" in the reactants input field to the right.
`
        },
        {
            title: "Drawing structures",
            element: "#reactants-edit-icon",
            placement: "bottom",
            content: `
Structures can also be drawn using the simple molecular editor, accessible through the edit (pencil) button.
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
Reactants and products have been prepopulated for you, and you can run the prediction by clicking submit (or click next and we'll pretend you clicked submit).
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
            orphan: true,
            backdropContainer: '#body',
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
If the product is found in the top 10 forward prediction results, a checkmark will appear next to the recommendation (with the rank assigned to the product).
Additionally, the reaction evaluator (which does not currently consider reaction conditions) will give a reaction score for the transformation from reactants to products.
`
        },
        {
            title: "Moving to the forward prediction",
            element: "#predict-conditions-0",
            placement: "top",
            content: `
When you've picked a set of conditions you'd like to explore further, you can click this button to make a forward prediction and see all of the possible results. 
Either click this button or click "Next >>" to continue.
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
When a major product has been selected, a list of possible impurities can be predicted by clicking on the arrow (or the "Next >>") button below. 
These impurities may come from a variety of impurity modes (i.e. - over reaction, dimerization).
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
The impurity prediction will start on the server, and the progress will be updated in the progress bar below. 
You may notice that several predictions are made based on various modes (i.e. - over reaction, dimerization). 
When the task finishes, the results will be shown below the progress bar. 
The results can be exported by clicking on the Export results button.
`
        },
        {
            title: "The end!",
            orphan: true,
            backdropContainer: '#body',
            content: `
That ends this tour of the interface!
`
        },
    ]
});
