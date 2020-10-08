solventScoreColorMap = {
    1: '#008000',
    2: '#3D9900',
    3: '#8FB300',
    4: '#CCA300',
    5: '#E65C00',
    6: '#FF0000',
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
        selectivityResults: [],
        reactionScore: null,
        mode: 'context',
        contextModel: 'neuralnetwork',
        forwardModel: 'wln',
        inspectionModel: 'fastFilter',
        atomMappingModel: 'wln',
        selectivityModel: 'qm_GNN',
        impurityTopk: 3,
        inspectionThreshold: 0.75,
        impurityCheckMapping: true,
        impurityProgress: {
            percent: 0,
            message: ''
        },
        // selectivity settings
        absoluteReagents: true,
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
        clearSelectivity() {
            this.selectivityResults = []
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
        constructForwardPostData(reagents, solvent) {
            var data = {
                reactants: this.reactants,
                num_results: this.numForwardResults
            }
            if (!!reagents) {
                data.reagents = reagents
            }
            if (!!solvent) {
                data.solvent = solvent
            }
            return data
        },
        constructContextPostData() {
            return {
                reactants: this.reactants,
                products: this.product,
                return_scores: true,
                num_results: this.numContextResults
            }
        },
        constructFastFilterPostData() {
            return {
                reactants: this.reactants,
                products: this.product
            }
        },
        constructSelectivityPostData() {
            var data= {
                reactants: this.reactants,
                product: this.product,
                mapper: this.atomMappingModel,
                no_map_reagents: this.absoluteReagents,
            }
            if (!!this.reagents) {
                data.reagents = this.reagents
            }
            if (!!this.solvent) {
                data.solvent = this.solvent
            }
            return data
        },
        constructImpurityPostData() {
            var data = {
                reactants: this.reactants,
                products: this.product,
                top_k: this.impurityTopk,
                threshold: this.inspectionThreshold,
                check_mapping: this.impurityCheckMapping
            }
            if (!!this.reagents) {
                data.reagents = this.reagents
            }
            if (!!this.solvent) {
                data.solvent = this.solvent
            }
            return data
        },
        apiAsyncPost(endpoint, postData, callback) {
            return fetch(endpoint, {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                    'X-CSRFToken': getCookie('csrftoken')
                },
                body: JSON.stringify(postData)
            })
                .then(resp => {
                    if (!resp.ok) {
                        throw Error(resp.statusText)
                    }
                    return resp
                })
                .then(resp => resp.json())
                .then(json => {
                    callback(json)
                })
                .catch(error => {
                    hideLoader();
                    alert('There was an error executing this POST request: '+error)
                })
        },
        celeryTaskAsyncPost(taskName, postData, callback) {
            var celeryCallback = (json) => {
                setTimeout(() => this.pollCeleryResult(json.task_id, callback), 1000)
            }
            return this.apiAsyncPost(`/api/v2/${taskName}/`, postData, celeryCallback)
        },
        pollCeleryResult: function(taskId, complete, progress, failed) {
            fetch(`/api/v2/celery/task/${taskId}/`)
            .then(resp => resp.json())
            .then(json => {
                if (json.complete) {
                    complete(json);
                    hideLoader();
                }
                else if (json.failed) {
                    if (!!failed) {
                        failed(json)
                    }
                    hideLoader();
                    throw Error('Celery task failed.');
                }
                else {
                    if (!!progress) {
                        progress(json)
                    }
                    setTimeout(() => {this.pollCeleryResult(taskId, complete, progress, failed)}, 1000)
                }
            })
            .catch(error => {
                if (error instanceof TypeError) {
                    console.log('Unable to fetch celery results due to connection error. Will keep trying.')
                    setTimeout(() => {this.pollCeleryResult(taskId, complete, progress, failed)}, 2000)
                } else {
                    console.error('There was a problem fetching results:', error);
                }
            });
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
                    case 'selectivity':
                        this.clearSelectivity()
                        this.selectivityPredict()
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
            var postData = this.constructForwardPostData(this.reagents, this.solvent)
            var callback = (json) => {
                this.forwardResults = json.output
            }
            this.celeryTaskAsyncPost('forward', postData, callback)
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
        goToSelectivity(smiles) {
            window.history.pushState({mode: 'selectivity'}, 'selectivity', '?mode=selectivity')
            this.canonicalizeAll()
                .then(() => {
                    this.product = smiles
                    this.mode = 'selectivity'
                    this.selectivityPredict()
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
            var postData = this.constructContextPostData()
            var callback = (json) => {
                this.contextResults = json.output
            }
            this.celeryTaskAsyncPost('context', postData, callback)
        },
        selectivityPredict() {
            showLoader()
            this.selectivityResults = []
            var postData = this.constructSelectivityPostData()
            var callback = (json) => {
                this.selectivityResults = json.output
            }
            this.celeryTaskAsyncPost('general-selectivity', postData, callback)
        },
        evaluateIndex(index) {
            this.$set(this.contextResults[index], 'evaluating', true)
            var reagents = this.contextResults[index]['reagent']
            if (!!this.contextResults[index]['catalyst']) {
                if (!!reagents) {
                    reagents += '.'
                }
                reagents += this.contextResults[index]['catalyst']
            }
            var solvent = this.contextResults[index]['solvent']

            var postData = this.constructForwardPostData(reagents, solvent)
            const app = this;
            var callback = (json) => {
                var outcomes = json.output
                for (var n in outcomes) {
                    var outcome = outcomes[n]
                    if (outcome.smiles == app.product) {
                        app.$set(app.contextResults[index], 'evaluation', Number(n)+1)
                        break
                    }
                }
                if (!app.contextResults[index]['evaluation']) {
                    app.$set(app.contextResults[index], 'evaluation', 0)
                }
                app.$set(app.contextResults[index], 'evaluating', false)
            }
            this.celeryTaskAsyncPost('forward', postData, callback)
        },
        canonicalize(smiles, input) {
            var postData = {smiles: smiles}
            var callback = (json) => {
                this[input] = json.smiles
            }
            return this.apiAsyncPost('/api/v2/rdkit/smiles/canonicalize/', postData, callback)
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
            var postData = this.constructFastFilterPostData()
            var callback = (json) => {
                this.reactionScore = json.output
            }
            this.celeryTaskAsyncPost('fast-filter', postData, callback)
            for (index in this.contextResults) {
                this.evaluateIndex(index)
            }
        },
        updateImpurityProgress(taskId) {
            hideLoader()
            var complete = (json) => {
                this.impurityProgress.percent = 1.0
                this.impurityProgress.message = 'Prediction complete!'
                this.impurityResults = json.output.predict_expand
            }
            var progress = (json) => {
                this.impurityProgress.percent = json['percent']
                this.impurityProgress.message = json['message']
            }
            var failed = (json) => {
                this.impurityProgress.percent = 0.0
                this.impurityProgress.message = 'impurity prediction failed!'
            }
            this.pollCeleryResult(taskId, complete, progress, failed)
        },
        impurityPredict() {
            showLoader()
            this.impurityResults = []
            var postData = this.constructImpurityPostData()
            var callback = (json) => {
                this.updateImpurityProgress(json.task_id)
            }
            this.apiAsyncPost('/api/v2/impurity/', postData, callback)
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
        downloadSelectivityResults() {
            if (!!!this.selectivityResults) {
                alert('There are no regio-selectivity results to download!')
            }
            var downloadData = "Rank,SMILES,Probability\n"
            this.selectivityResults.forEach((res) => {
                downloadData += `${res.rank},${res.smiles},${res.prob}\n`
            })
            var dataStr = "data:text/json;charset=utf-8," + downloadData
            var dlAnchorElem = document.getElementById('downloadSelectivityAnchorElem')
            dlAnchorElem.setAttribute("href",     dataStr     )
            dlAnchorElem.setAttribute("download", "askcos_regioselectivity_export.csv")
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
