function num2str(n, len) {
    if (len === undefined) {
        return n === undefined || isNaN(n) ? 'N/A' : n.toString();
    } else {
        return n === undefined || isNaN(n) ? 'N/A' : n.toFixed(len);
    }
}

var app = new Vue({
    el: "#app",
    data: {
        smiles: '',
        validSmiles: false,
        scscore: undefined,
        reactionScore: undefined,
        mappedSmiles: undefined,
        mapper: 'WLN atom mapper',
        mapperOptions: ['WLN atom mapper', 'Transformer'],
        showMappedSmiles: true,
        tbStatus: undefined,
    },
    methods: {
        updateSmilesFromKetcher() {
            // Not canonicalized
            this[drawBoxId] = ketcher.getSmiles()
        },
        canonicalize(smiles, field) {
            fetch(
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
            .then(resp => {
                if (!resp.ok) {
                    throw Error(resp.statusText)
                }
                return resp
            })
            .then(resp => resp.json())
            .then(json => {
                this[field] = json.smiles
            })
            .catch(error => {
                console.log('Could not canonicalize: '+error)
            })
        },
        getScscore(smiles) {
            this.scscore = 'evaluating'
            const url = '/api/v2/scscore/'
            const body = {
                smiles: smiles,
            }
            fetch(url, {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                    'X-CSRFToken': getCookie('csrftoken'),
                },
                body: JSON.stringify(body)
            })
                .then(resp => {
                    if (!resp.ok) {
                        throw resp.statusText
                    }
                    return resp
                })
                .then(resp => resp.json())
                .then(json => {
                    this.scscore = num2str(json.score, 3)
                })
                .catch(error => {
                    console.error('Could not evaluate SCScore:', error)
                })
        },
        getReactionScore(smiles) {
            this.reactionScore = 'evaluating'
            const url = '/api/v2/fast-filter/'
            const split = smiles.split('>')
            const body = {
                reactants: split[0],
                products: split[split.length-1],
            }
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
                        throw resp.statusText
                    }
                    return resp
                })
                .then(resp => resp.json())
                .then(json => {
                    let callback = (res) => {this.reactionScore = num2str(res, 4)}
                    setTimeout(() => this.pollCeleryResult(json.task_id, callback), 1000)
                })
                .catch(error => {
                    console.error('Could not evaluate reaction score:', error)
                })
        },
        getMappedSmiles(smiles) {
            this.mappedSmiles = 'evaluating'
            this.showMappedSmiles = true
            const url = '/api/v2/atom-mapper/'
            const body = {
                rxnsmiles: smiles,
                mapper: this.mapper,
            }
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
                        throw resp.statusText
                    }
                    return resp
                })
                .then(resp => resp.json())
                .then(json => {
                    let callback = (res) => {this.mappedSmiles = res}
                    setTimeout(() => this.pollCeleryResult(json.task_id, callback), 1000)
                })
                .catch(error => {
                    console.error('Could not generate atom mapping:', error)
                })
        },
        sendTreeBuilderJob(smiles) {
            this.tbStatus = 'pending'
            const url = '/api/v2/tree-builder/'
            const body = {
                description: smiles,
                smiles: smiles,
                template_set: 'reaxys',
                template_prioritizer_version: 1,
                max_depth: 5,
                max_branching: 20,
                expansion_time: 60,
                num_templates: 1000,
                max_cum_prob: 0.999,
                filter_threshold: 0.1,
                max_trees: 500,
                store_results: true,
                return_first: false,
                json_format: 'nodelink',
            }
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
                        throw resp.statusText
                    }
                    return resp
                })
                .then(resp => resp.json())
                .then(json => {
                    this.tbStatus = json.task_id
                })
                .catch(error => {
                    this.tbStatus = 'error'
                    console.error('Failed to submit tree builder job:', error)
                })
        },
        pollCeleryResult(taskId, callback) {
            fetch(`/api/v2/celery/task/${taskId}/`)
            .then(resp => resp.json())
            .then(json => {
                if (json.complete) {
                    callback(json.output);
                }
                else if (json.failed) {
                    throw 'Celery task failed'
                }
                else {
                    setTimeout(() => {this.pollCeleryResult(taskId, callback)}, 1000)
                }
            })
            .catch(error => {
                if (error instanceof TypeError) {
                    console.log('Unable to fetch celery results due to connection error. Will keep trying.')
                    setTimeout(() => {this.pollCeleryResult(taskId, callback)}, 2000)
                } else {
                    console.error('There was a problem fetching results:', error)
                }
            });
        },
    },
    watch: {
        smiles: function(newVal, oldVal) {
            if (newVal !== oldVal) {
                this.validSmiles = false
                this.scscore = undefined
                this.reactionScore = undefined
                this.mappedSmiles = undefined
                this.tbStatus = undefined
            }
        },
    },
    computed: {
        type: function() {
            if (this.smiles.includes('>')) {
                return 'rxn'
            } else {
                return 'mol'
            }
        },
    },
    delimiters: ['%%', '%%'],
})
