var app = new Vue({
    el: "#app",
    data: {
        smiles: '',
        drawMap: false,
        highlight: false,
    },
    methods: {
        updateSmilesFromKetcher() {
            let smiles = ketcher.getSmiles();
            this.canonicalize(smiles, drawBoxId)
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
    },
    delimiters: ['%%', '%%'],
})
