function getCookie(cname) {
    var name = cname + "=";
    var cookie_str = document.cookie;
    if (cookie_str && cookie_str !== '') {
        var cookie_parts = cookie_str.split(';');
        for ( var i = 0; i <cookie_parts.length; i++ ) {
            var c = cookie_parts[i].trim();
            if (c.indexOf(name) === 0) {
                return decodeURIComponent(c.substring(name.length, c.length));
            }
        }
    }
    return undefined;
}

var app = new Vue({
    el: "#app",
    data: {
        smiles: '',
    },
    methods: {
        updateSmilesFromJSME() {
            var smiles = jsmeApplet.smiles();
            this.canonicalize(smiles, drawBoxId)
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
    },
    delimiters: ['%%', '%%'],
})
