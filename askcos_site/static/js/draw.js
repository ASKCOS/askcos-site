function copyToClipboard(text) {
    const dummy = document.createElement("textarea");
    document.body.appendChild(dummy);
    dummy.value = text;
    dummy.select();
    document.execCommand("copy");
    document.body.removeChild(dummy);
}

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
        drawMap: false,
        highlight: false,
    },
    methods: {
        updateSmilesFromJSME() {
            this.canonicalize(jsmeApplet.smiles(), drawBoxId)
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
