function showLoader() {
    var loader = document.getElementsByClassName("loader")[0];
    loader.style.display = "block";
}

function hideLoader() {
    var loader = document.getElementsByClassName("loader")[0];
    loader.style.display = "none";
}

function getCookie(cname) {
    var name = cname + "=";
    var cookie_str = document.cookie;
    if (cookie_str && cookie_str != '') {
        var cookie_splitted = cookie_str.split(';');
        for (var i = 0; i < cookie_splitted.length; i++) {
            var c = cookie_splitted[i].trim();
            if (c.indexOf(name) == 0) {
                return decodeURIComponent(c.substring(name.length, c.length));
            }
        }
    }
    return undefined;
}

function copyToClipboard(text) {
    var dummy = document.createElement("textarea");
    document.body.appendChild(dummy);
    dummy.value = text;
    dummy.select();
    document.execCommand("copy");
    document.body.removeChild(dummy);
}

var app = new Vue({
    el: '#app',
    data: {
        templateId: null,
        templateInfo: {},
        templateReactions: []
    },
    mounted: function () {
        var urlParams = new URLSearchParams(window.location.search);
        this.templateId = urlParams.get('id')
        this.lookupTemplate()
    },
    methods: {
        lookupTemplate() {
            if (!!this.templateId) {
                fetch('/api/template/?id=' + this.templateId)
                    .then((resp) => resp.json())
                    .then(json => {
                        this.templateInfo = json.template
                        this.lookupReactions(json.template)
                    })
            }
        },
        lookupReactions() {
            var references = [...this.templateInfo.references]
            fetch(
                '/api/reactions/',
                {
                    method: 'POST',
                    headers: {
                        'Content-Type': 'application/json',
                        'X-CSRFToken': getCookie('csrftoken'),
                    },
                    body: JSON.stringify({
                        template_set: this.templateInfo.template_set,
                        ids: references.splice(0, 100)
                    })
                }
            )
            .then(resp => resp.json())
            .then(json => {
                this.templateReactions = json.reactions
            })
        },
        copyRxnIds() {
            var copyTooltip = document.querySelector('#copy-tooltip')
            copyToClipboard(this.reactionReferences)
            copyTooltip.innerHTML = 'Copied!'
            setTimeout(() => {copyTooltip.innerHTML = "Click to copy!"}, 2000)
        },
        downloadReactionQuery() {
            fetch(`/api/template/download/?id=${this.templateId}`)
            .then(resp => resp.json())
            .then(json => {
                var dataStr = "data:text/json;charset=utf-8," + JSON.stringify(json)
                var dlAnchorElem = document.getElementById('downloadAnchorElem')
                dlAnchorElem.setAttribute("href",     dataStr     )
                dlAnchorElem.setAttribute("download", "reaxys_query.json")
                dlAnchorElem.click()
            })
                
        },
    },
    computed: {
        reactionReferences() {
            if (!!this.templateInfo && !!this.templateInfo.references) {
                return this.templateInfo.references.join('; ')
            }
            else {
                return ''
            }
        },
        topReactionReferences() {
            if (!!this.templateInfo && !!this.templateInfo.references) {
                return this.templateInfo.references.splice(0, 100).join(';')
            }
            else {
                return ''
            }
        }
    },
    delimiters: ['%%', '%%'],
});
