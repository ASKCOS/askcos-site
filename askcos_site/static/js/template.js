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
                        ids: this.templateInfo.references
                    })
                }
            )
            .then(resp => resp.json())
            .then(json => {
                this.templateReactions = json.reactions
            })
        }
    },
    computed: {
        reactionReferences() {
            if (!!this.templateInfo && !!this.templateInfo.references) {
                return this.templateInfo.references.join('; ')
            }
            else {
                return ''
            }
        }
    },
    delimiters: ['%%', '%%'],
});
