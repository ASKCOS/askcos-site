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
        reagents: '',
        solvent: '',
        numResults: 100,
        results: []
    },
    methods: {
        constructQuery() {
            var reactants = encodeURIComponent(this.reactants)
            var reagents = encodeURIComponent(this.reagents)
            var solvent = encodeURIComponent(this.solvent)
            var query = [reactants, reagents, solvent].join('.')
            return `reactants=${reactants}&reagents=${reagents}&solvent=${solvent}&num_results=${this.numResults}`
        },
        predict() {
            showLoader()
            var query = this.constructQuery()
            fetch('/api/forward/?'+query)
            .then(resp => resp.json())
            .then(json => {
                this.results = json['outcomes']
                hideLoader()
            })
        }
    },
    delimiters: ['%%', '%%'],
});
