var app = new Vue({
    el: '#app',
    data: {results: []},
    created: function() {
        this.update();
    },
    methods: {
        update: function() {
            showLoader();
            fetch('/api/my-results/')
            .then(resp => resp.json())
            .then(json => {
                console.log(json['results']);
                this.results = json['results'];
            })
            .finally(() => hideLoader())
        },
        removeResult: function(id) {
            if (window.confirm('Click "OK" to confirm deleting this result')) {
                showLoader();
                fetch('/api/remove-result/?id='+id)
                .then(resp => resp.json())
                .then(json => {
                    if (json['status']==1) {
                        for (i in this.results) {
                            if (this.results[i]['id']==id) {
                                this.results.splice(i, 1)
                            }
                        }
                    }
                })
                .finally(() => hideLoader())
            }
        }
    },
    delimiters: ['%%', '%%'],
});
